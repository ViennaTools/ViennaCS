import math
import sys

import viennacs2d as vcs
import viennals2d as vls


def read_config_file(file_name: str) -> dict:
    params = {}

    with open(file_name, "r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue

            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            try:
                params[key] = float(value)
            except ValueError:
                params[key] = value

    return params


def make_bounds(x_extent: float, top_space: float, implant_depth: float, mask_height: float):
    return [-0.5 * x_extent, 0.5 * x_extent, -implant_depth, top_space + mask_height]


def make_domain(bounds, grid_delta: float):
    boundary_conditions = [
        vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
        vls.BoundaryConditionEnum.INFINITE_BOUNDARY,
    ]
    return vls.Domain(bounds, boundary_conditions, grid_delta)


def make_plane(domain, y_position: float):
    vls.MakeGeometry(domain, vls.Plane([0.0, y_position], [0.0, 1.0])).apply()


def make_box(domain, min_x: float, min_y: float, max_x: float, max_y: float):
    vls.MakeGeometry(domain, vls.Box([min_x, min_y], [max_x, max_y])).apply()


def add_level_set(level_sets, material_map, level_set, material_id: int, wrap_lower_level_set: bool = True):
    if level_sets and wrap_lower_level_set:
        vls.BooleanOperation(
            level_set, level_sets[-1], vls.BooleanOperationEnum.UNION
        ).apply()

    level_sets.append(level_set)
    material_map.insertNextMaterial(material_id)


def make_structure(
    x_extent: float,
    top_space: float,
    implant_depth: float,
    opening_width: float,
    mask_height: float,
    grid_delta: float,
):
    bounds = make_bounds(x_extent, top_space, implant_depth, mask_height)
    material_map = vls.MaterialMap()
    level_sets = []

    substrate_bottom = make_domain(bounds, grid_delta)
    make_plane(substrate_bottom, -implant_depth)
    add_level_set(level_sets, material_map, substrate_bottom, 1)

    substrate_top = make_domain(bounds, grid_delta)
    make_plane(substrate_top, 0.0)
    add_level_set(level_sets, material_map, substrate_top, 1)

    mask = make_domain(bounds, grid_delta)
    make_plane(mask, mask_height)

    opening = make_domain(bounds, grid_delta)
    make_box(
        opening,
        -0.5 * opening_width,
        -grid_delta,
        0.5 * opening_width,
        mask_height + grid_delta,
    )
    vls.BooleanOperation(mask, opening, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT).apply()
    add_level_set(level_sets, material_map, mask, 2)

    return level_sets, material_map


def build_cell_set(structure, top_space: float):
    cell_set = vcs.DenseCellSet()
    cell_set.setCellSetPosition(True)
    cell_set.setCoverMaterial(0)
    level_sets, material_map = structure
    cell_set.fromLevelSets(level_sets, material_map, top_space)
    cell_set.addScalarData("concentration", 0.0)
    cell_set.addScalarData("beamHits", 0.0)
    return cell_set


def frange(start: float, stop: float, step: float):
    value = start
    eps = 0.5 * step
    while value <= stop + eps:
        yield value
        value += step


def integrate_depth_profile(model, max_depth: float, step: float) -> float:
    area = 0.0
    depth = 0.0
    while depth <= max_depth:
        area += max(0.0, model.getDepthProfile(depth)) * step
        depth += step
    return area


def smoothstep(edge0: float, edge1: float, x: float) -> float:
    if edge1 <= edge0:
        return 1.0 if x >= edge0 else 0.0
    t = min(max((x - edge0) / (edge1 - edge0), 0.0), 1.0)
    return t * t * (3.0 - 2.0 * t)


def exponential_tail(
    depth: float, start_depth: float, decay_length: float, blend_width: float
) -> float:
    if decay_length <= 0.0:
        return 0.0

    onset = smoothstep(start_depth - 0.5 * blend_width, start_depth + 0.5 * blend_width, depth)
    if onset <= 0.0:
        return 0.0

    return onset * math.exp(-(depth - start_depth) / decay_length)


def integrate_exponential_tail(
    max_depth: float,
    step: float,
    start_depth: float,
    decay_length: float,
    blend_width: float,
) -> float:
    area = 0.0
    depth = 0.0
    while depth <= max_depth:
        area += exponential_tail(depth, start_depth, decay_length, blend_width) * step
        depth += step
    return area


def combined_depth_density(
    model,
    depth: float,
    random_norm: float,
    substrate_type: str,
    tail_fraction: float,
    tail_norm: float,
    tail_start_depth: float,
    tail_decay_length: float,
    tail_blend_width: float,
) -> float:
    random_density = max(0.0, model.getDepthProfile(depth)) / random_norm

    if substrate_type != "crystalline" or tail_fraction <= 0.0 or tail_norm <= 0.0:
        return random_density

    tail_density = (
        exponential_tail(depth, tail_start_depth, tail_decay_length, tail_blend_width)
        / tail_norm
    )
    return (1.0 - tail_fraction) * random_density + tail_fraction * tail_density


def apply_masked_implant(
    cell_set,
    model,
    opening_width: float,
    mask_height: float,
    angle_deg: float,
    beam_spacing: float,
    implant_depth: float,
    dose_cm2: float,
    substrate_type: str,
    tail_fraction: float,
    tail_start_depth: float,
    tail_decay_length: float,
    tail_blend_width: float,
):
    concentration = [0.0] * cell_set.getNumberOfCells()
    beam_hits = [0.0] * cell_set.getNumberOfCells()
    materials = cell_set.getScalarData("Material")
    angle = math.radians(angle_deg)

    # The Pearson IV implementation provides the profile shape, but not a
    # normalized probability density. Normalize it numerically over the depth
    # window of interest, then scale by an areal implant dose to obtain a
    # concentration field in cm^-3.
    depth_step = max(beam_spacing, 0.1)
    random_norm = integrate_depth_profile(model, implant_depth, depth_step)
    if random_norm <= 0.0:
        raise RuntimeError("Depth profile normalization failed.")

    tail_fraction = min(max(tail_fraction, 0.0), 1.0)
    tail_norm = integrate_exponential_tail(
        implant_depth, depth_step, tail_start_depth, tail_decay_length, tail_blend_width
    )

    dose_nm2 = dose_cm2 / 1.0e14
    dose_per_beam = dose_nm2 * beam_spacing

    # Only rays that pass through the aperture are launched.
    for x_mask in frange(-0.5 * opening_width, 0.5 * opening_width, beam_spacing):
        x_surface = x_mask - math.tan(angle) * mask_height

        for idx in range(cell_set.getNumberOfCells()):
            x_pos, y_pos, _ = cell_set.getCellCenter(idx)
            depth_axis = -y_pos

            if materials[idx] != 1.0 or depth_axis < 0.0:
                continue

            depth = math.cos(angle) * depth_axis + math.sin(angle) * (x_surface - x_pos)
            if depth < 0.0:
                continue

            lateral_displacement = abs(
                math.cos(angle) * (x_surface - x_pos) - math.sin(angle) * depth_axis
            )

            depth_density = combined_depth_density(
                model,
                depth,
                random_norm,
                substrate_type,
                tail_fraction,
                tail_norm,
                tail_start_depth,
                tail_decay_length,
                tail_blend_width,
            )
            lateral_density = model.getLateralProfile(lateral_displacement, depth)

            # Convert from atoms / nm^3 to atoms / cm^3 for a more useful field.
            concentration[idx] += dose_per_beam * depth_density * lateral_density * 1.0e21

            # Store a simple aperture-throughput field to make the mask footprint visible.
            if lateral_displacement <= 0.5 * beam_spacing:
                beam_hits[idx] += 1.0

    cell_set.setScalarData("concentration", concentration)
    cell_set.setScalarData("beamHits", beam_hits)


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <config file>")
        return 1

    params = read_config_file(sys.argv[1])

    grid_delta = params.get("gridDelta", 1.0)
    x_extent = params.get("xExtent", 40.0)
    top_space = params.get("topSpace", 6.0)
    implant_depth = params.get("implantDepth", 12.0)
    opening_width = params.get("openingWidth", 8.0)
    mask_height = params.get("maskHeight", 4.0)
    angle = params.get("angle", 7.0)
    beam_spacing = params.get("beamSpacing", grid_delta)
    implant_dose_cm2 = params.get("doseCm2", 1.0e15)
    substrate_type = str(params.get("substrateType", "amorphous")).strip().lower()

    # Representative low-energy B-in-Si moments for a shallow 10 keV implant
    # into Si(100). Rp is anchored to the commonly cited ~31 nm projected
    # range for B in Si at 10 keV, while skewness/kurtosis are chosen as a
    # realistic positively skewed Pearson IV profile for the random component.
    # Lateral spread is kept Gaussian as defined by ImplantPearsonIV.
    projected_range = params.get("projectedRange", 31.0)
    depth_sigma = params.get("depthSigma", 11.0)
    skewness = params.get("skewness", 1.15)
    kurtosis = params.get("kurtosis", 6.2)
    lateral_sigma = params.get("lateralSigma", 9.0)
    lateral_mu = params.get("lateralMu", 0.0)

    # A crystalline Si implant can show a channeling tail beyond the random
    # Pearson-IV peak. We represent that here with a second, normalized
    # exponential component. The defaults below target a practical
    # B / Si(100) / 10 keV / 7 degree tilt case, so the tail fraction is kept
    # deliberately small because the tilt suppresses channeling.
    tail_fraction = params.get("tailFraction", 0.03)
    tail_start_depth = params.get("tailStartDepth", 45.0)
    tail_decay_length = params.get("tailDecayLength", 55.0)
    tail_blend_width = params.get("tailBlendWidth", 20.0)

    structure = make_structure(
        x_extent,
        top_space,
        implant_depth,
        opening_width,
        mask_height,
        grid_delta,
    )
    cell_set = build_cell_set(structure, top_space)
    cell_set.writeVTU("initial.vtu")

    pearson = vcs.PearsonIVParameters()
    pearson.mu = projected_range
    pearson.sigma = depth_sigma
    pearson.gamma = skewness
    pearson.beta = kurtosis

    model = vcs.ImplantPearsonIV(pearson, lateral_mu, lateral_sigma)
    apply_masked_implant(
        cell_set,
        model,
        opening_width,
        mask_height,
        angle,
        beam_spacing,
        implant_depth,
        implant_dose_cm2,
        substrate_type,
        tail_fraction,
        tail_start_depth,
        tail_decay_length,
        tail_blend_width,
    )

    cell_set.writeVTU("final.vtu")
    print("Wrote initial.vtu and final.vtu")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
