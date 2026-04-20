import math
import sys

try:
    import viennacs as _vcs

    _vcs.setDimension(2)
    vcs = _vcs.d2
    PearsonIVParameters = _vcs.PearsonIVParameters
    import viennals as _vls

    _vls.setDimension(2)
    vls = _vls.d2
    MaterialMap = _vls.MaterialMap
    BoundaryConditionEnum = _vls.BoundaryConditionEnum
    BooleanOperationEnum = _vls.BooleanOperationEnum
except ImportError:
    import viennacs2d as vcs
    PearsonIVParameters = vcs.PearsonIVParameters
    import viennals2d as vls
    MaterialMap = vls.MaterialMap
    BoundaryConditionEnum = vls.BoundaryConditionEnum
    BooleanOperationEnum = vls.BooleanOperationEnum


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
        BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
        BoundaryConditionEnum.INFINITE_BOUNDARY,
    ]
    return vls.Domain(bounds, boundary_conditions, grid_delta)


def make_plane(domain, y_position: float):
    vls.MakeGeometry(domain, vls.Plane([0.0, y_position], [0.0, 1.0])).apply()


def make_box(domain, min_x: float, min_y: float, max_x: float, max_y: float):
    vls.MakeGeometry(domain, vls.Box([min_x, min_y], [max_x, max_y])).apply()


def add_level_set(level_sets, material_map, level_set, material_id: int, wrap_lower_level_set: bool = True):
    if level_sets and wrap_lower_level_set:
        vls.BooleanOperation(
            level_set, level_sets[-1], BooleanOperationEnum.UNION
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
    material_map = MaterialMap()
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
    vls.BooleanOperation(mask, opening, BooleanOperationEnum.RELATIVE_COMPLEMENT).apply()
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


def scale_concentration_to_cm3(cell_set, dose_cm2: float, beam_spacing: float):
    concentration = cell_set.getScalarData("concentration")
    dose_nm2 = dose_cm2 / 1.0e14
    dose_per_beam = dose_nm2 * beam_spacing
    scaled = [value * dose_per_beam * 1.0e21 for value in concentration]
    cell_set.setScalarData("concentration", scaled)


def compute_beam_hits(
    cell_set,
    opening_width: float,
    entry_height: float,
    angle_deg: float,
    beam_spacing: float,
):
    beam_hits = [0.0] * cell_set.getNumberOfCells()
    materials = cell_set.getScalarData("Material")
    angle = math.radians(angle_deg)

    for x_mask in frange(-0.5 * opening_width, 0.5 * opening_width, beam_spacing):
        for idx in range(cell_set.getNumberOfCells()):
            x_pos, y_pos, _ = cell_set.getCellCenter(idx)
            if materials[idx] != 1.0:
                continue

            depth_axis = -y_pos
            if depth_axis < 0.0:
                continue

            # Follow the centerline of the tilted beam from the top entry
            # plane down into the substrate so the beam-hit field matches the
            # apparent implant angle.
            x_beam = x_mask - math.tan(angle) * (entry_height + depth_axis)
            if abs(x_beam - x_pos) <= 0.5 * beam_spacing:
                beam_hits[idx] += 1.0

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

    # Crystalline Si is represented with a dual-Pearson model:
    # one Pearson-IV for the random component and one for the channeling tail.
    head_fraction = params.get("headFraction", 0.97)
    tail_projected_range = params.get("tailProjectedRange", 85.0)
    tail_depth_sigma = params.get("tailDepthSigma", 28.0)
    tail_skewness = params.get("tailSkewness", 0.35)
    tail_kurtosis = params.get("tailKurtosis", 3.2)

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

    pearson = PearsonIVParameters()
    pearson.mu = projected_range
    pearson.sigma = depth_sigma
    pearson.gamma = skewness
    pearson.beta = kurtosis

    tail_pearson = PearsonIVParameters()
    tail_pearson.mu = tail_projected_range
    tail_pearson.sigma = tail_depth_sigma
    tail_pearson.gamma = tail_skewness
    tail_pearson.beta = tail_kurtosis

    if substrate_type == "crystalline":
        model = vcs.ImplantDualPearsonIV(
            pearson,
            tail_pearson,
            head_fraction,
            lateral_mu,
            lateral_sigma,
        )
    else:
        model = vcs.ImplantPearsonIV(pearson, lateral_mu, lateral_sigma)

    implant = vcs.Implant()
    implant.setCellSet(cell_set)
    implant.setImplantModel(model)
    implant.setImplantAngle(angle)
    implant.setMaskMaterials([2])
    implant.apply()

    if abs(beam_spacing - grid_delta) > 1e-9:
        print(
            f"Note: ViennaCS core Implant uses gridDelta={grid_delta} as the beam spacing. "
            f"The config value beamSpacing={beam_spacing} is used only for the beamHits visualization."
        )

    scale_concentration_to_cm3(cell_set, implant_dose_cm2, grid_delta)
    compute_beam_hits(
        cell_set, opening_width, top_space + mask_height, angle, beam_spacing
    )

    cell_set.writeVTU("final.vtu")
    print("Wrote initial.vtu and final.vtu")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
