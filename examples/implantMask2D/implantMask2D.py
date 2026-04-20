import sys

try:
    import viennacs as _vcs

    _vcs.setDimension(2)
    vcs = _vcs.d2
    import viennals as _vls

    _vls.setDimension(2)
    vls = _vls.d2
    MaterialMap = _vls.MaterialMap
    BoundaryConditionEnum = _vls.BoundaryConditionEnum
    BooleanOperationEnum = _vls.BooleanOperationEnum
except ImportError:
    import viennacs2d as vcs
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
    return cell_set


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
    species = str(params.get("species", "B"))
    material = str(params.get("material", "Si"))
    energy_kev = params.get("energyKeV", 10.0)
    rotation_deg = params.get("rotationDeg", 0.0)
    screen_thickness = params.get("screenThickness", 0.0)
    damage_level = params.get("damageLevel", 0.0)
    dose_control = str(params.get("doseControl", "WaferDose"))
    preferred_model = str(
        params.get(
            "preferredModel",
            "DualPearsonIV" if substrate_type == "crystalline" else "PearsonIV",
        )
    )

    if "_vcs" not in globals() or not hasattr(vcs, "RecipeDrivenImplantModel"):
        print("This example requires the unified viennacs package with RecipeDrivenImplantModel.")
        return 2

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
    recipe = _vcs.ImplantRecipe()
    recipe.species = species
    recipe.material = material
    recipe.substrateType = substrate_type
    recipe.preferredModel = preferred_model
    recipe.energyKeV = energy_kev
    recipe.tiltDeg = angle
    recipe.rotationDeg = rotation_deg
    recipe.screenThickness = screen_thickness
    recipe.damageLevel = damage_level
    recipe.tableFileName = str(
        params.get("tablePath", _vcs.getDefaultImplantTablePath())
    )
    recipe.useTableLookup = "projectedRange" not in params

    if not recipe.useTableLookup:
        entry = recipe.entry
        entry.modelType = preferred_model
        entry.headFraction = params.get("headFraction", 0.97)
        entry.screenDecayLength = params.get("screenDecayLength", 0.0)
        entry.damageDecay = params.get("damageDecay", 0.0)
        entry.tiltDecayDeg = params.get("tiltDecayDeg", 0.0)
        entry.reshapeStrength = params.get("reshapeStrength", 1.0)

        entry.headParams.mu = params.get("projectedRange", 31.0)
        entry.headParams.sigma = params.get("depthSigma", 11.0)
        entry.headParams.gamma = params.get("skewness", 1.15)
        entry.headParams.beta = params.get("kurtosis", 6.2)
        entry.headLateralMu = params.get("lateralMu", 0.0)
        entry.headLateralSigma = params.get("lateralSigma", 5.0)
        entry.headLateralModel = str(params.get("headLateralModel", params.get("lateralModel", "constant")))
        entry.headLateralScale = params.get("headLateralScale", 1.0)
        entry.headLateralLv = params.get("headLateralLv", 1.0)
        entry.headLateralDeltaSigma = params.get("headLateralDeltaSigma", 0.0)
        entry.headLateralP1 = params.get("headLateralP1", 0.0)
        entry.headLateralP2 = params.get("headLateralP2", 0.0)
        entry.headLateralP3 = params.get("headLateralP3", 0.0)
        entry.headLateralP4 = params.get("headLateralP4", 0.0)
        entry.headLateralP5 = params.get("headLateralP5", 0.0)

        entry.tailParams.mu = params.get("tailProjectedRange", 85.0)
        entry.tailParams.sigma = params.get("tailDepthSigma", 28.0)
        entry.tailParams.gamma = params.get("tailSkewness", 0.35)
        entry.tailParams.beta = params.get("tailKurtosis", 4.0)
        entry.tailLateralMu = params.get("tailLateralMu", entry.headLateralMu)
        entry.tailLateralSigma = params.get("tailLateralSigma", 6.0)
        entry.tailLateralModel = str(params.get("tailLateralModel", entry.headLateralModel))
        entry.tailLateralScale = params.get("tailLateralScale", entry.headLateralScale)
        entry.tailLateralLv = params.get("tailLateralLv", entry.headLateralLv)
        entry.tailLateralDeltaSigma = params.get(
            "tailLateralDeltaSigma", entry.headLateralDeltaSigma
        )
        entry.tailLateralP1 = params.get("tailLateralP1", entry.headLateralP1)
        entry.tailLateralP2 = params.get("tailLateralP2", entry.headLateralP2)
        entry.tailLateralP3 = params.get("tailLateralP3", entry.headLateralP3)
        entry.tailLateralP4 = params.get("tailLateralP4", entry.headLateralP4)
        entry.tailLateralP5 = params.get("tailLateralP5", entry.headLateralP5)

    model = vcs.RecipeDrivenImplantModel(recipe)
    if recipe.useTableLookup:
        print(f"Using implant table defaults from {recipe.tableFileName}")

    implant = vcs.Implant()
    implant.setCellSet(cell_set)
    implant.setImplantModel(model)
    implant.setImplantAngle(angle)
    if "_vcs" in globals() and hasattr(_vcs, "ImplantDoseControl"):
        dose_modes = {
            "off": _vcs.ImplantDoseControl.Off,
            "waferdose": _vcs.ImplantDoseControl.WaferDose,
            "beamdose": _vcs.ImplantDoseControl.BeamDose,
        }
        implant.setDoseControl(
            dose_modes.get(dose_control.strip().lower(), _vcs.ImplantDoseControl.WaferDose)
        )
    if hasattr(implant, "setDose"):
        implant.setDose(implant_dose_cm2)
    if hasattr(implant, "setLengthUnitInCm"):
        implant.setLengthUnitInCm(1.0e-7)
    if hasattr(implant, "enableBeamHits"):
        implant.enableBeamHits(True)
    if hasattr(implant, "setOutputConcentrationInCm3"):
        implant.setOutputConcentrationInCm3(True)
    implant.setMaskMaterials([2])
    implant.apply()

    if abs(beam_spacing - grid_delta) > 1e-9:
        print(
            f"Note: ViennaCS core Implant uses gridDelta={grid_delta} as the beam spacing. "
            f"The config value beamSpacing={beam_spacing} is currently not used by the core implant."
        )

    cell_set.writeVTU("final.vtu")
    print("Wrote initial.vtu and final.vtu")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
