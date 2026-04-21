import sys
import math
import re
from pathlib import Path

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


def parse_float_list(value) -> list[float]:
    if value is None:
        return []
    if isinstance(value, (float, int)):
        return [float(value)]
    if isinstance(value, str):
        return [
            float(tok.strip())
            for tok in value.split(",")
            if tok.strip()
        ]
    return [float(v) for v in value]


def parse_bool(value, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value != 0
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off"}:
        return False
    return default


def _canonical_species_token(species: str) -> str:
    lut = {
        "b": "B",
        "boron": "B",
        "p": "P",
        "phosphorus": "P",
        "as": "As",
        "arsenic": "As",
        "sb": "Sb",
        "antimony": "Sb",
        "in": "In",
        "indium": "In",
    }
    return lut.get(species.strip().lower(), species.strip())


def _canonical_material_token(material: str) -> str:
    lut = {
        "si": "Si",
        "silicon": "Si",
        "ge": "Ge",
        "germanium": "Ge",
    }
    return lut.get(material.strip().lower(), material.strip())


def _eval_arrhenius_cm2_per_s(prefactor: float, activation_eV: float, temperatureK: float) -> float:
    kB_eV_per_K = 8.617333262145e-5
    t = max(float(temperatureK), 1.0)
    return float(prefactor) * math.exp(-float(activation_eV) / (kB_eV_per_K * t))


def load_advcal_anneal_defaults(
    file_name: str,
    material: str,
    species: str,
    temperatureK: float,
    length_unit_in_cm: float,
) -> dict:
    file_path = Path(file_name)
    if not file_path.exists():
        return {}

    text = file_path.read_text(encoding="utf-8", errors="ignore")
    mat = _canonical_material_token(material)
    sp = _canonical_species_token(species)

    result = {
        "tableFileName": str(file_path),
        "material": mat,
        "species": sp,
    }

    m_int = re.search(
        rf"pdbSetDouble\s+{re.escape(mat)}\s+Int\s+Di\s+\{{\s*\[Arr\s+([0-9eE+.\-]+)\s+([0-9eE+.\-]+)\]\s*\}}",
        text,
    )
    m_vac = re.search(
        rf"pdbSetDouble\s+{re.escape(mat)}\s+Vac\s+Dv\s+\{{\s*\[Arr\s+([0-9eE+.\-]+)\s+([0-9eE+.\-]+)\]\s*\}}",
        text,
    )

    # For dopant diffusivity, use the neutral BI pair Arrhenius entry:
    # pdbSetDoubleArray Si B Int D { 0 {[Arr D0 Ea]} ... }.
    d0_ea = re.search(
        rf"pdbSetDoubleArray\s+{re.escape(mat)}\s+{re.escape(sp)}\s+Int\s+D\s*\{{.*?\n\s*0\s+\{{\s*\[Arr\s+([0-9eE+.\-]+)\s+([0-9eE+.\-]+)\]\s*\}}",
        text,
        flags=re.DOTALL,
    )

    ikfi_factor = re.search(
        rf"pdbSet\s+{re.escape(mat)}\s+ICluster\s+Ikfi\s+\{{\s*\[expr\s+([0-9eE+.\-]+)\s*\*\s*\[pdbGet\s+{re.escape(mat)}\s+I\s+Di\]\]\s*\}}",
        text,
    )
    ikfc_arr = re.search(
        rf"pdbSet\s+{re.escape(mat)}\s+ICluster\s+Ikfc\s+\{{\s*\[expr\s+\[Arr\s+([0-9eE+.\-]+)\s+([0-9eE+.\-]+)\]\s*\*\s*\[pdbGet\s+{re.escape(mat)}\s+I\s+Di\]\]\s*\}}",
        text,
    )
    ikr_arr = re.search(
        rf"pdbSet\s+{re.escape(mat)}\s+ICluster\s+Ikr\s+\{{\s*\[Arr\s+([0-9eE+.\-]+)\s+([0-9eE+.\-]+)\]\s*\}}",
        text,
    )
    icluster_init_percent = re.search(
        rf"pdbSet\s+{re.escape(mat)}\s+ICluster\s+InitPercent\s+([0-9eE+.\-]+)",
        text,
    )

    unit_scale = 1.0 / max(float(length_unit_in_cm) ** 2, 1.0e-30)

    if d0_ea:
        d0_cm2_s = float(d0_ea.group(1))
        ea_eV = float(d0_ea.group(2))
        result["annealD0_cm2_per_s"] = d0_cm2_s
        result["annealEa_eV"] = ea_eV
        result["annealD0"] = d0_cm2_s * unit_scale
        result["annealEa"] = ea_eV

    if m_int:
        di_cm2_s = _eval_arrhenius_cm2_per_s(
            float(m_int.group(1)),
            float(m_int.group(2)),
            temperatureK,
        )
        result["annealInterstitialDiffusivity_cm2_per_s"] = di_cm2_s
        result["annealInterstitialDiffusivity"] = di_cm2_s * unit_scale

    if m_vac:
        dv_cm2_s = _eval_arrhenius_cm2_per_s(
            float(m_vac.group(1)),
            float(m_vac.group(2)),
            temperatureK,
        )
        result["annealVacancyDiffusivity_cm2_per_s"] = dv_cm2_s
        result["annealVacancyDiffusivity"] = dv_cm2_s * unit_scale

    # Approximate bulk I-V recombination coefficient from AdvCal expression:
    # k_bulk ~ 4*pi*Rcapture*(Di + Dv), with Rcapture = 5e-8 cm in AdvCal.
    if "annealInterstitialDiffusivity_cm2_per_s" in result and "annealVacancyDiffusivity_cm2_per_s" in result:
        r_capture_cm = 5.0e-8
        krec_cm3_s = 4.0 * math.pi * r_capture_cm * (
            result["annealInterstitialDiffusivity_cm2_per_s"]
            + result["annealVacancyDiffusivity_cm2_per_s"]
        )
        # Raw AdvCal value is too stiff for the current reduced I/V state model.
        # Use a bounded surrogate as default while keeping raw for traceability.
        result["annealDefectRecombinationRateAdvCalRaw"] = krec_cm3_s
        result["annealDefectRecombinationRate"] = min(krec_cm3_s, 1.0e-25)
        result["advcalBulkCaptureRadius_cm"] = r_capture_cm

    # One-moment interstitial cluster kinetics (surrogate mapping).
    if "annealInterstitialDiffusivity_cm2_per_s" in result:
        di_cm2_s = result["annealInterstitialDiffusivity_cm2_per_s"]
        if ikfi_factor:
            result["annealIClusterIkfi"] = float(ikfi_factor.group(1)) * di_cm2_s
        if ikfc_arr:
            ikfc_pref = _eval_arrhenius_cm2_per_s(
                float(ikfc_arr.group(1)),
                float(ikfc_arr.group(2)),
                temperatureK,
            )
            result["annealIClusterIkfc"] = ikfc_pref * di_cm2_s
        if ikr_arr:
            result["annealIClusterIkr"] = _eval_arrhenius_cm2_per_s(
                float(ikr_arr.group(1)),
                float(ikr_arr.group(2)),
                temperatureK,
            )
        if icluster_init_percent:
            result["annealIClusterInitFraction"] = max(
                0.0, min(1.0, float(icluster_init_percent.group(1)) / 100.0)
            )

    return result


def make_bounds(
    x_extent: float,
    top_space: float,
    substrate_depth: float,
    mask_height: float,
    oxide_thickness: float,
):
    return [
        -0.5 * x_extent,
        0.5 * x_extent,
        -substrate_depth,
        top_space + oxide_thickness + mask_height,
    ]


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
    substrate_depth: float,
    opening_width: float,
    mask_height: float,
    oxide_thickness: float,
    grid_delta: float,
):
    bounds = make_bounds(
        x_extent, top_space, substrate_depth, mask_height, oxide_thickness
    )
    material_map = MaterialMap()
    level_sets = []

    substrate_bottom = make_domain(bounds, grid_delta)
    make_plane(substrate_bottom, -substrate_depth)
    add_level_set(level_sets, material_map, substrate_bottom, 1)

    substrate_top = make_domain(bounds, grid_delta)
    make_plane(substrate_top, 0.0)
    add_level_set(level_sets, material_map, substrate_top, 1)

    if oxide_thickness > 0.0:
        oxide = make_domain(bounds, grid_delta)
        make_plane(oxide, oxide_thickness)
        add_level_set(level_sets, material_map, oxide, 3)

    mask = make_domain(bounds, grid_delta)
    make_plane(mask, oxide_thickness + mask_height)

    opening = make_domain(bounds, grid_delta)
    make_box(
        opening,
        -0.5 * opening_width,
        oxide_thickness - grid_delta,
        0.5 * opening_width,
        oxide_thickness + mask_height + grid_delta,
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
    config_dir = Path(sys.argv[1]).resolve().parent

    grid_delta = params.get("gridDelta", 1.0)
    x_extent = params.get("xExtent", 40.0)
    top_space = params.get("topSpace", 6.0)
    substrate_depth = params.get("substrateDepth", 12.0)
    opening_width = params.get("openingWidth", 8.0)
    mask_height = params.get("maskHeight", 4.0)
    oxide_thickness = params.get("oxideThickness", 0.0)
    angle = params.get("angle", 7.0)
    implant_dose_cm2 = params.get("doseCm2", 1.0e15)
    substrate_type = str(params.get("substrateType", "amorphous")).strip().lower()
    species = str(params.get("species", "B"))
    material = str(params.get("material", "Si"))
    energy_kev = params.get("energyKeV", 10.0)
    rotation_deg = params.get("rotationDeg", 0.0)
    screen_thickness = params.get("screenThickness", 0.0)
    damage_level = params.get("damageLevel", 0.0)
    length_unit_in_cm = params.get("lengthUnitInCm", 1.0e-7)
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
        substrate_depth,
        opening_width,
        mask_height,
        oxide_thickness,
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
    recipe.dosePerCm2 = implant_dose_cm2
    recipe.screenThickness = screen_thickness
    recipe.damageLevel = damage_level
    table_file = params.get("tablePath", _vcs.getDefaultImplantTablePath())
    table_file = Path(str(table_file))
    if not table_file.is_absolute():
        candidate = (config_dir / table_file).resolve()
        if candidate.exists():
            table_file = candidate
    recipe.tableFileName = str(table_file)
    recipe.useTableLookup = "projectedRange" not in params

    if not recipe.useTableLookup:
        entry = recipe.entry
        entry.modelType = preferred_model
        entry.headFraction = params.get("headFraction", 0.97)

        entry.headParams.mu = params.get("projectedRange", 31.0)
        entry.headParams.sigma = params.get("depthSigma", 11.0)
        entry.headParams.gamma = params.get("skewness", 1.15)
        entry.headParams.beta = params.get("kurtosis", 6.2)
        entry.headLateralMu = params.get("lateralMu", 0.0)
        entry.headLateralSigma = params.get("lateralSigma", 5.0)

        entry.tailParams.mu = params.get("tailProjectedRange", 85.0)
        entry.tailParams.sigma = params.get("tailDepthSigma", 28.0)
        entry.tailParams.gamma = params.get("tailSkewness", 0.35)
        entry.tailParams.beta = params.get("tailKurtosis", 4.0)
        entry.tailLateralMu = params.get("tailLateralMu", entry.headLateralMu)
        entry.tailLateralSigma = params.get("tailLateralSigma", 6.0)

    model = vcs.RecipeDrivenImplantModel(recipe)
    if recipe.useTableLookup:
        print(f"Using implant table defaults from {recipe.tableFileName}")

    damage_recipe = _vcs.DamageRecipe()
    damage_recipe.species = species
    damage_recipe.material = material
    damage_recipe.energyKeV = energy_kev
    damage_recipe.tiltDeg = angle
    damage_recipe.rotationDeg = rotation_deg
    damage_recipe.dosePerCm2 = implant_dose_cm2
    damage_recipe.screenThickness = screen_thickness

    default_damage_table = table_file
    if (
        default_damage_table.exists()
        and "_in_" in default_damage_table.name
        and "_damage_in_" not in default_damage_table.name
    ):
        candidate_name = default_damage_table.name.replace("_in_", "_damage_in_", 1)
        candidate = default_damage_table.with_name(candidate_name)
        if candidate.exists():
            default_damage_table = candidate

    damage_table_file = params.get("damageTablePath", str(default_damage_table))
    damage_table_file = Path(str(damage_table_file))
    if not damage_table_file.is_absolute():
        candidate = (config_dir / damage_table_file).resolve()
        if candidate.exists():
            damage_table_file = candidate
    damage_recipe.tableFileName = str(damage_table_file)
    damage_recipe.useTableLookup = "damageProjectedRange" not in params

    if not damage_recipe.useTableLookup:
        damage_entry = damage_recipe.entry
        damage_entry.projectedRange = params.get("damageProjectedRange", 30.0)
        damage_entry.verticalSigma = params.get("damageVerticalSigma", 10.0)
        setattr(damage_entry, "lambda", params.get("damageLambda", 20.0))
        damage_entry.defectsPerIon = params.get("damageDefectsPerIon", 100.0)
        damage_entry.lateralMu = params.get("damageLateralMu", 0.0)
        damage_entry.lateralSigma = params.get("damageLateralSigma", 10.0)
        damage_entry.lateralModel = str(params.get("damageLateralModel", "taurus"))
        damage_entry.lateralScale = params.get("damageLateralScale", 1.0)
        damage_entry.lateralLv = params.get("damageLateralLv", 1.0)
        damage_entry.lateralDeltaSigma = params.get("damageLateralDeltaSigma", 0.0)
        damage_entry.lateralP1 = params.get("damageLateralP1", 0.0)
        damage_entry.lateralP2 = params.get("damageLateralP2", 0.0)
        damage_entry.lateralP3 = params.get("damageLateralP3", 0.0)
        damage_entry.lateralP4 = params.get("damageLateralP4", 0.0)
        damage_entry.lateralP5 = params.get("damageLateralP5", 0.0)

    damage_model = vcs.RecipeDrivenDamageModel(damage_recipe)
    if damage_recipe.useTableLookup:
        print(f"Using damage table defaults from {damage_recipe.tableFileName}")

    implant = vcs.Implant()
    implant_species_label = str(params.get("implantSpeciesLabel", "concentration_annealed"))
    implant.setCellSet(cell_set)
    implant.setImplantModel(model)
    implant.setDamageModel(damage_model)
    if hasattr(implant, "setConcentrationLabel"):
        implant.setConcentrationLabel(implant_species_label)
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
        implant.setLengthUnitInCm(length_unit_in_cm)
    if hasattr(implant, "enableBeamHits"):
        implant.enableBeamHits(True)
    if hasattr(implant, "setOutputConcentrationInCm3"):
        implant.setOutputConcentrationInCm3(True)
    if hasattr(implant, "setDamageFactor"):
        implant.setDamageFactor(params.get("damageFactor", 1.0))
    implant.setMaskMaterials([2])
    if oxide_thickness > 0.0 and hasattr(implant, "setScreenMaterials"):
        implant.setScreenMaterials([3])
    implant.apply()

    # Preserve the as-implanted profile before any anneal step.
    as_implanted = list(cell_set.getScalarData(implant_species_label))
    cell_set.addScalarData("concentration_as_implanted", 0.0)
    cell_set.setScalarData("concentration_as_implanted", as_implanted)

    anneal_step_durations = parse_float_list(params.get("annealStepDurations"))
    anneal_duration = params.get(
        "annealDuration",
        sum(anneal_step_durations) if anneal_step_durations else 0.0,
    )
    if (anneal_duration > 0.0 or anneal_step_durations) and hasattr(vcs, "Anneal"):
        anneal_reference_temperature = params.get("annealTemperature", 1273.15)
        parsed_schedule_temperatures = parse_float_list(params.get("annealTemperatures"))
        if not parsed_schedule_temperatures:
            parsed_schedule_temperatures = parse_float_list(params.get("annealRampTemperatures"))
        if parsed_schedule_temperatures:
            anneal_reference_temperature = parsed_schedule_temperatures[0]

        anneal_table_defaults = {}
        if parse_bool(params.get("annealUseTableLookup", 1), True):
            anneal_table_file = params.get(
                "annealTablePath",
                _vcs.getDefaultAnnealTablePath()
                if hasattr(_vcs, "getDefaultAnnealTablePath")
                else "data/AnnealData/AdvCal_2023.12.fps",
            )
            anneal_table_file = Path(str(anneal_table_file))
            if not anneal_table_file.is_absolute():
                candidate = (config_dir / anneal_table_file).resolve()
                if candidate.exists():
                    anneal_table_file = candidate
                else:
                    repo_candidate = (Path(__file__).resolve().parents[2] / anneal_table_file).resolve()
                    if repo_candidate.exists():
                        anneal_table_file = repo_candidate
            anneal_table_defaults = load_advcal_anneal_defaults(
                str(anneal_table_file),
                material=material,
                species=species,
                temperatureK=anneal_reference_temperature,
                length_unit_in_cm=length_unit_in_cm,
            )
            if anneal_table_defaults:
                print(
                    f"Using anneal defaults from {anneal_table_defaults.get('tableFileName')}"
                )
                raw_krec = anneal_table_defaults.get("annealDefectRecombinationRateAdvCalRaw")
                used_krec = anneal_table_defaults.get("annealDefectRecombinationRate")
                if raw_krec is not None and used_krec is not None and raw_krec > used_krec:
                    print(
                        f"Using bounded defect recombination surrogate {used_krec:.3e} "
                        f"(raw AdvCal value {raw_krec:.3e})."
                    )

        anneal = vcs.Anneal()
        anneal.setCellSet(cell_set)
        anneal.setSpeciesLabel(str(params.get("annealSpeciesLabel", implant_species_label)))
        anneal.setDuration(anneal_duration)
        anneal_mode = str(params.get("annealMode", "explicit")).strip().lower()
        if hasattr(_vcs, "AnnealMode"):
            if anneal_mode == "implicit":
                anneal.setMode(_vcs.AnnealMode.Implicit)
                if "annealImplicitMaxIterations" in params or "annealImplicitTolerance" in params:
                    anneal.setImplicitSolverOptions(
                        int(params.get("annealImplicitMaxIterations", 400)),
                        params.get("annealImplicitTolerance", 1.0e-6),
                    )
            else:
                anneal.setMode(_vcs.AnnealMode.Explicit)
        if "annealTimeStep" in params:
            anneal.setTimeStep(params.get("annealTimeStep"))
        if "annealStabilityFactor" in params:
            anneal.setStabilityFactor(params.get("annealStabilityFactor"))

        if "annealDiffusionCoefficient" in params:
            anneal.setDiffusionCoefficient(params.get("annealDiffusionCoefficient"))
        elif "annealD0" in params and "annealEa" in params:
            anneal.setArrheniusParameters(params.get("annealD0"), params.get("annealEa"))
            anneal.setTemperature(params.get("annealTemperature", 1273.15))
        elif "annealD0" in anneal_table_defaults and "annealEa" in anneal_table_defaults:
            anneal.setArrheniusParameters(
                anneal_table_defaults.get("annealD0"),
                anneal_table_defaults.get("annealEa"),
            )
            anneal.setTemperature(params.get("annealTemperature", anneal_reference_temperature))
        else:
            print("Anneal requested but no diffusion model parameters found; skipping anneal.")
            anneal_duration = 0.0

        if anneal_duration > 0.0 or anneal_step_durations:
            step_durations = anneal_step_durations
            step_temperatures = parsed_schedule_temperatures

            if step_durations and hasattr(anneal, "clearTemperatureSchedule"):
                anneal.clearTemperatureSchedule()
                if step_temperatures and len(step_temperatures) == len(step_durations) + 1:
                    for i, step_duration in enumerate(step_durations):
                        anneal.addRampStep(
                            step_duration,
                            step_temperatures[i],
                            step_temperatures[i + 1],
                        )
                    print(
                        f"Using {len(step_durations)}-step temperature ramp anneal schedule."
                    )
                elif step_temperatures and len(step_temperatures) == len(step_durations):
                    for step_duration, temp in zip(step_durations, step_temperatures):
                        anneal.addIsothermalStep(step_duration, temp)
                    print(
                        f"Using {len(step_durations)}-step isothermal anneal schedule."
                    )
                else:
                    fallback_temp = params.get("annealTemperature", 1273.15)
                    for step_duration in step_durations:
                        anneal.addIsothermalStep(step_duration, fallback_temp)
                    print(
                        "Using multi-step anneal schedule with shared annealTemperature."
                    )

            anneal.setDiffusionMaterials([1])
            anneal.setBlockingMaterials([2])
            if parse_bool(params.get("annealDefectCoupling", 0.0), False):
                anneal.enableDefectCoupling(True)
                anneal.setDamageLabels(
                    str(params.get("annealDamageLabel", "Damage")),
                    str(params.get("annealDamageLastImpLabel", "Damage_LastImp")),
                )
                anneal.setDefectLabels(
                    str(params.get("annealInterstitialLabel", "Interstitial")),
                    str(params.get("annealVacancyLabel", "Vacancy")),
                )
                anneal.setDefectSourceWeights(
                    params.get("annealDefectFromDamageHistoryWeight", 0.0),
                    params.get("annealDefectFromDamageLastImpWeight", 1.0),
                )
                anneal.setDefectPartition(
                    params.get("annealDefectToInterstitial", 0.5),
                    params.get("annealDefectToVacancy", 0.5),
                )
                anneal.setDefectDiffusivities(
                    params.get(
                        "annealInterstitialDiffusivity",
                        anneal_table_defaults.get("annealInterstitialDiffusivity", 0.0),
                    ),
                    params.get(
                        "annealVacancyDiffusivity",
                        anneal_table_defaults.get("annealVacancyDiffusivity", 0.0),
                    ),
                )
                anneal.setDefectReactionRates(
                    params.get(
                        "annealDefectRecombinationRate",
                        anneal_table_defaults.get("annealDefectRecombinationRate", 0.0),
                    ),
                    params.get("annealInterstitialSinkRate", 0.0),
                    params.get("annealVacancySinkRate", 0.0),
                )
                anneal.setDefectEnhancedDiffusion(
                    params.get("annealTEDCoefficient", 0.0),
                    params.get("annealTEDNormalization", 1.0e20),
                )
                if parse_bool(params.get("annealDefectClustering", 0), False):
                    ikfi = params.get(
                        "annealIClusterIkfi",
                        anneal_table_defaults.get("annealIClusterIkfi", 0.0),
                    )
                    ikfc = params.get(
                        "annealIClusterIkfc",
                        anneal_table_defaults.get("annealIClusterIkfc", 0.0),
                    )
                    ikr = params.get(
                        "annealIClusterIkr",
                        anneal_table_defaults.get("annealIClusterIkr", 0.0),
                    )
                    if any(v > 0.0 for v in (ikfi, ikfc, ikr)):
                        anneal.enableDefectClustering(True)
                        anneal.setDefectClusterLabel(
                            str(params.get("annealIClusterLabel", "ICluster"))
                        )
                        anneal.setDefectClusterKinetics(ikfi, ikfc, ikr)
                        anneal.setDefectClusterInitFraction(
                            params.get(
                                "annealIClusterInitFraction",
                                anneal_table_defaults.get(
                                    "annealIClusterInitFraction", 0.0
                                ),
                            )
                        )
            anneal.apply()
            print("Anneal step completed.")

    # Store annealed profile explicitly as a separate scalar.
    annealed = list(cell_set.getScalarData(implant_species_label))
    cell_set.setScalarData("concentration_annealed", annealed)

    cell_set.writeVTU("final.vtu")
    print("Wrote initial.vtu and final.vtu")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
