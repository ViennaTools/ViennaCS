import sys
import math
import re
import csv
import json
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


def _load_json_file(file_name: str) -> dict | None:
    file_path = Path(file_name)
    if not file_path.exists():
        return None
    try:
        return json.loads(file_path.read_text(encoding="utf-8"))
    except Exception:
        return None


def resolve_existing_path(path_like, config_dir: Path) -> Path:
    path = Path(str(path_like))
    if path.is_absolute():
        return path

    script_dir = Path(__file__).resolve().parent
    repo_root = Path(__file__).resolve().parents[2]
    candidates = [
        (config_dir / path).resolve(),
        (script_dir / path).resolve(),
        (repo_root / path).resolve(),
    ]
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]


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


def _score_material_dir(material: str) -> str:
    mat = _canonical_material_token(material)
    lut = {
        "Si": "Silicon",
        "Ge": "Germanium",
    }
    return lut.get(mat, mat)


def _score_species_dir(species: str) -> str:
    sp = _canonical_species_token(species)
    lut = {
        "B": "Boron",
        "P": "Phosphorus",
        "As": "Arsenic",
        "Sb": "Antimony",
        "In": "Indium",
    }
    return lut.get(sp, sp)


def _parse_score_info_scalar(info_text: str, key: str) -> float | None:
    m = re.search(
        rf"array set \$Base \{{{re.escape(key)}\s+\{{Double\s+\{{([0-9eE+.\-]+)\}}\}}\}}",
        info_text,
    )
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def load_score_implant_factors(
    score_params_root: str,
    material: str,
    species: str,
) -> dict:
    mat_dir = _score_material_dir(material)
    sp_dir = _score_species_dir(species)
    info_file = Path(score_params_root) / mat_dir / sp_dir / "Info"
    if not info_file.exists():
        return {}

    text = info_file.read_text(encoding="utf-8", errors="ignore")
    out = {
        "scoreInfoFile": str(info_file),
    }
    for key in ("IFactor", "VFactor", "DFactor", "MCIFactor", "MCVFactor"):
        value = _parse_score_info_scalar(text, key)
        if value is not None:
            out[key] = value
    return out


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
    m_int_cstar = re.search(
        rf"pdbSet\s+{re.escape(mat)}\s+Int\s+Cstar\s+\{{\s*\[Arr\s+([0-9eE+.\-]+)\s+([0-9eE+.\-]+)\]\s*\}}",
        text,
    )
    m_vac_cstar = re.search(
        rf"pdbSet\s+{re.escape(mat)}\s+Vac\s+Cstar\s+\{{\s*\[Arr\s+([0-9eE+.\-]+)\s+([0-9eE+.\-]+)\]\s*\}}",
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
        di_d0 = float(m_int.group(1))
        di_ea = float(m_int.group(2))
        di_cm2_s = _eval_arrhenius_cm2_per_s(
            di_d0,
            di_ea,
            temperatureK,
        )
        result["annealInterstitialDiffusivityD0_cm2_per_s"] = di_d0
        result["annealInterstitialDiffusivityEa_eV"] = di_ea
        result["annealInterstitialDiffusivity_cm2_per_s"] = di_cm2_s
        result["annealInterstitialDiffusivity"] = di_cm2_s * unit_scale

    if m_vac:
        dv_d0 = float(m_vac.group(1))
        dv_ea = float(m_vac.group(2))
        dv_cm2_s = _eval_arrhenius_cm2_per_s(
            dv_d0,
            dv_ea,
            temperatureK,
        )
        result["annealVacancyDiffusivityD0_cm2_per_s"] = dv_d0
        result["annealVacancyDiffusivityEa_eV"] = dv_ea
        result["annealVacancyDiffusivity_cm2_per_s"] = dv_cm2_s
        result["annealVacancyDiffusivity"] = dv_cm2_s * unit_scale

    if m_int_cstar:
        ci_c0 = float(m_int_cstar.group(1))
        ci_ea = float(m_int_cstar.group(2))
        cstar_i_cm3 = _eval_arrhenius_cm2_per_s(
            ci_c0,
            ci_ea,
            temperatureK,
        )
        result["annealInterstitialEquilibriumC0_cm3"] = ci_c0
        result["annealInterstitialEquilibriumEa_eV"] = ci_ea
        result["annealInterstitialEquilibrium_cm3"] = cstar_i_cm3
        result["annealInterstitialEquilibrium"] = cstar_i_cm3

    if m_vac_cstar:
        cv_c0 = float(m_vac_cstar.group(1))
        cv_ea = float(m_vac_cstar.group(2))
        cstar_v_cm3 = _eval_arrhenius_cm2_per_s(
            cv_c0,
            cv_ea,
            temperatureK,
        )
        result["annealVacancyEquilibriumC0_cm3"] = cv_c0
        result["annealVacancyEquilibriumEa_eV"] = cv_ea
        result["annealVacancyEquilibrium_cm3"] = cstar_v_cm3
        result["annealVacancyEquilibrium"] = cstar_v_cm3

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

    # Approximate interface sink rates from AdvCal Ksurf = D / Lrec.
    # For reduced 0D sink terms: k_sink ~ D / Lrec^2.
    m_ksurf_i = re.search(
        rf"pdbSet\s+Ox_{re.escape(mat)}\s+I\s+Ksurf\s+\{{\s*\[expr\s+\[pdbGet\s+{re.escape(mat)}\s+I\s+Di\]\*([0-9eE+.\-]+)\]\s*\}}",
        text,
    )
    m_ksurf_v = re.search(
        rf"pdbSet\s+Ox_{re.escape(mat)}\s+V\s+Ksurf\s+\{{\s*\[expr\s+\[pdbGet\s+{re.escape(mat)}\s+V\s+Dv\]\*([0-9eE+.\-]+)\]\s*\}}",
        text,
    )
    if (
        m_ksurf_i
        and m_ksurf_v
        and "annealInterstitialDiffusivity_cm2_per_s" in result
        and "annealVacancyDiffusivity_cm2_per_s" in result
    ):
        mi = float(m_ksurf_i.group(1))
        mv = float(m_ksurf_v.group(1))
        li_cm = 1.0 / max(mi, 1.0e-30)
        lv_cm = 1.0 / max(mv, 1.0e-30)
        ki_raw = result["annealInterstitialDiffusivity_cm2_per_s"] / max(li_cm * li_cm, 1.0e-30)
        kv_raw = result["annealVacancyDiffusivity_cm2_per_s"] / max(lv_cm * lv_cm, 1.0e-30)
        cap = 1.0e-4
        result["annealInterstitialSinkRateFromKsurfRaw"] = ki_raw
        result["annealVacancySinkRateFromKsurfRaw"] = kv_raw
        result["annealInterstitialSinkRateFromKsurf"] = min(ki_raw, cap)
        result["annealVacancySinkRateFromKsurf"] = min(kv_raw, cap)
        result["advcalKsurfRecombinationLengthI_cm"] = li_cm
        result["advcalKsurfRecombinationLengthV_cm"] = lv_cm

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


def load_internal_implant_defaults(
    file_name: str,
    species: str,
    material: str,
    substrate_type: str,
) -> dict:
    db = _load_json_file(file_name)
    if not db:
        return {}

    sp = _canonical_species_token(species)
    mat = _canonical_material_token(material)
    sub = str(substrate_type).strip().lower()
    entries = db.get("entries", [])
    if not isinstance(entries, list):
        return {}

    for entry in entries:
        if not isinstance(entry, dict):
            continue
        if _canonical_species_token(str(entry.get("species", ""))) != sp:
            continue
        if _canonical_material_token(str(entry.get("material", ""))) != mat:
            continue
        allowed_sub = [
            str(v).strip().lower() for v in entry.get("substrateTypes", [])
        ]
        if allowed_sub and sub not in allowed_sub:
            continue
        return entry
    return {}


def load_internal_anneal_defaults(
    file_name: str,
    species: str,
    material: str,
    temperatureK: float,
    length_unit_in_cm: float,
) -> dict:
    db = _load_json_file(file_name)
    if not db:
        return {}

    sp = _canonical_species_token(species)
    mat = _canonical_material_token(material)
    entries = db.get("entries", [])
    if not isinstance(entries, list):
        return {}

    selected = None
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        if _canonical_material_token(str(entry.get("material", ""))) != mat:
            continue
        species_list = [
            _canonical_species_token(str(s)) for s in entry.get("species", [])
        ]
        if species_list and sp not in species_list:
            continue
        selected = entry
        break

    if selected is None:
        return {}

    t = max(float(temperatureK), 1.0)
    unit_scale_area = 1.0 / max(float(length_unit_in_cm) ** 2, 1.0e-30)
    out = {
        "tableFileName": str(file_name),
        "material": mat,
        "species": sp,
        "entryId": str(selected.get("id", "")),
    }

    dopant = selected.get("dopantDiffusion", {})
    if isinstance(dopant, dict):
        d0 = float(dopant.get("D0_cm2_per_s", 0.0))
        ea = float(dopant.get("Ea_eV", 0.0))
        if d0 > 0.0:
            out["annealD0_cm2_per_s"] = d0
            out["annealD0"] = d0 * unit_scale_area
            out["annealEa_eV"] = ea
            out["annealEa"] = ea

    defects = selected.get("defects", {})
    if isinstance(defects, dict):
        interstitial = defects.get("interstitial", {})
        vacancy = defects.get("vacancy", {})

        if isinstance(interstitial, dict):
            di = interstitial.get("Di", {})
            if isinstance(di, dict):
                di_d0 = float(di.get("D0_cm2_per_s", 0.0))
                di_ea = float(di.get("Ea_eV", 0.0))
                di_cm2_s = _eval_arrhenius_cm2_per_s(
                    di_d0,
                    di_ea,
                    t,
                )
                out["annealInterstitialDiffusivityD0_cm2_per_s"] = di_d0
                out["annealInterstitialDiffusivityEa_eV"] = di_ea
                out["annealInterstitialDiffusivity_cm2_per_s"] = di_cm2_s
                out["annealInterstitialDiffusivity"] = di_cm2_s * unit_scale_area
            cstar = interstitial.get("Cstar", {})
            if isinstance(cstar, dict):
                ci_c0 = float(cstar.get("C0_cm3", 0.0))
                ci_ea = float(cstar.get("Ea_eV", 0.0))
                ci_cm3 = _eval_arrhenius_cm2_per_s(
                    ci_c0,
                    ci_ea,
                    t,
                )
                out["annealInterstitialEquilibriumC0_cm3"] = ci_c0
                out["annealInterstitialEquilibriumEa_eV"] = ci_ea
                out["annealInterstitialEquilibrium_cm3"] = ci_cm3
                out["annealInterstitialEquilibrium"] = ci_cm3

        if isinstance(vacancy, dict):
            dv = vacancy.get("Dv", {})
            if isinstance(dv, dict):
                dv_d0 = float(dv.get("D0_cm2_per_s", 0.0))
                dv_ea = float(dv.get("Ea_eV", 0.0))
                dv_cm2_s = _eval_arrhenius_cm2_per_s(
                    dv_d0,
                    dv_ea,
                    t,
                )
                out["annealVacancyDiffusivityD0_cm2_per_s"] = dv_d0
                out["annealVacancyDiffusivityEa_eV"] = dv_ea
                out["annealVacancyDiffusivity_cm2_per_s"] = dv_cm2_s
                out["annealVacancyDiffusivity"] = dv_cm2_s * unit_scale_area
            cstar = vacancy.get("Cstar", {})
            if isinstance(cstar, dict):
                cv_c0 = float(cstar.get("C0_cm3", 0.0))
                cv_ea = float(cstar.get("Ea_eV", 0.0))
                cv_cm3 = _eval_arrhenius_cm2_per_s(
                    cv_c0,
                    cv_ea,
                    t,
                )
                out["annealVacancyEquilibriumC0_cm3"] = cv_c0
                out["annealVacancyEquilibriumEa_eV"] = cv_ea
                out["annealVacancyEquilibrium_cm3"] = cv_cm3
                out["annealVacancyEquilibrium"] = cv_cm3

        rec = defects.get("recombination", {})
        if (
            isinstance(rec, dict)
            and rec.get("mode", "").strip().lower() == "capture_radius"
            and "annealInterstitialDiffusivity_cm2_per_s" in out
            and "annealVacancyDiffusivity_cm2_per_s" in out
        ):
            r_capture_cm = float(rec.get("radius_cm", 5.0e-8))
            krec_cm3_s = 4.0 * math.pi * r_capture_cm * (
                out["annealInterstitialDiffusivity_cm2_per_s"]
                + out["annealVacancyDiffusivity_cm2_per_s"]
            )
            cap = float(rec.get("capForReducedModel_cm3_per_s", 1.0e-25))
            out["annealDefectRecombinationRateAdvCalRaw"] = krec_cm3_s
            out["annealDefectRecombinationRate"] = min(krec_cm3_s, cap)

    clusters = selected.get("interstitialCluster", {})
    if isinstance(clusters, dict):
        ikfi = clusters.get("Ikfi", {})
        ikfc = clusters.get("Ikfc", {})
        ikr = clusters.get("Ikr", {})
        if isinstance(ikfi, dict):
            out["annealIClusterIkfi"] = _eval_arrhenius_cm2_per_s(
                float(ikfi.get("D0", 0.0)),
                float(ikfi.get("Ea_eV", 0.0)),
                t,
            )
        if isinstance(ikfc, dict):
            out["annealIClusterIkfc"] = _eval_arrhenius_cm2_per_s(
                float(ikfc.get("D0", 0.0)),
                float(ikfc.get("Ea_eV", 0.0)),
                t,
            )
        if isinstance(ikr, dict):
            out["annealIClusterIkr"] = _eval_arrhenius_cm2_per_s(
                float(ikr.get("D0", 0.0)),
                float(ikr.get("Ea_eV", 0.0)),
                t,
            )
        if "InitPercent" in clusters:
            out["annealIClusterInitFraction"] = max(
                0.0, min(1.0, float(clusters.get("InitPercent", 0.0)) / 100.0)
            )

    implant_factors = selected.get("implantDamageFactors", {})
    if isinstance(implant_factors, dict):
        for key in ("IFactor", "VFactor", "DFactor", "MCIFactor", "MCVFactor"):
            if key in implant_factors:
                out[f"score{key}"] = float(implant_factors.get(key, 0.0))
        if "source" in implant_factors:
            out["scoreFactorsSource"] = str(implant_factors.get("source", ""))

    return out


def collect_defect_diagnostics(
    cell_set,
    material_label: str,
    substrate_material_id: int,
    interstitial_label: str,
    vacancy_label: str,
    ieq: float,
    veq: float,
    time_s: float,
    step_index: int,
    temperatureK: float,
) -> dict | None:
    interstitial = cell_set.getScalarData(interstitial_label)
    vacancy = cell_set.getScalarData(vacancy_label)
    materials = cell_set.getScalarData(material_label)
    if interstitial is None or vacancy is None:
        return None

    idxs = range(len(interstitial))
    if materials is not None:
        idxs = [
            i
            for i, m in enumerate(materials)
            if int(round(m)) == substrate_material_id
        ]
    if not idxs:
        return None

    sum_i = 0.0
    sum_v = 0.0
    sum_i_over_ieq = 0.0
    sum_v_over_veq = 0.0
    sum_iv_over_eq = 0.0
    min_i = float("inf")
    max_i = 0.0
    min_v = float("inf")
    max_v = 0.0

    for i in idxs:
        i_val = max(float(interstitial[i]), 0.0)
        v_val = max(float(vacancy[i]), 0.0)
        min_i = min(min_i, i_val)
        max_i = max(max_i, i_val)
        min_v = min(min_v, v_val)
        max_v = max(max_v, v_val)
        sum_i += i_val
        sum_v += v_val
        sum_i_over_ieq += (i_val / ieq) if ieq > 0.0 else float("nan")
        sum_v_over_veq += (v_val / veq) if veq > 0.0 else float("nan")
        if ieq > 0.0 and veq > 0.0:
            sum_iv_over_eq += (i_val * v_val) / (ieq * veq)
        else:
            sum_iv_over_eq += float("nan")

    count = float(len(idxs))
    row = {
        "step": step_index,
        "time_s": time_s,
        "temperature_K": temperatureK,
        "I_mean": sum_i / count,
        "V_mean": sum_v / count,
        "I_min": min_i,
        "I_max": max_i,
        "V_min": min_v,
        "V_max": max_v,
        "I_over_Ieq_mean": (sum_i_over_ieq / count) if ieq > 0.0 else float("nan"),
        "V_over_Veq_mean": (sum_v_over_veq / count) if veq > 0.0 else float("nan"),
        "IV_over_IeqVeq_mean": (sum_iv_over_eq / count)
        if ieq > 0.0 and veq > 0.0
        else float("nan"),
        "Ieq": ieq,
        "Veq": veq,
    }
    return row


def resolve_defect_coupling_defaults(
    params: dict,
    anneal_table_defaults: dict,
    score_implant_factors: dict,
) -> dict:
    resolved = {}

    # Damage-to-defect source weights
    resolved["annealDefectFromDamageHistoryWeight"] = params.get(
        "annealDefectFromDamageHistoryWeight", 0.0
    )
    resolved["annealDefectFromDamageLastImpWeight"] = params.get(
        "annealDefectFromDamageLastImpWeight", 1.0
    )

    # I/V partition from SCORE implant generation factors.
    i_fac = score_implant_factors.get("MCIFactor", score_implant_factors.get("IFactor", 1.0))
    v_fac = score_implant_factors.get("MCVFactor", score_implant_factors.get("VFactor", 1.0))
    i_fac = max(float(i_fac), 0.0)
    v_fac = max(float(v_fac), 0.0)
    if i_fac + v_fac <= 0.0:
        i_part = 0.5
        v_part = 0.5
    else:
        i_part = i_fac / (i_fac + v_fac)
        v_part = v_fac / (i_fac + v_fac)
    resolved["annealDefectToInterstitial"] = params.get("annealDefectToInterstitial", i_part)
    resolved["annealDefectToVacancy"] = params.get("annealDefectToVacancy", v_part)

    # Defect transport
    resolved["annealInterstitialDiffusivity"] = params.get(
        "annealInterstitialDiffusivity",
        anneal_table_defaults.get("annealInterstitialDiffusivity", 0.0),
    )
    resolved["annealVacancyDiffusivity"] = params.get(
        "annealVacancyDiffusivity",
        anneal_table_defaults.get("annealVacancyDiffusivity", 0.0),
    )
    resolved["annealDefectRecombinationRate"] = params.get(
        "annealDefectRecombinationRate",
        anneal_table_defaults.get("annealDefectRecombinationRate", 0.0),
    )

    # Sink rates from Ksurf-derived defaults if available; bounded for reduced model stability.
    sink_cap = float(params.get("annealSinkRateCap", 1.0e-4))
    if (
        "annealInterstitialSinkRateFromKsurfRaw" not in anneal_table_defaults
        and "annealVacancySinkRateFromKsurfRaw" not in anneal_table_defaults
    ):
        # If no explicit Ksurf metadata exists, use AdvCal-like 1 nm recombination
        # length as a reasonable default for Ox/Si interfaces.
        lrec_cm = float(params.get("annealKsurfRecombinationLengthCm", 1.0e-7))
        lrec_cm = max(lrec_cm, 1.0e-12)
        di = float(anneal_table_defaults.get("annealInterstitialDiffusivity_cm2_per_s", 0.0))
        dv = float(anneal_table_defaults.get("annealVacancyDiffusivity_cm2_per_s", 0.0))
        anneal_table_defaults["annealInterstitialSinkRateFromKsurfRaw"] = di / (lrec_cm * lrec_cm)
        anneal_table_defaults["annealVacancySinkRateFromKsurfRaw"] = dv / (lrec_cm * lrec_cm)
        anneal_table_defaults["advcalKsurfRecombinationLengthI_cm"] = lrec_cm
        anneal_table_defaults["advcalKsurfRecombinationLengthV_cm"] = lrec_cm
    default_ki = anneal_table_defaults.get(
        "annealInterstitialSinkRateFromKsurf",
        min(float(anneal_table_defaults.get("annealInterstitialSinkRateFromKsurfRaw", 0.0)), sink_cap),
    )
    default_kv = anneal_table_defaults.get(
        "annealVacancySinkRateFromKsurf",
        min(float(anneal_table_defaults.get("annealVacancySinkRateFromKsurfRaw", 0.0)), sink_cap),
    )
    resolved["annealInterstitialSinkRate"] = params.get("annealInterstitialSinkRate", default_ki)
    resolved["annealVacancySinkRate"] = params.get("annealVacancySinkRate", default_kv)

    # TED defaults from implant damage generation factor.
    d_factor = float(score_implant_factors.get("DFactor", 1.0))
    ted_coeff_default = max(0.0, 0.5 * d_factor)
    resolved["annealTEDCoefficient"] = params.get("annealTEDCoefficient", ted_coeff_default)
    resolved["annealTEDNormalization"] = params.get("annealTEDNormalization", 1.0e20)

    return resolved


def resolve_equilibrium_at_temperature(
    params: dict,
    anneal_table_defaults: dict,
    temperatureK: float,
    fallback_ieq: float,
    fallback_veq: float,
) -> tuple[float, float]:
    if "annealInterstitialEquilibrium" in params:
        ieq = float(params.get("annealInterstitialEquilibrium"))
    elif (
        "annealInterstitialEquilibriumC0_cm3" in anneal_table_defaults
        and "annealInterstitialEquilibriumEa_eV" in anneal_table_defaults
    ):
        ieq = _eval_arrhenius_cm2_per_s(
            float(anneal_table_defaults.get("annealInterstitialEquilibriumC0_cm3", 0.0)),
            float(anneal_table_defaults.get("annealInterstitialEquilibriumEa_eV", 0.0)),
            temperatureK,
        )
    else:
        ieq = float(fallback_ieq)

    if "annealVacancyEquilibrium" in params:
        veq = float(params.get("annealVacancyEquilibrium"))
    elif (
        "annealVacancyEquilibriumC0_cm3" in anneal_table_defaults
        and "annealVacancyEquilibriumEa_eV" in anneal_table_defaults
    ):
        veq = _eval_arrhenius_cm2_per_s(
            float(anneal_table_defaults.get("annealVacancyEquilibriumC0_cm3", 0.0)),
            float(anneal_table_defaults.get("annealVacancyEquilibriumEa_eV", 0.0)),
            temperatureK,
        )
    else:
        veq = float(fallback_veq)

    return ieq, veq


def write_diagnostics_csv(file_name: Path, rows: list[dict]):
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with file_name.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_diagnostics_plot(file_name: Path, rows: list[dict]):
    if not rows:
        return
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    t = [float(r["time_s"]) for r in rows]
    iv = [float(r["IV_over_IeqVeq_mean"]) for r in rows]
    ii = [float(r["I_over_Ieq_mean"]) for r in rows]
    vv = [float(r["V_over_Veq_mean"]) for r in rows]

    plt.figure(figsize=(8, 4.5))
    plt.plot(t, ii, label="I/Ieq")
    plt.plot(t, vv, label="V/Veq")
    plt.plot(t, iv, label="I*V/(Ieq*Veq)")
    plt.xlabel("Anneal Time [s]")
    plt.ylabel("Dimensionless ratio")
    plt.yscale("log")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(file_name, dpi=150)
    plt.close()


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
    score_implant_factors = {}

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
    implant_internal_defaults = {}
    if parse_bool(params.get("implantUseInternalDb", 1), True):
        implant_db_file = params.get(
            "implantDbFile",
            "vsclib/implant/si_taurus_defaults.json",
        )
        implant_db_file = resolve_existing_path(implant_db_file, config_dir)
        implant_internal_defaults = load_internal_implant_defaults(
            str(implant_db_file),
            species=species,
            material=material,
            substrate_type=substrate_type,
        )
        if implant_internal_defaults:
            print(
                f"Using internal implant defaults from {implant_db_file} "
                f"({implant_internal_defaults.get('id', 'entry')})"
            )

    if "preferredModel" not in params and implant_internal_defaults.get("preferredModel"):
        recipe.preferredModel = str(implant_internal_defaults.get("preferredModel"))

    table_path_default = (
        implant_internal_defaults.get("tableFile")
        if implant_internal_defaults.get("tableFile")
        else _vcs.getDefaultImplantTablePath()
    )
    table_file = params.get("tablePath", table_path_default)
    table_file = resolve_existing_path(table_file, config_dir)
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

    damage_table_file = params.get(
        "damageTablePath",
        implant_internal_defaults.get("damageTableFile", str(default_damage_table)),
    )
    damage_table_file = resolve_existing_path(damage_table_file, config_dir)
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

        # Simple one-ramp option:
        # annealRampTime=<s>, annealRampT1=<K>, annealRampT2=<K>
        # If provided, this overrides list-based annealStepDurations/annealTemperatures.
        simple_ramp_time = params.get("annealRampTime")
        simple_ramp_t1 = params.get(
            "annealRampT1",
            params.get("annealRampStartTemperature", params.get("annealRampStart")),
        )
        simple_ramp_t2 = params.get(
            "annealRampT2",
            params.get("annealRampEndTemperature", params.get("annealRampEnd")),
        )
        if (
            simple_ramp_time is not None
            and simple_ramp_t1 is not None
            and simple_ramp_t2 is not None
        ):
            anneal_step_durations = [float(simple_ramp_time)]
            parsed_schedule_temperatures = [float(simple_ramp_t1), float(simple_ramp_t2)]
            anneal_duration = float(simple_ramp_time)
            print(
                "Using simple anneal ramp from annealRampT1/annealRampT2/annealRampTime."
            )

        if parsed_schedule_temperatures:
            anneal_reference_temperature = parsed_schedule_temperatures[0]

        anneal_table_defaults = {}
        anneal_internal_defaults = {}
        if parse_bool(params.get("annealUseInternalDb", 1), True):
            anneal_db_file = params.get(
                "annealDbFile",
                "vsclib/anneal/si_score_defaults.json",
            )
            anneal_db_file = resolve_existing_path(anneal_db_file, config_dir)
            anneal_internal_defaults = load_internal_anneal_defaults(
                str(anneal_db_file),
                species=species,
                material=material,
                temperatureK=anneal_reference_temperature,
                length_unit_in_cm=length_unit_in_cm,
            )
            if anneal_internal_defaults:
                anneal_table_defaults.update(anneal_internal_defaults)
                print(
                    f"Using internal anneal defaults from {anneal_db_file}"
                )

        if parse_bool(params.get("annealUseTableLookup", 0), False):
            anneal_table_file = params.get(
                "annealTablePath",
                _vcs.getDefaultAnnealTablePath()
                if "_vcs" in globals() and hasattr(_vcs, "getDefaultAnnealTablePath")
                else "data/AnnealData/AdvCal_2023.12.fps",
            )
            anneal_table_file = resolve_existing_path(anneal_table_file, config_dir)
            anneal_table_defaults = load_advcal_anneal_defaults(
                str(anneal_table_file),
                material=material,
                species=species,
                temperatureK=anneal_reference_temperature,
                length_unit_in_cm=length_unit_in_cm,
            )
            if anneal_table_defaults:
                # Explicit AdvCal table can override internal database defaults.
                if anneal_internal_defaults:
                    merged_defaults = dict(anneal_internal_defaults)
                    merged_defaults.update(anneal_table_defaults)
                    anneal_table_defaults = merged_defaults
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

        # Prefer implant damage factors from internal vsclib anneal entry.
        for key in ("IFactor", "VFactor", "DFactor", "MCIFactor", "MCVFactor"):
            k = f"score{key}"
            if k in anneal_table_defaults:
                score_implant_factors[key] = float(anneal_table_defaults[k])
        if score_implant_factors:
            src = anneal_table_defaults.get("scoreFactorsSource", "vsclib")
            entry_id = anneal_table_defaults.get("entryId", "")
            if entry_id:
                print(f"Using internal SCORE factors from vsclib entry {entry_id} ({src})")
            else:
                print(f"Using internal SCORE factors from vsclib ({src})")
        elif parse_bool(params.get("annealUseScoreDirectFallback", 0), False):
            score_params_root = resolve_existing_path(
                params.get("scoreParamsRoot", "data/lib/score/Params"),
                config_dir,
            )
            score_implant_factors = load_score_implant_factors(
                str(score_params_root),
                material=material,
                species=species,
            )
            if score_implant_factors:
                print(
                    f"Using SCORE implant factors from {score_implant_factors.get('scoreInfoFile')}"
                )

        anneal = vcs.Anneal()
        anneal.setCellSet(cell_set)
        anneal.setSpeciesLabel(str(params.get("annealSpeciesLabel", implant_species_label)))
        anneal.setDuration(anneal_duration)
        anneal_mode = str(params.get("annealMode", "explicit")).strip().lower()
        if "_vcs" in globals() and hasattr(_vcs, "AnnealMode"):
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

        defect_coupling_enabled = parse_bool(params.get("annealDefectCoupling", 1.0), True)
        defect_use_equilibrium = parse_bool(params.get("annealDefectUseEquilibrium", 1), True)
        interstitial_label = str(params.get("annealInterstitialLabel", "Interstitial"))
        vacancy_label = str(params.get("annealVacancyLabel", "Vacancy"))
        ieq_value = params.get(
            "annealInterstitialEquilibrium",
            anneal_table_defaults.get("annealInterstitialEquilibrium", 0.0),
        )
        veq_value = params.get(
            "annealVacancyEquilibrium",
            anneal_table_defaults.get("annealVacancyEquilibrium", 0.0),
        )
        equilibrium_from_arrhenius = (
            "annealInterstitialEquilibrium" not in params
            and "annealVacancyEquilibrium" not in params
            and "annealInterstitialEquilibriumC0_cm3" in anneal_table_defaults
            and "annealVacancyEquilibriumC0_cm3" in anneal_table_defaults
        )

        if anneal_duration > 0.0 or anneal_step_durations:
            step_durations = anneal_step_durations
            step_temperatures = parsed_schedule_temperatures

            anneal.setDiffusionMaterials([1])
            anneal.setBlockingMaterials([2])
            if defect_coupling_enabled:
                defect_defaults = resolve_defect_coupling_defaults(
                    params=params,
                    anneal_table_defaults=anneal_table_defaults,
                    score_implant_factors=score_implant_factors,
                )
                anneal.enableDefectCoupling(True)
                anneal.setDamageLabels(
                    str(params.get("annealDamageLabel", "Damage")),
                    str(params.get("annealDamageLastImpLabel", "Damage_LastImp")),
                )
                anneal.setDefectLabels(
                    interstitial_label,
                    vacancy_label,
                )
                anneal.setDefectSourceWeights(
                    defect_defaults["annealDefectFromDamageHistoryWeight"],
                    defect_defaults["annealDefectFromDamageLastImpWeight"],
                )
                anneal.setDefectPartition(
                    defect_defaults["annealDefectToInterstitial"],
                    defect_defaults["annealDefectToVacancy"],
                )
                anneal.setDefectDiffusivities(
                    defect_defaults["annealInterstitialDiffusivity"],
                    defect_defaults["annealVacancyDiffusivity"],
                )
                anneal.setDefectReactionRates(
                    defect_defaults["annealDefectRecombinationRate"],
                    defect_defaults["annealInterstitialSinkRate"],
                    defect_defaults["annealVacancySinkRate"],
                )
                print(
                    "Defect reaction parameters: "
                    f"krec={defect_defaults['annealDefectRecombinationRate']:.3e}, "
                    f"kI={defect_defaults['annealInterstitialSinkRate']:.3e}, "
                    f"kV={defect_defaults['annealVacancySinkRate']:.3e}"
                )
                if "annealInterstitialSinkRateFromKsurfRaw" in anneal_table_defaults:
                    print(
                        "Ksurf-derived sink defaults (raw, bounded): "
                        f"kI_raw={anneal_table_defaults.get('annealInterstitialSinkRateFromKsurfRaw', 0.0):.3e}, "
                        f"kV_raw={anneal_table_defaults.get('annealVacancySinkRateFromKsurfRaw', 0.0):.3e}, "
                        f"cap={params.get('annealSinkRateCap', 1.0e-4):.3e}"
                    )
                if score_implant_factors:
                    print(
                        "Derived defect partition/TED from SCORE factors: "
                        f"IFactor={score_implant_factors.get('IFactor', float('nan')):.3g}, "
                        f"VFactor={score_implant_factors.get('VFactor', float('nan')):.3g}, "
                        f"DFactor={score_implant_factors.get('DFactor', float('nan')):.3g}, "
                        f"toI={defect_defaults['annealDefectToInterstitial']:.3f}, "
                        f"toV={defect_defaults['annealDefectToVacancy']:.3f}, "
                        f"TED={defect_defaults['annealTEDCoefficient']:.3g}"
                    )
                if defect_use_equilibrium:
                    anneal.enableDefectEquilibrium(True)
                    if equilibrium_from_arrhenius:
                        print(
                            "Defect equilibrium enabled with temperature-dependent Ieq/Veq "
                            "from Arrhenius Cstar."
                        )
                    else:
                        anneal.setDefectEquilibrium(
                            ieq_value,
                            veq_value,
                        )
                        print(
                            f"Defect equilibrium enabled: Ieq={ieq_value:.3e}, Veq={veq_value:.3e} "
                            "(simulation concentration units)."
                        )
                anneal.setDefectEnhancedDiffusion(
                    defect_defaults["annealTEDCoefficient"],
                    defect_defaults["annealTEDNormalization"],
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
                anneal.resetDefectInitialization()

            use_diagnostics = defect_coupling_enabled and parse_bool(
                params.get("annealDiagnostics", 1), True
            )
            diagnostics_rows = []

            # Execute anneal in stepwise mode to export defect diagnostics over time.
            # If no explicit schedule is given, use one isothermal segment.
            segments = []
            if step_durations:
                if step_temperatures and len(step_temperatures) == len(step_durations) + 1:
                    segments = [
                        {
                            "duration": dt,
                            "kind": "ramp",
                            "t0": step_temperatures[i],
                            "t1": step_temperatures[i + 1],
                            "teval": 0.5 * (step_temperatures[i] + step_temperatures[i + 1]),
                        }
                        for i, dt in enumerate(step_durations)
                    ]
                    print(
                        f"Using {len(step_durations)}-step temperature ramp anneal schedule."
                    )
                elif step_temperatures and len(step_temperatures) == len(step_durations):
                    segments = [
                        {
                            "duration": dt,
                            "kind": "iso",
                            "t0": t,
                            "t1": t,
                            "teval": t,
                        }
                        for dt, t in zip(step_durations, step_temperatures)
                    ]
                    print(
                        f"Using {len(step_durations)}-step isothermal anneal schedule."
                    )
                else:
                    fallback_temp = params.get("annealTemperature", 1273.15)
                    segments = [
                        {
                            "duration": dt,
                            "kind": "iso",
                            "t0": fallback_temp,
                            "t1": fallback_temp,
                            "teval": fallback_temp,
                        }
                        for dt in step_durations
                    ]
                    print(
                        "Using multi-step anneal schedule with shared annealTemperature."
                    )
            else:
                fallback_temp = params.get("annealTemperature", anneal_reference_temperature)
                segments = [
                    {
                        "duration": anneal_duration,
                        "kind": "iso",
                        "t0": fallback_temp,
                        "t1": fallback_temp,
                        "teval": fallback_temp,
                    }
                ]

            elapsed = 0.0
            for i, seg in enumerate(segments, start=1):
                seg_ieq = float(ieq_value)
                seg_veq = float(veq_value)
                if defect_coupling_enabled and defect_use_equilibrium:
                    seg_ieq, seg_veq = resolve_equilibrium_at_temperature(
                        params=params,
                        anneal_table_defaults=anneal_table_defaults,
                        temperatureK=float(seg["teval"]),
                        fallback_ieq=seg_ieq,
                        fallback_veq=seg_veq,
                    )
                    anneal.setDefectEquilibrium(seg_ieq, seg_veq)

                if hasattr(anneal, "clearTemperatureSchedule"):
                    anneal.clearTemperatureSchedule()
                    if seg["kind"] == "ramp":
                        anneal.addRampStep(seg["duration"], seg["t0"], seg["t1"])
                    else:
                        anneal.addIsothermalStep(seg["duration"], seg["t0"])
                    anneal.setDuration(0.0)
                else:
                    anneal.setDuration(seg["duration"])
                    anneal.setTemperature(seg["t0"])

                anneal.apply()
                elapsed += seg["duration"]
                if use_diagnostics:
                    row = collect_defect_diagnostics(
                        cell_set,
                        material_label="Material",
                        substrate_material_id=1,
                        interstitial_label=interstitial_label,
                        vacancy_label=vacancy_label,
                        ieq=seg_ieq,
                        veq=seg_veq,
                        time_s=elapsed,
                        step_index=i,
                        temperatureK=seg["teval"],
                    )
                    if row is not None:
                        diagnostics_rows.append(row)
                print(f"Anneal step {i}/{len(segments)} completed.")

            if use_diagnostics and diagnostics_rows:
                csv_file = Path(str(params.get("annealDiagnosticsCsv", "anneal_defect_diagnostics.csv")))
                if not csv_file.is_absolute():
                    csv_file = (config_dir / csv_file).resolve()
                write_diagnostics_csv(csv_file, diagnostics_rows)
                print(f"Wrote anneal defect diagnostics: {csv_file}")

                if parse_bool(params.get("annealDiagnosticsPlot", 1), True):
                    plot_file = Path(
                        str(params.get("annealDiagnosticsPlotFile", "anneal_defect_diagnostics.png"))
                    )
                    if not plot_file.is_absolute():
                        plot_file = (config_dir / plot_file).resolve()
                    write_diagnostics_plot(plot_file, diagnostics_rows)
                    if plot_file.exists():
                        print(f"Wrote anneal defect diagnostics plot: {plot_file}")

    # Store annealed profile explicitly as a separate scalar.
    annealed = list(cell_set.getScalarData(implant_species_label))
    cell_set.setScalarData("concentration_annealed", annealed)

    if (
        parse_bool(params.get("annealDefectCoupling", 1.0), True)
        and parse_bool(params.get("annealDefectUseEquilibrium", 1), True)
    ):
        interstitial = cell_set.getScalarData(str(params.get("annealInterstitialLabel", "Interstitial")))
        vacancy = cell_set.getScalarData(str(params.get("annealVacancyLabel", "Vacancy")))
        ieq_value = params.get(
            "annealInterstitialEquilibrium",
            anneal_table_defaults.get("annealInterstitialEquilibrium", 0.0)
            if "anneal_table_defaults" in locals()
            else 0.0,
        )
        veq_value = params.get(
            "annealVacancyEquilibrium",
            anneal_table_defaults.get("annealVacancyEquilibrium", 0.0)
            if "anneal_table_defaults" in locals()
            else 0.0,
        )
        if interstitial is not None and vacancy is not None and ieq_value > 0.0 and veq_value > 0.0:
            i_over_ieq = [max(float(v), 0.0) / ieq_value for v in interstitial]
            v_over_veq = [max(float(v), 0.0) / veq_value for v in vacancy]
            iv_over_eq = [
                (max(float(i), 0.0) * max(float(v), 0.0)) / (ieq_value * veq_value)
                for i, v in zip(interstitial, vacancy)
            ]
            cell_set.addScalarData("I_over_Ieq", 0.0)
            cell_set.addScalarData("V_over_Veq", 0.0)
            cell_set.addScalarData("IV_over_IeqVeq", 0.0)
            cell_set.setScalarData("I_over_Ieq", i_over_ieq)
            cell_set.setScalarData("V_over_Veq", v_over_veq)
            cell_set.setScalarData("IV_over_IeqVeq", iv_over_eq)

    cell_set.writeVTU("final.vtu")
    print("Wrote initial.vtu and final.vtu")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
