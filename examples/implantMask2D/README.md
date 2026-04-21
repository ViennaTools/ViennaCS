# implantMask2D: Masked Implant + Damage + Anneal

This example builds a 2D masked structure, performs implantation with optional table-driven calibration, computes implant damage, and optionally anneals dopant and defects (explicit or implicit, isothermal or ramped multi-step).

The script is [implantMask2D.py](/home/filipov/Software/IonImplantation/ViennaCS/examples/implantMask2D/implantMask2D.py).  
Core physics lives in:
- [csImplant.hpp](/home/filipov/Software/IonImplantation/ViennaCS/include/viennacs/csImplant.hpp)
- [csImplantPearson.hpp](/home/filipov/Software/IonImplantation/ViennaCS/include/viennacs/models/csImplantPearson.hpp)
- [csImplantDamage.hpp](/home/filipov/Software/IonImplantation/ViennaCS/include/viennacs/models/csImplantDamage.hpp)
- [csImplantTable.hpp](/home/filipov/Software/IonImplantation/ViennaCS/include/viennacs/models/csImplantTable.hpp)
- [csAnneal.hpp](/home/filipov/Software/IonImplantation/ViennaCS/include/viennacs/csAnneal.hpp)

## 1. How to run

From repository root:

```bash
python -m pip install .
cd examples/implantMask2D
python implantMask2D.py config_table.txt
```

Outputs:
- `initial.vtu`: geometry and initial fields
- `final.vtu`: implanted + optional annealed fields

## 2. Structure and materials

The example creates:
- Material `1`: Si substrate
- Material `2`: mask
- Material `3`: optional SiO2 screen/cap (`oxideThickness > 0`)
- Material `0`: cover/gas

Mask opening is rectangular (`openingWidth`), with top-space above the structure.

## 3. Implant model equations

ViennaCS applies an empirical kernel:

```text
P(d, l) = Pd(d) * Pl(l, d)
```

where:
- `d`: depth along beam direction
- `l`: lateral displacement

Per beam/cell contribution:

```text
DeltaC = w_dose * P(d, l)
```

with dose weight:
- `Off`: `w_dose = 1`
- `WaferDose`: `w_dose = D_wafer * uL^2 * DeltaS`
- `BeamDose`: `w_dose = D_wafer * cos(theta) * uL^2 * DeltaS`

Here `uL = lengthUnitInCm` (example sets `1e-7` cm per nm), and `DeltaS` is beam area/length sampling from `gridDelta`.

### 3.1 Depth profile (Pearson IV)

`ImplantPearsonIV` uses `PearsonIV(mu, sigma, beta, gamma)` from [csConstants.hpp](/home/filipov/Software/IonImplantation/ViennaCS/include/viennacs/csConstants.hpp), normalized numerically over depth:

```text
Pd(d) = PearsonIV(d; mu, sigma, beta, gamma)
        / Integral[0..dmax] PearsonIV(x) dx
```

### 3.2 Crystalline tail options

Two tail-capable options exist:
- `ImplantPearsonIVChanneling`: Pearson + exponential channeling tail with smooth onset
```text
Pd = (1 - ft) * Ppearson + ft * Ptail
Ptail ~ s(d) * exp(-(d - d0) / lambda_t)
```
- `ImplantDualPearsonIV`: weighted sum of head/tail Pearson profiles
```text
Pd = fh * Phead + (1 - fh) * Ptail
```

Table-driven crystalline mode additionally applies channeling suppression scaling vs tilt/screen/damage in [csImplantTable.hpp](/home/filipov/Software/IonImplantation/ViennaCS/include/viennacs/models/csImplantTable.hpp).

### 3.3 Lateral profile

Base lateral kernel is Gaussian with depth-dependent `sigma_l(d)`:

```text
Pl(l, d) = 1 / (sigma_l(d) * sqrt(2*pi))
           * exp(-0.5 * ((l - mu_l) / sigma_l(d))^2)
```

Supported `sigma_l(d)` modes:
- `Constant`: `sigma_l = sigma0`
- `SentaurusS3`: `sigma_l = sigma0 * exp(-d / (sigma0 * lv))`
- `Taurus`: `sigma_l = sigma0 * (1 + DeltaSigma * (d / Rp - 1))`
- `Dios`: `p1..p5`-based log-sum-exp form (see code)

## 4. Damage model equations

Damage uses `ImplantDamageHobler`:
- Gaussian core around `Rp` with `sigma`
- Piecewise exponential branch controlled by `lambda`
- `lambda > 0`: near-surface exponential + Gaussian deeper
- `lambda < 0`: Gaussian near `Rp` + exponential bulk tail

Raw depth shape is normalized and scaled by `defectsPerIon`:

```text
Pd_damage(d) = Ndef_per_ion * SHobler(d)
               / Integral[0..dmax] SHobler(x) dx
```

Lateral damage spread uses the same lateral framework as dopant.

Accumulation in core:
- `Damage_LastImp`: damage from this implant call
- `Damage`: cumulative field, updated as

```text
Damage <- Damage + damageFactor * Damage_LastImp
```

## 5. Visibility / mask / screen behavior

For each beam:
- gas (`Material=0`) is skipped
- mask materials terminate that beam path
- screen materials advance the beam without depositing dose in that cell
- deposition starts at first non-gas, non-mask, non-screen material

`beamHits` is a kernel-footprint indicator (set to 1 where near-center contribution occurs), not a strict binary line-of-sight map.

## 6. Anneal model equations

Anneal solves diffusion on the DenseCellSet neighborhood graph.

### 6.1 Dopant diffusion

```text
dC/dt = div(D * grad(C))
```

with either:
- constant `D`: `annealDiffusionCoefficient`
- Arrhenius:

```text
D(T) = D0 * exp(-Ea / (kB * T))
kB   = 8.617333262145e-5 eV/K
```

Numerics:
- `explicit`: forward update with stability-limited time step
- `implicit`: backward-Euler solved by Gauss-Seidel iterations on sparse stencil

### 6.2 Defect-coupled anneal (I/V + TED)

With defect-coupled anneal enabled (default: `annealDefectCoupling=1`), interstitial (`I`) and vacancy (`V`) fields are initialized from damage:

```text
S  = wh * Damage + wl * Damage_LastImp
I0 = fI * S
V0 = fV * S
```

Then per timestep:

```text
dI/dt = DI * Lap(I) - krec * (I*V - Ieq*Veq) - kI * (I - Ieq)
dV/dt = DV * Lap(V) - krec * (I*V - Ieq*Veq) - kV * (V - Veq)
```

With `annealDefectUseEquilibrium=1` (default), equilibrium values can be taken from AdvCal `Cstar` lookup or provided explicitly.  
With `annealDefectUseEquilibrium=0`, `Ieq=Veq=0` and fields represent excess-defect-like concentrations.

TED coupling modifies dopant diffusivity locally:

```text
Deff = D(T) * (1 + kTED * max(0, (I - V) / NTED))
```

### 6.3 Temperature ramps and multi-step anneal

Core anneal supports schedules:
- isothermal steps: `(duration_i, Ti)`
- ramp steps: `(duration_i, Ti -> T(i+1))`

Config encoding:
- `annealStepDurations=t1,t2,...,tN`
- isothermal: `annealTemperatures=T1,...,TN`
- ramp: `annealTemperatures=T0,T1,...,TN` (N+1 temperatures)

If `annealDuration` is omitted, script defaults it to `sum(annealStepDurations)`.

## 7. Table lookup behavior

When manual Pearson/damage moments are omitted, recipe-driven models read tables (for example Taurus under `data/ImplData`).

Selection process:
1. Filter by species/material/substrate/model compatibility.
2. Restrict to local grid neighbors in each axis (`energyKeV`, `tiltDeg`, `rotationDeg`, `dosePerCm2`, `screenThickness`).
3. If multiple remain, inverse-distance blending:

```text
wi = 1 / distance_i
p  = (sum_i wi * pi) / (sum_i wi)
```

For crystalline dual-Pearson entries, additional reshaping/tail suppression is applied from `tilt`, `screenThickness`, and `damageLevel`.

## 8. Configuration reference (script-level)

### Geometry/process
- `gridDelta`, `xExtent`, `topSpace`, `substrateDepth`
- `openingWidth`, `maskHeight`, `oxideThickness`
- `angle`, `rotationDeg`, `doseCm2`, `doseControl`

### Implant recipe / table
- `species`, `material`, `energyKeV`, `substrateType`, `screenThickness`, `damageLevel`
- `preferredModel` (`PearsonIV`, `DualPearsonIV`, or `auto`)
- `tablePath`
- Internal implant DB: `implantUseInternalDb`, `implantDbFile` (default: `vsclib/implant/si_taurus_defaults.json`)
- `implantSpeciesLabel` (default: `concentration_annealed`)

### Manual implant moments (used when `projectedRange` is provided)
- Head: `projectedRange`, `depthSigma`, `skewness`, `kurtosis`, `lateralMu`, `lateralSigma`
- Tail: `headFraction`, `tailProjectedRange`, `tailDepthSigma`, `tailSkewness`, `tailKurtosis`, `tailLateralMu`, `tailLateralSigma`

### Damage recipe / table
- `damageTablePath`
- Manual damage fields: `damageProjectedRange`, `damageVerticalSigma`, `damageLambda`, `damageDefectsPerIon`, `damageLateralMu`, `damageLateralSigma`, `damageLateralModel`, `damageLateralScale`, `damageLateralLv`, `damageLateralDeltaSigma`, `damageLateralP1..P5`
- `damageFactor`

### Anneal (core)
- Activation: `annealDuration` or `annealStepDurations`
- Anneal species field: `annealSpeciesLabel` (default: `implantSpeciesLabel`)
- Mode: `annealMode` (`explicit`/`implicit`), `annealTimeStep`, `annealStabilityFactor`
- Implicit solver: `annealImplicitMaxIterations`, `annealImplicitTolerance`
- Dopant diffusion (constant): `annealDiffusionCoefficient`
- Dopant diffusion (Arrhenius): `annealD0`, `annealEa`, `annealTemperature`
- Anneal table lookup: `annealUseTableLookup`, `annealTablePath` (default: `data/AnnealData/AdvCal_2023.12.fps` in the package)
  - Current lookup maps AdvCal `Si B Int D`, `Si Int Di`, `Si Vac Dv`, `Si Int Cstar`, `Si Vac Cstar`, bulk I-V recombination surrogate, and Si ICluster (`Ikfi`, `Ikfc`, `Ikr`, `InitPercent`)
- Internal anneal DB: `annealUseInternalDb`, `annealDbFile` (default: `vsclib/anneal/si_score_defaults.json`)
- SCORE implant factor fallback (off by default): `annealUseScoreDirectFallback=1`, `scoreParamsRoot` (default: `data/lib/score/Params`)
- Temperature schedule keys: `annealStepDurations`, `annealTemperatures`, `annealRampTemperatures` (alias fallback in script)
- Defect coupling enable: `annealDefectCoupling` (default: `1`; set `0` to disable)
- Defect labels: `annealDamageLabel`, `annealDamageLastImpLabel`, `annealInterstitialLabel`, `annealVacancyLabel`
- Defect source and partition: `annealDefectFromDamageHistoryWeight`, `annealDefectFromDamageLastImpWeight`, `annealDefectToInterstitial`, `annealDefectToVacancy`
- Defect transport and reactions: `annealInterstitialDiffusivity`, `annealVacancyDiffusivity`, `annealDefectRecombinationRate`, `annealInterstitialSinkRate`, `annealVacancySinkRate`
- Sink surrogate controls: `annealSinkRateCap`, `annealKsurfRecombinationLengthCm`
- Defect equilibrium control: `annealDefectUseEquilibrium`, `annealInterstitialEquilibrium`, `annealVacancyEquilibrium`
- Interstitial clustering (one-moment surrogate): `annealDefectClustering`, `annealIClusterLabel`, `annealIClusterIkfi`, `annealIClusterIkfc`, `annealIClusterIkr`, `annealIClusterInitFraction`
  - This mapping is experimental in the reduced ViennaCS defect model and is off by default.
- TED coupling: `annealTEDCoefficient`, `annealTEDNormalization`
- Diagnostics: `annealDiagnostics`, `annealDiagnosticsCsv`, `annealDiagnosticsPlot`, `annealDiagnosticsPlotFile`

### 8.1 Defect parameter meaning, origin, and dopant dependence

- `annealDefectCoupling`: turns on solved I/V fields and TED coupling (default `1`); if set to `0`, anneal is pure dopant diffusion.
- `annealDefectFromDamageHistoryWeight`: weight of accumulated `Damage` when seeding initial I/V fields at anneal start.
- `annealDefectFromDamageLastImpWeight`: weight of `Damage_LastImp` when seeding I/V; in single-implant examples this is usually the dominant source.
- `annealDefectToInterstitial`: partition factor from defect source to initial interstitial field.
- `annealDefectToVacancy`: partition factor from defect source to initial vacancy field.
- `annealInterstitialSinkRate`: linear relaxation/removal coefficient for interstitials; reduced-model surrogate, not directly looked up from AdvCal.
- `annealVacancySinkRate`: linear relaxation/removal coefficient for vacancies; reduced-model surrogate, not directly looked up from AdvCal.
- `annealDefectUseEquilibrium`: switches recombination/sinks from decay-to-zero (excess-defect mode) to relaxation toward `(Ieq,Veq)`.
- `annealTEDCoefficient`: scales how strongly `(I-V)` boosts dopant diffusivity in this reduced model.
- `annealTEDNormalization`: reference concentration scale used to nondimensionalize `(I-V)` inside TED factor.
- `annealDefectClustering`: enables one-moment interstitial-cluster surrogate (`ICluster`) in the reduced model.

Origin and calibration notes:
- `annealInterstitialDiffusivity`, `annealVacancyDiffusivity`, `annealD0`, `annealEa`, `annealInterstitialEquilibrium`, `annealVacancyEquilibrium` are looked up from `AdvCal_2023.12.fps` (or overridden explicitly).
- `annealDefectRecombinationRate` defaults from an AdvCal-inspired capture-radius formula, then is bounded for stability in this reduced PDE model.
- `annealIClusterIkfi`, `annealIClusterIkfc`, `annealIClusterIkr`, `annealIClusterInitFraction` can be mapped from AdvCal but remain experimental here.
- If not explicitly provided, `annealDefectToInterstitial`/`annealDefectToVacancy` are now derived from SCORE implant factors (`MCIFactor`/`MCVFactor`, fallback `IFactor`/`VFactor`) stored in internal `vsclib` anneal entries.
- Direct parsing from `data/lib/score/Params/.../Info` is available only via `annealUseScoreDirectFallback=1`.
- If not explicitly provided, sink rates are derived from Ksurf-style defaults via `k_sink ~ D/Lrec^2` and then bounded by `annealSinkRateCap` (default `1e-4`).
- If not explicitly provided, `annealTEDCoefficient` defaults to `0.5 * DFactor` from SCORE `Info`; `annealTEDNormalization` remains user-controlled (default `1e20`).

Dopant dependence:
- Physically, most of these parameters are not universal constants; they depend on dopant species, damage state, and anneal temperature.
- In this example, only parameters explicitly table-mapped from AdvCal are species/material-aware by default; the remaining knobs are user calibration parameters.

## 9. Output fields in `final.vtu`

Common fields:
- `concentration_as_implanted`
- `concentration_annealed`
- `beamHits`
- `Damage`
- `Damage_LastImp`

If defect-coupled anneal is enabled:
- `Interstitial`
- `Vacancy`
- `ICluster` (when defect clustering is enabled and kinetics are nonzero)

If defect equilibrium mode is enabled:
- `I_over_Ieq`
- `V_over_Veq`
- `IV_over_IeqVeq`

Units:
- geometric lengths in this example are nm (because `setLengthUnitInCm(1e-7)`)
- concentration and damage outputs are scaled to `cm^-3` when `setOutputConcentrationInCm3(True)` is active

## 10. Practical notes

- For robust TCAD calibration, prefer table-driven mode (`config_table.txt`) and keep manual moment mode for controlled experiments.
- For large anneal timesteps, use `annealMode=implicit`.
- Defect kinetics are stiff in `cm^-3` units. If `annealDefectRecombinationRate`, sink rates, or `annealDefectToInterstitial/Vacancy` are too large, `Interstitial/Vacancy` can collapse quickly; this is expected in excess-defect mode, but may be too aggressive for long thermal budgets.
- The AdvCal bulk I-V recombination term is much stiffer than this reduced model can use directly; the example applies a bounded surrogate default and prints both bounded and raw values.
