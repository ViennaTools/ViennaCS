# implantMask2D: Masked Implant + Damage + Anneal (ViennaCS Core)

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

When `annealDefectCoupling=1`, interstitial (`I`) and vacancy (`V`) fields are initialized from damage:

```text
S  = wh * Damage + wl * Damage_LastImp
I0 = fI * S
V0 = fV * S
```

Then per timestep:

```text
dI/dt = DI * Lap(I) - krec * I * V - kI * I
dV/dt = DV * Lap(V) - krec * I * V - kV * V
```

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
  - Current lookup maps AdvCal `Si B Int D`, `Si Int Di`, `Si Vac Dv`, bulk I-V recombination surrogate, and Si ICluster (`Ikfi`, `Ikfc`, `Ikr`, `InitPercent`)
- Temperature schedule keys: `annealStepDurations`, `annealTemperatures`, `annealRampTemperatures` (alias fallback in script)
- Defect coupling enable: `annealDefectCoupling`
- Defect labels: `annealDamageLabel`, `annealDamageLastImpLabel`, `annealInterstitialLabel`, `annealVacancyLabel`
- Defect source and partition: `annealDefectFromDamageHistoryWeight`, `annealDefectFromDamageLastImpWeight`, `annealDefectToInterstitial`, `annealDefectToVacancy`
- Defect transport and reactions: `annealInterstitialDiffusivity`, `annealVacancyDiffusivity`, `annealDefectRecombinationRate`, `annealInterstitialSinkRate`, `annealVacancySinkRate`
- Interstitial clustering (one-moment surrogate): `annealDefectClustering`, `annealIClusterLabel`, `annealIClusterIkfi`, `annealIClusterIkfc`, `annealIClusterIkr`, `annealIClusterInitFraction`
  - This mapping is experimental in the reduced ViennaCS defect model and is off by default.
- TED coupling: `annealTEDCoefficient`, `annealTEDNormalization`

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

Units:
- geometric lengths in this example are nm (because `setLengthUnitInCm(1e-7)`)
- concentration and damage outputs are scaled to `cm^-3` when `setOutputConcentrationInCm3(True)` is active

## 10. Practical notes

- For robust TCAD calibration, prefer table-driven mode (`config_table.txt`) and keep manual moment mode for controlled experiments.
- For large anneal timesteps, use `annealMode=implicit`.
- Defect kinetics are stiff in `cm^-3` units. If `annealDefectRecombinationRate`, sink rates, or `annealDefectToInterstitial/Vacancy` are too large, `Interstitial/Vacancy` can collapse in high-damage regions and produce hollow/ring-like artifacts.
- The AdvCal bulk I-V recombination term is much stiffer than this reduced model can use directly; the example applies a bounded surrogate default and prints both bounded and raw values.
