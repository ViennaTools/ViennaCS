#pragma once

// csSheetResistance — computes the sheet resistance (Rsh) of a doped
// semiconductor layer from a 1-D concentration profile stored in a
// DenseCellSet scalar field.
//
// Usage (C++)
// -----------
//   auto sr = SheetResistance<double, 2>{};
//   sr.setCellSet(cellSet);
//   sr.setConcentrationLabel("P_active");  // nm⁻³ field from Anneal
//   double rsh = sr.computeElectron();     // Ω/□, Masetti n-type model
//
// Usage (Python, ViennaCS or ViennaPS module)
// --------------------------------------------
//   sr = vps.SheetResistance()
//   sr.setCellSet(domain.getCellSet())
//   sr.setConcentrationLabel("P_active")
//   rsh = sr.computeElectron()

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <csDenseCellSet.hpp>

namespace viennacs {

// ── Masetti-Severi mobility models
// ───────────────────────────────────────────── Masetti et al., IEEE Trans.
// Electron Devices 30, 764 (1983) — Table I.
//
// General form:
//   μ = μ_min * exp(-Pc/N) + (μ_max - μ_min)/(1 + (N/Cr)^α) + μ1/(1+(Cs/N)^β)
//
// N is the majority-carrier impurity concentration in cm⁻³.

template <typename T> struct MasettiElectronMobility {
  // Electrons in phosphorus-doped silicon.
  T operator()(T N_cm3) const noexcept {
    const T N = std::max(N_cm3, T(1e10));
    return (T(68.5) * std::exp(-T(9.23e16) / N) +
            T(1345.5) / (T(1) + std::pow(N / T(9.68e16), T(0.680))) +
            T(56.1) / (T(1) + std::pow(T(3.43e20) / N, T(2.0))));
  }
  // where 1345.5 = μ_max(1414) − μ_min(68.5)
};

template <typename T> struct MasettiHoleMobility {
  // Holes in boron-doped silicon.
  T operator()(T N_cm3) const noexcept {
    const T N = std::max(N_cm3, T(1e10));
    return (T(44.9) * std::exp(-T(2.23e17) / N) +
            T(425.6) / (T(1) + std::pow(N / T(2.35e17), T(0.719))) +
            T(29.0) / (T(1) + std::pow(T(1.80e20) / N, T(2.0))));
  }
  // where 425.6 = μ_max(470.5) − μ_min(44.9)
};

// ── SheetResistance
// ────────────────────────────────────────────────────────────
//
// Computes Rsh [Ω/□] by extracting a 1-D depth profile from a DenseCellSet
// scalar field and integrating:
//
//   Rsh = 1 / (q · ∫ μ(C_active(z)) · C_active(z) dz_cm)
//
// Cells are binned by rounded depth coordinate.  Lateral cells at the same
// depth are averaged (appropriate for blanket / 1-D profiles).
//
// Default configuration targets ViennaPS domains (length unit = 1 nm):
//   • lengthUnit  = 1e-7  (nm → cm)
//   • concUnit    = 1e21  (nm⁻³ → cm⁻³)
//   • depthAxis   = D−1   (y for 2-D, z for 3-D)
//   • surfacePosition = 0  (substrate at negative coordinates; depth = surface
//   − y)

template <typename NumericType, int D> class SheetResistance {
  using CellSet = DenseCellSet<NumericType, D>;

  SmartPointer<CellSet> cellSet_;
  std::string concentrationLabel_ = "active_concentration";
  int depthAxis_ = D - 1;
  NumericType surfacePosition_ = NumericType(0);
  NumericType lengthUnitCm_ = NumericType(1e-7); // nm → cm
  NumericType concUnitCm3_ = NumericType(1e21);  // nm⁻³ → cm⁻³

public:
  SheetResistance() = default;
  explicit SheetResistance(SmartPointer<CellSet> cs)
      : cellSet_(std::move(cs)) {}

  // ── configuration ─────────────────────────────────────────────────────────

  void setCellSet(SmartPointer<CellSet> cs) { cellSet_ = std::move(cs); }

  // Cell-set scalar field containing the active concentration (default:
  // "active_concentration" — the field written by the Anneal solid-activation
  // model; for total concentration use the total/species label directly).
  void setConcentrationLabel(const std::string &label) {
    concentrationLabel_ = label;
  }

  // Index of the depth axis in the cell-centre coordinate vector.
  // Default: D−1  (y for 2-D, z for 3-D).
  void setDepthAxis(int axis) { depthAxis_ = axis; }

  // Coordinate of the wafer surface along depthAxis_. Depth is computed as
  // surfacePosition - coordinate, matching SIMS' positive-into-substrate
  // convention for ViennaPS domains with the substrate at y < 0.
  void setSurfacePosition(NumericType surfacePosition) {
    surfacePosition_ = surfacePosition;
  }

  // Conversion factor from the domain's length unit to cm.
  // Setting this also updates concUnit to be consistent (1/lu³).
  // Default: 1e-7  (1 nm = 1e-7 cm).
  void setLengthUnit(NumericType lu_cm) {
    lengthUnitCm_ = lu_cm;
    const NumericType lu3 = lu_cm * lu_cm * lu_cm;
    concUnitCm3_ =
        (lu3 > NumericType(0)) ? (NumericType(1) / lu3) : NumericType(1e21);
  }

  // Override the concentration unit independently.
  // Default: 1e21  (1 nm⁻³ = 1e21 cm⁻³).
  void setConcentrationUnit(NumericType unit) { concUnitCm3_ = unit; }

  // ── computation ───────────────────────────────────────────────────────────

  // Compute Rsh [Ω/□] using the supplied mobility functor.
  // The functor must accept a concentration in cm⁻³ and return mobility in
  // cm²/(V·s).
  template <typename MobilityFunctor>
  NumericType compute(MobilityFunctor mobility) const {
    if (!cellSet_)
      return std::numeric_limits<NumericType>::infinity();

    const auto *conc_data = cellSet_->getScalarData(concentrationLabel_);
    if (!conc_data || conc_data->empty())
      return std::numeric_limits<NumericType>::infinity();

    const NumericType delta =
        static_cast<NumericType>(cellSet_->getGridDelta());
    const int n = static_cast<int>(cellSet_->getNumberOfCells());

    // Bin concentrations by rounded depth coordinate.
    // std::map keeps bins sorted by depth automatically.
    std::map<NumericType, std::pair<NumericType, int>> bins;
    for (int i = 0; i < n; ++i) {
      const NumericType C = (*conc_data)[i];
      if (C <= NumericType(0))
        continue;
      const auto center = cellSet_->getCellCenter(i);
      const NumericType z =
          surfacePosition_ - static_cast<NumericType>(center[depthAxis_]);
      if (z < NumericType(0))
        continue;
      const NumericType key =
          std::round(static_cast<double>(z / delta)) * delta;
      auto &slot = bins[key];
      slot.first += C;
      slot.second += 1;
    }

    if (bins.size() < 2)
      return std::numeric_limits<NumericType>::infinity();

    // Build depth/concentration arrays (average across lateral direction).
    std::vector<NumericType> depths, concs;
    depths.reserve(bins.size());
    concs.reserve(bins.size());
    for (const auto &[depth, val] : bins) {
      depths.push_back(depth);
      concs.push_back(val.first / static_cast<NumericType>(val.second));
    }

    // Trapezoidal integration: q · ∫ μ(C_cm3) · C_cm3 dz_cm
    constexpr NumericType q = NumericType(1.602176634e-19);
    NumericType conductance = NumericType(0);
    for (std::size_t i = 0; i + 1 < depths.size(); ++i) {
      const NumericType dz_cm = (depths[i + 1] - depths[i]) * lengthUnitCm_;
      const NumericType C0 = concs[i] * concUnitCm3_;
      const NumericType C1 = concs[i + 1] * concUnitCm3_;
      conductance +=
          NumericType(0.5) * dz_cm * (mobility(C0) * C0 + mobility(C1) * C1);
    }
    conductance *= q;

    return (conductance > NumericType(0))
               ? NumericType(1) / conductance
               : std::numeric_limits<NumericType>::infinity();
  }

  // Convenience: Masetti electron mobility (n-type, e.g. P-doped Si).
  NumericType computeElectron() const {
    return compute(MasettiElectronMobility<NumericType>{});
  }

  // Convenience: Masetti hole mobility (p-type, e.g. B-doped Si).
  NumericType computeHole() const {
    return compute(MasettiHoleMobility<NumericType>{});
  }
};

} // namespace viennacs
