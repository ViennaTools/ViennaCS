#pragma once

// csScreenEnergyLoss — analytic screen-oxide energy-loss (range-reduction) model.
//
// An ion that passes through a screen layer (e.g. a screen oxide) before
// entering the substrate loses energy, so its projected range in the substrate
// is reduced. This model returns the range-scale factor
//     k = Rp_substrate(E1(t)) / Rp_substrate(E1(ref))
// relative to a calibrated reference screen thickness, so an implant profile's
// projected range mu can be scaled mu -> k*mu (the straggle is left unchanged,
// which keeps the Pearson-IV shape valid).
//
// Power-law range-energy model, Rp(E) = a * E^p (a few constants per species,
// fitted once to tabulated Rp(E) data). Only the ratio to the reference is used,
// so the substrate prefactor a_substrate cancels and k == 1 exactly at the
// reference thickness (the calibrated shapes are untouched at the reference).
//
//   R0    = a_screen * E0^p_screen                  range in screen at incoming E0
//   E1(t) = ((R0 - t) / a_screen)^(1/p_screen)      residual energy after t of screen
//   k(t)  = ( E1(t) / E1(ref) )^p_substrate

#include <algorithm>
#include <cmath>

namespace viennacs {

template <class NumericType> class ScreenEnergyLoss {
public:
  ScreenEnergyLoss(NumericType screenRangePrefactor,
                   NumericType screenRangeExponent,
                   NumericType substrateRangeExponent, NumericType energyKeV,
                   NumericType referenceThicknessNm)
      : aScreen_(screenRangePrefactor), pScreen_(screenRangeExponent),
        pSubstrate_(substrateRangeExponent), energy_(energyKeV),
        refThickness_(referenceThicknessNm) {}

  // Residual ion energy [keV] after passing through `thickness` nm of screen.
  NumericType residualEnergy(NumericType thickness) const {
    const NumericType r0 = aScreen_ * std::pow(energy_, pScreen_);
    const NumericType residualRange =
        std::max(r0 - thickness, NumericType(1e-6));
    return std::pow(residualRange / aScreen_, NumericType(1) / pScreen_);
  }

  // Projected-range scale factor k for `thickness` nm of screen, relative to the
  // calibrated reference. k == 1 exactly at the reference thickness.
  NumericType rangeScale(NumericType thickness) const {
    if (std::abs(thickness - refThickness_) < NumericType(1e-9))
      return NumericType(1);
    const NumericType num = std::pow(residualEnergy(thickness), pSubstrate_);
    const NumericType den = std::pow(residualEnergy(refThickness_), pSubstrate_);
    return den > NumericType(0) ? num / den : NumericType(1);
  }

  NumericType getReferenceThickness() const { return refThickness_; }

private:
  NumericType aScreen_;      // screen range prefactor   [nm / keV^p]
  NumericType pScreen_;      // screen range exponent
  NumericType pSubstrate_;   // substrate range exponent
  NumericType energy_;       // incoming ion energy      [keV]
  NumericType refThickness_; // calibrated reference screen thickness [nm]
};

} // namespace viennacs
