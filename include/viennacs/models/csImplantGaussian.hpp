#pragma once

#include "../csImplantModel.hpp"

#include <algorithm>
#include <cmath>

namespace viennacs {

template <class NumericType, int D>
class ImplantGaussian final : public ImplantModel<NumericType, D> {
public:
  ImplantGaussian(NumericType mu, NumericType sigma, NumericType lateralSigma,
                  NumericType lateralMu)
      : mu_(mu), sigma_(sigma), lateralMu_(lateralMu),
        lateralSigma_(lateralSigma) {}

  NumericType getDepthProfile(NumericType depth) override {
    return (1.0 / (sigma_ * std::sqrt(2 * M_PI))) *
           std::exp(-0.5 * std::pow((depth - mu_) / sigma_, 2));
  }

  NumericType getLateralProfile(NumericType offset,
                                NumericType depth) override {
    // profile of empirical implant model in lateral directions
    return (1.0 / (lateralSigma_ * std::sqrt(2 * M_PI))) *
           std::exp(-0.5 * std::pow((offset - lateralMu_) / lateralSigma_, 2));
  }

  NumericType getMaxDepth() override { return std::max(mu_ + 6 * sigma_, mu_); }

  NumericType getMaxLateralRange() override {
    return std::abs(lateralMu_) + 6 * lateralSigma_;
  }

private:
  const NumericType mu_;
  const NumericType sigma_;

  const NumericType lateralMu_;
  const NumericType lateralSigma_;
};
} // namespace viennacs
