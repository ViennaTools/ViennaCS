#pragma once

#include "../csConstants.hpp"
#include "../csImplantModel.hpp"

#include <cmath>

namespace viennacs {

template <class NumericType, int D>
class ImplantPearsonIV final : public ImplantModel<NumericType, D> {
public:
  ImplantPearsonIV(const constants::PearsonIVParameters<NumericType> &params,
                   NumericType lateralMu, NumericType lateralSigma)
      : params_(params), lateralMu_(lateralMu), lateralSigma_(lateralSigma) {}

  NumericType getDepthProfile(NumericType depth) override {
    return constants::PearsonIV(depth, params_);
  }

  NumericType getLateralProfile(NumericType offset,
                                NumericType depth) override {
    return (1.0 / (lateralSigma_ * std::sqrt(2 * M_PI))) *
           std::exp(-0.5 * std::pow((offset - lateralMu_) / lateralSigma_, 2));
  }

private:
  const constants::PearsonIVParameters<NumericType> params_;

  const NumericType lateralMu_;
  const NumericType lateralSigma_;
};

// template <class NumericType, int D>
// class DualPearsonModel : public cs::ImplantModel<NumericType, D> {
// public:
//   DualPearsonModel(const NumericType a) : a_(a) {};
//   NumericType getDepthProfile(NumericType depth) override {
//     NumericType mu1 = params.get("mu1");
//     NumericType mu2 = params.get("mu2");
//     NumericType sigma1 = params.get("sigma1");
//     NumericType sigma2 = params.get("sigma2");
//     NumericType gamma1 = params.get("gamma1");
//     NumericType gamma2 = params.get("gamma2");
//     NumericType beta1 = params.get("beta1");
//     NumericType beta2 = params.get("beta2");
//     NumericType r = params.get("r");
//     return (1 - r) * singlePearsonIV(depth, mu1, sigma1, gamma1, beta1) +
//            r * singlePearsonIV(depth, mu2, sigma2, gamma2, beta2);
//   }
//
//   NumericType singlePearsonIV(NumericType depth, NumericType mu,
//                               NumericType sigma, NumericType gamma,
//                               NumericType beta) {
//     depth = depth - mu;
//     NumericType A = 10 * beta - 12 * gamma * gamma - 18;
//     NumericType a = -gamma * sigma * (beta + 3) / A;
//     NumericType b_0 = -sigma * sigma * (4 * beta - 3 * gamma * gamma) / A;
//     NumericType b_1 = a;
//     NumericType b_2 = -(2 * beta - 3 * gamma * gamma - 6) / A;
//     NumericType m = 1 / (2 * b_2);
//     return pow(abs(b_0 + b_1 * depth + b_2 * depth * depth), m) *
//            exp((-(b_1 / b_2 + 2 * a) / sqrt(4 * b_0 * b_2 - b_1 * b_1)) *
//                atan((2 * b_2 * depth + b_1) / sqrt(4 * b_0 * b_2 - b_1 *
//                b_1)));
//   }
//
// private:
//   const NumericType a_ = 1;
// };

} // namespace viennacs
