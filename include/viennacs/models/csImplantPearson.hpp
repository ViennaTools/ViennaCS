#pragma once

#include "../csConstants.hpp"
#include "../csImplantModel.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

namespace viennacs {

namespace impl {

template <typename NumericType>
inline NumericType smoothstep(const NumericType edge0, const NumericType edge1,
                              const NumericType x) {
  if (edge1 <= edge0)
    return x >= edge0 ? NumericType(1) : NumericType(0);

  const auto t = std::clamp((x - edge0) / (edge1 - edge0), NumericType(0),
                            NumericType(1));
  return t * t * (NumericType(3) - NumericType(2) * t);
}

template <typename NumericType>
inline NumericType integrateTrapezoidal(
    const std::function<NumericType(NumericType)> &func, NumericType start,
    NumericType stop, NumericType step) {
  if (stop <= start || step <= NumericType(0))
    return NumericType(0);

  NumericType integral = 0;
  NumericType x0 = start;
  NumericType y0 = func(x0);
  for (NumericType x1 = start + step; x1 <= stop + step / 2; x1 += step) {
    const auto clippedX1 = std::min(x1, stop);
    const auto y1 = func(clippedX1);
    integral += (y0 + y1) * (clippedX1 - x0) / NumericType(2);
    x0 = clippedX1;
    y0 = y1;
    if (clippedX1 >= stop)
      break;
  }
  return integral;
}

} // namespace impl

template <class NumericType, int D>
class ImplantPearsonIV final : public ImplantModel<NumericType, D> {
public:
  ImplantPearsonIV(const constants::PearsonIVParameters<NumericType> &params,
                   NumericType lateralMu, NumericType lateralSigma)
      : params_(params), lateralMu_(lateralMu), lateralSigma_(lateralSigma),
        maxDepth_(std::max(params.mu + NumericType(8) * params.sigma,
                           NumericType(0))),
        maxLateralRange_(std::abs(lateralMu) + NumericType(6) * lateralSigma) {
    const auto integrationStep =
        std::max(params_.sigma / NumericType(50), NumericType(1e-3));
    depthNormalization_ =
        impl::integrateTrapezoidal<NumericType>(
            [this](NumericType depth) { return rawDepthProfile(depth); },
            NumericType(0), maxDepth_, integrationStep);
    if (depthNormalization_ <= NumericType(0))
      depthNormalization_ = NumericType(1);
  }

  NumericType getDepthProfile(NumericType depth) override {
    return rawDepthProfile(depth) / depthNormalization_;
  }

  NumericType getLateralProfile(NumericType offset,
                                NumericType depth) override {
    return (1.0 / (lateralSigma_ * std::sqrt(2 * M_PI))) *
           std::exp(-0.5 * std::pow((offset - lateralMu_) / lateralSigma_, 2));
  }

  NumericType getMaxDepth() override { return maxDepth_; }

  NumericType getMaxLateralRange() override { return maxLateralRange_; }

protected:
  NumericType rawDepthProfile(NumericType depth) const {
    if (depth < NumericType(0))
      return NumericType(0);
    return constants::PearsonIV(depth, params_);
  }

private:
  const constants::PearsonIVParameters<NumericType> params_;

  const NumericType lateralMu_;
  const NumericType lateralSigma_;
  NumericType depthNormalization_ = NumericType(1);
  NumericType maxDepth_ = NumericType(0);
  NumericType maxLateralRange_ = NumericType(0);
};

template <class NumericType, int D>
class ImplantPearsonIVChanneling final : public ImplantModel<NumericType, D> {
public:
  ImplantPearsonIVChanneling(
      const constants::PearsonIVParameters<NumericType> &params,
      NumericType lateralMu, NumericType lateralSigma,
      NumericType tailFraction, NumericType tailStartDepth,
      NumericType tailDecayLength, NumericType tailBlendWidth = NumericType(0))
      : randomImplant_(params, lateralMu, lateralSigma),
        lateralMu_(lateralMu), lateralSigma_(lateralSigma),
        tailFraction_(std::clamp(tailFraction, NumericType(0), NumericType(1))),
        tailStartDepth_(tailStartDepth), tailDecayLength_(tailDecayLength),
        tailBlendWidth_(tailBlendWidth),
        maxDepth_(std::max(randomImplant_.getMaxDepth(),
                           tailStartDepth + NumericType(10) * tailDecayLength)),
        maxLateralRange_(std::abs(lateralMu) + NumericType(6) * lateralSigma) {
    const auto integrationStep =
        std::max(params.sigma / NumericType(50), NumericType(1e-3));
    tailNormalization_ = impl::integrateTrapezoidal<NumericType>(
        [this](NumericType depth) { return rawTailProfile(depth); },
        NumericType(0), maxDepth_, integrationStep);
    if (tailNormalization_ <= NumericType(0))
      tailNormalization_ = NumericType(1);
  }

  NumericType getDepthProfile(NumericType depth) override {
    const auto randomDensity = randomImplant_.getDepthProfile(depth);
    if (tailFraction_ <= NumericType(0) || tailDecayLength_ <= NumericType(0))
      return randomDensity;

    const auto tailDensity = rawTailProfile(depth) / tailNormalization_;
    return (NumericType(1) - tailFraction_) * randomDensity +
           tailFraction_ * tailDensity;
  }

  NumericType getLateralProfile(NumericType offset,
                                NumericType depth) override {
    return (NumericType(1) / (lateralSigma_ * std::sqrt(2 * M_PI))) *
           std::exp(-NumericType(0.5) *
                    std::pow((offset - lateralMu_) / lateralSigma_, 2));
  }

  NumericType getMaxDepth() override { return maxDepth_; }

  NumericType getMaxLateralRange() override { return maxLateralRange_; }

private:
  NumericType rawTailProfile(NumericType depth) const {
    if (depth < NumericType(0) || tailDecayLength_ <= NumericType(0))
      return NumericType(0);
    const auto onset = impl::smoothstep(
        tailStartDepth_ - NumericType(0.5) * tailBlendWidth_,
        tailStartDepth_ + NumericType(0.5) * tailBlendWidth_, depth);
    if (onset <= NumericType(0))
      return NumericType(0);
    return onset *
           std::exp(-(depth - tailStartDepth_) / tailDecayLength_);
  }

  ImplantPearsonIV<NumericType, D> randomImplant_;
  NumericType lateralMu_;
  NumericType lateralSigma_;
  NumericType tailFraction_;
  NumericType tailStartDepth_;
  NumericType tailDecayLength_;
  NumericType tailBlendWidth_;
  NumericType tailNormalization_ = NumericType(1);
  NumericType maxDepth_ = NumericType(0);
  NumericType maxLateralRange_ = NumericType(0);
};

template <class NumericType, int D>
class ImplantDualPearsonIV final : public ImplantModel<NumericType, D> {
public:
  ImplantDualPearsonIV(
      const constants::PearsonIVParameters<NumericType> &headParams,
      const constants::PearsonIVParameters<NumericType> &tailParams,
      NumericType headFraction, NumericType lateralMu,
      NumericType lateralSigma)
      : headImplant_(headParams, lateralMu, lateralSigma),
        tailImplant_(tailParams, lateralMu, lateralSigma),
        headFraction_(
            std::clamp(headFraction, NumericType(0), NumericType(1))),
        maxDepth_(std::max(headImplant_.getMaxDepth(),
                           tailImplant_.getMaxDepth())),
        maxLateralRange_(
            std::max(headImplant_.getMaxLateralRange(),
                     tailImplant_.getMaxLateralRange())) {}

  NumericType getDepthProfile(NumericType depth) override {
    return headFraction_ * headImplant_.getDepthProfile(depth) +
           (NumericType(1) - headFraction_) *
               tailImplant_.getDepthProfile(depth);
  }

  NumericType getLateralProfile(NumericType offset,
                                NumericType depth) override {
    return headFraction_ * headImplant_.getLateralProfile(offset, depth) +
           (NumericType(1) - headFraction_) *
               tailImplant_.getLateralProfile(offset, depth);
  }

  NumericType getMaxDepth() override { return maxDepth_; }

  NumericType getMaxLateralRange() override { return maxLateralRange_; }

private:
  ImplantPearsonIV<NumericType, D> headImplant_;
  ImplantPearsonIV<NumericType, D> tailImplant_;
  NumericType headFraction_;
  NumericType maxDepth_ = NumericType(0);
  NumericType maxLateralRange_ = NumericType(0);
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
