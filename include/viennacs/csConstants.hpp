#pragma once

#include <cmath>

namespace viennacs::constants {

template <typename NumericType> struct PearsonIVParameters {
  NumericType mu;
  NumericType sigma;
  NumericType beta;
  NumericType gamma;
};

template <typename NumericType>
NumericType PearsonIV(NumericType x,
                      const PearsonIVParameters<NumericType> &params) {
  if (params.sigma <= NumericType(0))
    return NumericType(0);

  x = x - params.mu;
  NumericType A = 10 * params.beta - 12 * params.gamma * params.gamma - 18;
  if (std::abs(A) <= NumericType(1e-12))
    return NumericType(0);

  NumericType a = -params.gamma * params.sigma * (params.beta + 3) / A;
  NumericType b_0 = -params.sigma * params.sigma *
                    (4 * params.beta - 3 * params.gamma * params.gamma) / A;
  NumericType b_1 = a;
  NumericType b_2 =
      -(2 * params.beta - 3 * params.gamma * params.gamma - 6) / A;
  if (std::abs(b_2) <= NumericType(1e-12))
    return NumericType(0);

  const NumericType discriminant = 4 * b_0 * b_2 - b_1 * b_1;
  if (discriminant <= NumericType(0))
    return NumericType(0);

  NumericType m = 1 / (2 * b_2);
  const NumericType sqrtDisc = std::sqrt(discriminant);
  const NumericType polynomial = std::abs(b_0 + b_1 * x + b_2 * x * x);
  if (polynomial <= NumericType(0))
    return NumericType(0);

  const NumericType result =
      std::pow(polynomial, m) *
      std::exp((-(b_1 / b_2 + 2 * a) / sqrtDisc) *
               std::atan((2 * b_2 * x + b_1) / sqrtDisc));

  if (!std::isfinite(result))
    return NumericType(0);

  return result;
}

} // namespace viennacs::constants
