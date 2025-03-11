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
  x = x - params.mu;
  NumericType A = 10 * params.beta - 12 * params.gamma * params.gamma - 18;
  NumericType a = -params.gamma * params.sigma * (params.beta + 3) / A;
  NumericType b_0 = -params.sigma * params.sigma *
                    (4 * params.beta - 3 * params.gamma * params.gamma) / A;
  NumericType b_1 = a;
  NumericType b_2 =
      -(2 * params.beta - 3 * params.gamma * params.gamma - 6) / A;
  NumericType m = 1 / (2 * b_2);
  return pow(abs(b_0 + b_1 * x + b_2 * x * x), m) *
         exp((-(b_1 / b_2 + 2 * a) / sqrt(4 * b_0 * b_2 - b_1 * b_1)) *
             atan((2 * b_2 * x + b_1) / sqrt(4 * b_0 * b_2 - b_1 * b_1)));
}

} // namespace viennacs::constants
