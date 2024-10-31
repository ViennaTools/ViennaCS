#pragma once

#include "vcUtil.hpp"

namespace viennacs {

template <typename NumericType, int D> class ImplantModel {
public:
  virtual ~ImplantModel() = default;
  virtual NumericType getDepthProfile(NumericType depth, util::Parameters params) {
    // profile of empirical implant model in z direction
    return 0.;
  }

  virtual NumericType getLateralProfile(NumericType offset, NumericType depth, util::Parameters params) {
    // profile of empirical implant model in lateral directions
    NumericType stddev = params.get("lateral_sigma");
    NumericType mean = params.get("lateral_mu");
    return (1.0 / (stddev * std::sqrt(2 * M_PI))) *
       std::exp(-0.5 * std::pow((offset - mean) / stddev, 2));

  }

  virtual NumericType getMaxDepth() {
    // maximum depth of the implant
    return 0.;
  }
};

} // namespace viennacs