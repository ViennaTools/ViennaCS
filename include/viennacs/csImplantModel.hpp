#pragma once

namespace viennacs {

template <typename NumericType, int D> class ImplantModel {
public:
  virtual ~ImplantModel() = default;

  virtual NumericType getDepthProfile(NumericType depth) {
    // profile of empirical implant model in z direction
    return 0.;
  }

  virtual NumericType getLateralProfile(NumericType offset, NumericType depth) {
    // profile of empirical implant model in lateral directions
    NumericType stddev = 0.1;
    NumericType mean = 0;
    const double pi = 3.14159265358979323846;
    return (1.0 / (stddev * std::sqrt(2 * pi))) *
        std::exp(-0.5 * std::pow((offset - mean) / stddev, 2));
  }

  virtual NumericType getMaxDepth() {
    // maximum depth of the implant
    return 0.;
  }
};

} // namespace viennacs