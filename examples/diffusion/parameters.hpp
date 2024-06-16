#pragma once

#include <vcUtil.hpp>

template <typename T> struct Parameters {
  // Domain
  T gridDelta = 1.5; // nm
  T xExtent = 80.0;  // nm
  T yExtent = 80.0;  // nm

  // Geometry
  T substrateHeight = 50.; // nm
  T coverHeight = 30.;
  T maskHeight = 10.;
  T holeRadius = 13.3;

  // Process
  T duration = 20.;             // s
  T diffusionCoefficient = 1.;  // nm²/s
  T velocity = 0.;              // Advection (not implemented)
  T boundaryValue = 1.;
  T timeStabilityFactor = 0.95; // Stable below 1 (explicit scheme)

  int numThreads = 4;

  Parameters() {}

  void fromMap(std::unordered_map<std::string, std::string> &m) {
    using namespace viennacore;
    util::AssignItems(
        m,
        util::Item{"gridDelta", gridDelta},
        util::Item{"xExtent", xExtent},
        util::Item{"yExtent", yExtent},
        util::Item{"substrateHeight", substrateHeight},
        util::Item{"coverHeight", coverHeight},
        util::Item{"maskHeight", maskHeight},
        util::Item{"holeRadius", holeRadius},
        util::Item{"duration", duration},
        util::Item{"diffusionCoefficient", diffusionCoefficient},
        util::Item{"velocity", velocity},
        util::Item{"boundaryValue", boundaryValue},
        util::Item{"timeStabilityFactor", timeStabilityFactor},
        util::Item{"numThreads", numThreads}
    );
  }
};
