#pragma once

#include <raySource.hpp>

#include <vcVectorUtil.hpp>

namespace viennacs {

using namespace viennacore;

template <typename NumericType, int D>
class PointSource : public viennaray::Source<NumericType> {
  const unsigned mNumPoints;
  const Vec3D<NumericType> origin;
  const Vec3D<NumericType> direction;

public:
  PointSource(Vec3D<NumericType> passedOrigin,
              Vec3D<NumericType> passedDirection,
              std::array<int, 5> &pTraceSettings, const size_t pNumPoints)
      : origin(passedOrigin), direction(passedDirection),
        mNumPoints(pNumPoints) {}

  Vec2D<Vec3D<NumericType>>
  getOriginAndDirection(const size_t idx, RNG &rngState) const override {
    return {origin, direction};
  }

  size_t getNumPoints() const override { return mNumPoints; }

  NumericType getSourceArea() const override { return 1.; }
};

} // namespace viennacs
