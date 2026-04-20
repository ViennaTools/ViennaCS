#pragma once

#include "csBoundingVolume.hpp"

namespace viennacs {

using namespace viennacore;

/// Helper class to quickly determine the cell in which a given point resides
/// in. To do so, an octree is built around the cell set structure.
template <class T, int D> class BVH {
private:
  using BVPtrType = SmartPointer<BoundingVolume<T, D>>;
  using BoundsType = std::array<VectorType<T, D>, 2>;
  using CellIdsPtr = std::set<unsigned> *;

  unsigned numLayers = 1;
  BVPtrType BV = nullptr;

public:
  BVH(const BoundsType &domainBounds, unsigned layers = 1) : numLayers(layers) {
    BV = BVPtrType::New(domainBounds, numLayers - 1);
  }

  BVPtrType getTopBV() { return BV; }

  void getLowestBVBounds(const Vec3D<T> &point) {
    BV->getBoundingVolumeBounds(point);
  }

  CellIdsPtr getCellIds(const Vec3D<T> &point) { return BV->getCellIds(point); }

  void clearCellIds() { BV->clear(); }

  size_t getTotalCellCount() { return BV->getTotalCellCounts(); }
};

} // namespace viennacs
