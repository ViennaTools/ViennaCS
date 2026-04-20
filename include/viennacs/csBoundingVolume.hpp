#pragma once

#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

#include <array>
#include <set>

namespace viennacs {

using namespace viennacore;

template <class T, int D> class BoundingVolume {
  using BVPtrType = SmartPointer<BoundingVolume>;
  using BoundsType = std::array<VectorType<T, D>, 2>;
  using CellIdsPtr = std::set<unsigned> *;

  static constexpr int numCells = 1 << D;
  std::array<std::set<unsigned>, numCells> cellIds;
  std::array<BoundsType, numCells> bounds;
  std::array<BVPtrType, numCells> links;
  int layer = -1;

public:
  BoundingVolume() = default;
  BoundingVolume(const BoundsType &outerBound, int thisLayer)
      : layer(thisLayer) {

    if constexpr (D == 3)
      buildBounds3D(outerBound);
    else
      buildBounds2D(outerBound);

    if (layer > 0) {
      for (size_t i = 0; i < numCells; i++) {
        links[i] = BVPtrType::New(bounds[i], layer - 1);
      }
    }
  }

  CellIdsPtr getCellIds(const Vec3D<T> &point) {
    auto vid = getVolumeIndex(point);
    // assert(vid < numCells && "Point in invalid BV");

    if (vid == numCells)
      return nullptr;

    if (layer == 0) {
      return &cellIds[vid];
    }

    return links[vid]->getCellIds(point);
  }

  void getBoundingVolumeBounds(const std::array<T, D> &point) {
    auto vid = getVolumeIndex(point);
    assert(vid < numCells && "Point in invalid BV");

    if (layer == 0) {
      printBound(vid);
      return;
    }

    links[vid]->getBoundingVolumeBounds(point);
  }

  size_t getTotalCellCounts() {
    size_t count = 0;
    if (layer == 0) {
      for (size_t i = 0; i < numCells; i++) {
        count += cellIds[i].size();
      }
      return count;
    } else {
      for (size_t i = 0; i < numCells; i++) {
        count += links[i]->getTotalCellCounts();
      }
    }
    return count;
  }

  void clear() {
    if (layer == 0) {
      for (size_t i = 0; i < numCells; i++) {
        cellIds[i].clear();
      }
    } else {
      for (size_t i = 0; i < numCells; i++) {
        links[i]->clear();
      }
    }
  }

  BVPtrType getLink(const Vec3D<T> &point) {
    auto vid = getVolumeIndex(point);
    return getLink(vid);
  }

  BVPtrType getLink(size_t vid) { return links[vid]; }

  size_t getVolumeIndex(const Vec3D<T> &point) {
    size_t vid = numCells;
    for (size_t idx = 0; idx < numCells; idx++) {
      if (insideVolume(point, idx)) {
        vid = idx;
        break;
      }
    }
    return vid;
  }

  bool insideVolume(const Vec3D<T> &p, const size_t idx) const {
    if constexpr (D == 3) {
      return p[0] > bounds[idx][0][0] && p[0] <= bounds[idx][1][0] &&
             p[1] > bounds[idx][0][1] && p[1] <= bounds[idx][1][1] &&
             p[2] > bounds[idx][0][2] && p[2] <= bounds[idx][1][2];
    } else {
      return p[0] > bounds[idx][0][0] && p[0] <= bounds[idx][1][0] &&
             p[1] > bounds[idx][0][1] && p[1] <= bounds[idx][1][1];
    }
  }

  void printBound(size_t vid) {
    std::cout << "Bounding volume span: [";
    for (int i = 0; i < D; i++)
      std::cout << bounds[vid][0][i] << ", ";
    std::cout << "] - [";
    for (int i = 0; i < D; i++)
      std::cout << bounds[vid][1][i] << ", ";
    std::cout << "]\n";
  }

private:
  void buildBounds2D(const BoundsType &outerBound) {
    T xExt = (outerBound[1][0] - outerBound[0][0]) / 2.;
    T yExt = (outerBound[1][1] - outerBound[0][1]) / 2.;

    const auto BVH1 = std::array<Vec2D<T>, 2>{
        Vec2D<T>{outerBound[0][0], outerBound[0][1]},
        Vec2D<T>{outerBound[0][0] + xExt, outerBound[0][1] + yExt}};
    const auto BVH2 = std::array<Vec2D<T>, 2>{
        Vec2D<T>{outerBound[0][0] + xExt, outerBound[0][1]},
        Vec2D<T>{outerBound[0][0] + 2 * xExt, outerBound[0][1] + yExt}};
    const auto BVH3 = std::array<Vec2D<T>, 2>{
        Vec2D<T>{outerBound[0][0] + xExt, outerBound[0][1] + yExt},
        Vec2D<T>{outerBound[0][0] + 2 * xExt, outerBound[0][1] + 2 * yExt}};
    const auto BVH4 = std::array<Vec2D<T>, 2>{
        Vec2D<T>{outerBound[0][0], outerBound[0][1] + yExt},
        Vec2D<T>{outerBound[0][0] + xExt, outerBound[0][1] + 2 * yExt}};

    bounds =
        std::array<std::array<Vec2D<T>, 2>, numCells>{BVH1, BVH2, BVH3, BVH4};
  }

  void buildBounds3D(const BoundsType &outerBound) {
    auto xExt = (outerBound[1][0] - outerBound[0][0]) / T(2);
    auto yExt = (outerBound[1][1] - outerBound[0][1]) / T(2);
    auto zExt = (outerBound[1][2] - outerBound[0][2]) / T(2);

    const auto BVH1 = std::array<Vec3D<T>, 2>{
        Vec3D<T>{outerBound[0][0], outerBound[0][1], outerBound[0][2]},
        Vec3D<T>{outerBound[0][0] + xExt, outerBound[0][1] + yExt,
                 outerBound[0][2] + zExt}};
    const auto BVH2 = std::array<Vec3D<T>, 2>{
        Vec3D<T>{outerBound[0][0] + xExt, outerBound[0][1], outerBound[0][2]},
        Vec3D<T>{outerBound[0][0] + 2 * xExt, outerBound[0][1] + yExt,
                 outerBound[0][2] + zExt}};
    const auto BVH3 = std::array<Vec3D<T>, 2>{
        Vec3D<T>{outerBound[0][0] + xExt, outerBound[0][1] + yExt,
                 outerBound[0][2]},
        Vec3D<T>{outerBound[0][0] + 2 * xExt, outerBound[0][1] + 2 * yExt,
                 outerBound[0][2] + zExt}};
    const auto BVH4 = std::array<Vec3D<T>, 2>{
        Vec3D<T>{outerBound[0][0], outerBound[0][1] + yExt, outerBound[0][2]},
        Vec3D<T>{outerBound[0][0] + xExt, outerBound[0][1] + 2 * yExt,
                 outerBound[0][2] + zExt}};

    // top
    const auto BVH5 = std::array<Vec3D<T>, 2>{
        Vec3D<T>{outerBound[0][0], outerBound[0][1], outerBound[0][2] + zExt},
        Vec3D<T>{outerBound[0][0] + xExt, outerBound[0][1] + yExt,
                 outerBound[0][2] + 2 * zExt}};
    const auto BVH6 = std::array<Vec3D<T>, 2>{
        Vec3D<T>{outerBound[0][0] + xExt, outerBound[0][1],
                 outerBound[0][2] + zExt},
        Vec3D<T>{outerBound[0][0] + 2 * xExt, outerBound[0][1] + yExt,
                 outerBound[0][2] + 2 * zExt}};
    const auto BVH7 = std::array<Vec3D<T>, 2>{
        Vec3D<T>{outerBound[0][0] + xExt, outerBound[0][1] + yExt,
                 outerBound[0][2] + zExt},
        Vec3D<T>{outerBound[0][0] + 2 * xExt, outerBound[0][1] + 2 * yExt,
                 outerBound[0][2] + 2 * zExt}};
    const auto BVH8 = std::array<Vec3D<T>, 2>{
        Vec3D<T>{outerBound[0][0], outerBound[0][1] + yExt,
                 outerBound[0][2] + zExt},
        Vec3D<T>{outerBound[0][0] + xExt, outerBound[0][1] + 2 * yExt,
                 outerBound[0][2] + 2 * zExt}};
    bounds = std::array<std::array<Vec3D<T>, 2>, numCells>{
        BVH1, BVH2, BVH3, BVH4, BVH5, BVH6, BVH7, BVH8};
  }
};

} // namespace viennacs