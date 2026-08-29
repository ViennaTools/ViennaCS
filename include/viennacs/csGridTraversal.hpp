#pragma once

#include "csDenseCellSet.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace viennacs {

using namespace viennacore;

/// Maps a lattice coordinate to a cell index of a DenseCellSet.
///
/// A cell set stores its cells compactly: `fromLevelSets` appends a cell only
/// where material was found, in the order the level-set iterator walked, so
/// the cell array carries no lattice structure and a point has to be looked up
/// through the BVH. A ray walking the grid asks the opposite question -- "which
/// cell is at (i,j,k)" -- once per cell it crosses, which is far too often for
/// a tree descent.
///
/// This builds the dense table once: one int per lattice site, -1 where the
/// cell set has no cell. The lattice never changes during a simulation (only
/// the filling fractions do), so it is built once and read for the rest of the
/// run.
template <class T, int D> class LatticeMap {
  std::array<int, D> dims_{};
  std::array<T, D> min_{};
  T gridDelta_ = 0;
  std::vector<int> cellIds_;

public:
  LatticeMap() = default;

  explicit LatticeMap(SmartPointer<DenseCellSet<T, D>> cellSet) {
    build(*cellSet);
  }

  explicit LatticeMap(DenseCellSet<T, D> &cellSet) { build(cellSet); }

  void build(SmartPointer<DenseCellSet<T, D>> cellSet) { build(*cellSet); }

  void build(DenseCellSet<T, D> &cellSet) {
    gridDelta_ = cellSet.getGridDelta();
    const auto bb = cellSet.getBoundingBox();

    size_t total = 1;
    for (int d = 0; d < D; ++d) {
      min_[d] = bb[0][d];
      // The bounding box spans the cell corners, so the number of cells along
      // an axis is the span divided by the spacing. Round rather than
      // truncate: the division is exact in exact arithmetic and a hair short
      // of it in floating point.
      dims_[d] = static_cast<int>(std::round((bb[1][d] - bb[0][d]) / gridDelta_));
      if (dims_[d] < 1)
        dims_[d] = 1;
      total *= static_cast<size_t>(dims_[d]);
    }

    cellIds_.assign(total, -1);

    // Node 0 of a cell is its minimum corner -- the same convention
    // DenseCellSet uses when it recovers a cell centre.
    const auto &elements = cellSet.getElements();
    const auto &nodes = cellSet.getNodes();
    for (size_t c = 0; c < elements.size(); ++c) {
      const auto &corner = nodes[elements[c][0]];
      std::array<int, D> idx{};
      bool inside = true;
      for (int d = 0; d < D; ++d) {
        idx[d] = static_cast<int>(std::round((corner[d] - min_[d]) / gridDelta_));
        if (idx[d] < 0 || idx[d] >= dims_[d])
          inside = false;
      }
      if (inside)
        cellIds_[flatten(idx)] = static_cast<int>(c);
    }
  }

  const std::array<int, D> &dims() const { return dims_; }
  const std::array<T, D> &minCorner() const { return min_; }
  T gridDelta() const { return gridDelta_; }

  /// The cell at a lattice coordinate, or -1 where the cell set has none.
  int cellId(const std::array<int, D> &idx) const {
    for (int d = 0; d < D; ++d)
      if (idx[d] < 0 || idx[d] >= dims_[d])
        return -1;
    return cellIds_[flatten(idx)];
  }

private:
  /// Private on purpose: which axis varies fastest is an internal detail, and
  /// the only requirement on it is that `build` and `cellId` agree. A caller
  /// that indexed its own array this way would be depending on something this
  /// class does not promise.
  size_t flatten(const std::array<int, D> &idx) const {
    size_t flat = 0;
    for (int d = D - 1; d >= 0; --d)
      flat = flat * static_cast<size_t>(dims_[d]) + static_cast<size_t>(idx[d]);
    return flat;
  }
};

namespace detail {
/// Wraps the lattice coordinate so that D is deduced from the LatticeMap
/// alone: std::array takes a size_t extent while the lattice is templated on
/// int D, and deducing from both at once is a mismatch.
template <int D> struct LatticeIndexOf {
  using type = std::array<int, static_cast<std::size_t>(D)>;
};
} // namespace detail

template <int D> using LatticeIndex = typename detail::LatticeIndexOf<D>::type;

/// THE FILL FIELD, read the two ways this method needs it.
///
/// These are free functions on purpose. The gradient below is the single
/// definition of the interface geometry: VoxelInteraction normalises it into
/// the surface normal a rate law sees, and VoxelAdvance takes its magnitude
/// as the interface area that turns a velocity into a volume. The claim that
/// "a mechanism sees an area and a normal that agree with each other" is only
/// structurally true while both come from here; as two hand-copied stencils
/// it was one edit away from a silent physics divergence.

/// The fill at a lattice coordinate, zero off the lattice: a ray that leaves
/// the grid must find nothing to hit.
template <class T, int D>
inline T fillFieldAt(const LatticeMap<T, D> &lattice,
                     const std::vector<T> &fill,
                     const LatticeIndex<D> &idx) {
  const int id = lattice.cellId(idx);
  return id < 0 ? T(0) : fill[id];
}

/// The same, for the purpose of a DERIVATIVE, with the opposite boundary
/// convention: the lattice edge is a cut through the material, not a surface,
/// so the field continues across it with zero gradient. Reading zero instead
/// would give every cell on the edge of the domain a normal pointing out of
/// it, and the bottom row of a solid block would grow downwards out of the
/// grid.
template <class T, int D>
inline T fillFieldClamped(const LatticeMap<T, D> &lattice,
                          const std::vector<T> &fill,
                          const LatticeIndex<D> &idx) {
  auto inside = idx;
  for (int d = 0; d < D; ++d)
    inside[d] = std::min(std::max(inside[d], 0), lattice.dims()[d] - 1);
  return fillFieldAt(lattice, fill, inside);
}

/// -grad(f) in physical units, unnormalised: its direction is the outward
/// normal, and its magnitude is the interfacial area per unit volume of a
/// volume-fraction field.
///
/// `wide` sweeps the whole 3^D neighbourhood with a (1,2,1) weight across
/// each axis other than the one being differenced -- Youngs' stencil, as in
/// volume of fluid. The narrow form differences only the two face
/// neighbours: cheaper, and markedly more anisotropic.
template <class T, int D>
inline Vec3D<T> fillFieldGradient(const LatticeMap<T, D> &lattice,
                                  const std::vector<T> &fill,
                                  const LatticeIndex<D> &idx,
                                  bool wide = true) {
  const T delta = lattice.gridDelta();
  Vec3D<T> g{0, 0, 0};

  if (!wide) {
    for (int d = 0; d < D; ++d) {
      auto lo = idx, hi = idx;
      lo[d] -= 1;
      hi[d] += 1;
      g[d] = (fillFieldClamped(lattice, fill, lo) -
              fillFieldClamped(lattice, fill, hi)) /
             (T(2) * delta);
    }
    return g;
  }

  int span = 1;
  for (int d = 0; d < D; ++d)
    span *= 3;
  // the normalisation is the total weight on one side, 2^(D-1), times the
  // two-cell span of the central difference
  T total = 1;
  for (int d = 0; d < D - 1; ++d)
    total *= 4;
  for (int s = 0; s < span; ++s) {
    std::array<int, D> offset{};
    int rem = s;
    for (int d = 0; d < D; ++d) {
      offset[d] = rem % 3 - 1;
      rem /= 3;
    }
    auto probe = idx;
    for (int d = 0; d < D; ++d)
      probe[d] += offset[d];
    const T f = fillFieldClamped(lattice, fill, probe);
    for (int d = 0; d < D; ++d) {
      if (offset[d] == 0)
        continue;
      T weight = 1;
      for (int k = 0; k < D; ++k)
        if (k != d)
          weight *= (offset[k] == 0) ? T(2) : T(1);
      g[d] -= static_cast<T>(offset[d]) * weight * f;
    }
  }
  for (int d = 0; d < D; ++d)
    g[d] /= (total * T(2) * delta);
  return g;
}

/// One cell of a ray's passage through the grid.
template <class T, int D> struct GridStep {
  std::array<int, D> index{}; ///< lattice coordinate
  T tEntry = 0;               ///< distance along the ray where it enters
  T tExit = 0;                ///< and where it leaves
};

/// Walks the cells a ray crosses, in order, by the method of Amanatides and
/// Woo, "A Fast Voxel Traversal Algorithm for Ray Tracing" (Eurographics '87).
///
/// The traversal is deliberately ignorant of what the cells contain. It
/// enumerates them with the distances at which the ray enters and leaves each,
/// and a visitor decides what that means -- whether a partially filled cell
/// interacts, how a path length accumulates, where a surface is taken to be.
/// That decision is the physics of a voxel method and does not belong here.
///
/// Distances are in world units, so `direction` must be normalised.
template <class T, int D> class GridTraversal {
  std::array<int, D> dims_{};
  std::array<T, D> min_{};
  T gridDelta_ = 0;

public:
  GridTraversal() = default;

  GridTraversal(const std::array<T, D> &minCorner,
                const std::array<int, D> &dims, T gridDelta)
      : dims_(dims), min_(minCorner), gridDelta_(gridDelta) {}

  explicit GridTraversal(const LatticeMap<T, D> &lattice)
      : dims_(lattice.dims()), min_(lattice.minCorner()),
        gridDelta_(lattice.gridDelta()) {}

  /// Clips the ray to the grid. Returns false if it misses entirely; otherwise
  /// `tNear`/`tFar` bracket the part inside, with tNear >= 0.
  bool clip(const std::array<T, D> &origin, const std::array<T, D> &direction,
            T &tNear, T &tFar) const {
    tNear = 0;
    tFar = std::numeric_limits<T>::max();

    for (int d = 0; d < D; ++d) {
      const T lo = min_[d];
      const T hi = min_[d] + gridDelta_ * static_cast<T>(dims_[d]);
      if (std::abs(direction[d]) < eps_) {
        // Parallel to this pair of faces: inside for all t, or never.
        if (origin[d] < lo || origin[d] > hi)
          return false;
        continue;
      }
      T t0 = (lo - origin[d]) / direction[d];
      T t1 = (hi - origin[d]) / direction[d];
      if (t0 > t1)
        std::swap(t0, t1);
      tNear = std::max(tNear, t0);
      tFar = std::min(tFar, t1);
      if (tNear > tFar)
        return false;
    }
    return tFar > tNear;
  }

  /// Visits every cell the ray crosses. `visitor(GridStep)` returns false to
  /// stop the walk -- which is what an interaction does. Returns the number of
  /// cells visited.
  template <class Visitor>
  size_t traverse(const std::array<T, D> &origin,
                  const std::array<T, D> &direction, Visitor &&visitor) const {
    T tNear, tFar;
    if (!clip(origin, direction, tNear, tFar))
      return 0;

    std::array<int, D> index{};
    std::array<int, D> step{};
    std::array<T, D> tMax{};
    std::array<T, D> tDelta{};
    const T inf = std::numeric_limits<T>::max();

    for (int d = 0; d < D; ++d) {
      // Enter the grid, then locate the cell. Nudging in by a fraction of a
      // cell keeps a ray that enters exactly on a face out of the neighbour.
      const T entry = origin[d] + direction[d] * (tNear + eps_ * gridDelta_);
      int i = static_cast<int>(std::floor((entry - min_[d]) / gridDelta_));
      index[d] = std::min(std::max(i, 0), dims_[d] - 1);

      if (std::abs(direction[d]) < eps_) {
        step[d] = 0;
        tMax[d] = inf;
        tDelta[d] = inf;
        continue;
      }
      step[d] = direction[d] > 0 ? 1 : -1;
      tDelta[d] = gridDelta_ / std::abs(direction[d]);
      const T boundary =
          min_[d] + gridDelta_ * static_cast<T>(index[d] + (step[d] > 0 ? 1 : 0));
      tMax[d] = (boundary - origin[d]) / direction[d];
    }

    size_t visited = 0;
    T tEntry = tNear;
    while (true) {
      // The next face crossed is the nearest of the per-axis candidates.
      int axis = 0;
      T tNext = tMax[0];
      for (int d = 1; d < D; ++d) {
        if (tMax[d] < tNext) {
          tNext = tMax[d];
          axis = d;
        }
      }
      const T tOut = std::min(tNext, tFar);

      ++visited;
      if (!visitor(GridStep<T, D>{index, tEntry, tOut}))
        break;

      if (tNext >= tFar)
        break;

      index[axis] += step[axis];
      if (index[axis] < 0 || index[axis] >= dims_[axis])
        break;
      tEntry = tNext;
      tMax[axis] += tDelta[axis];
    }
    return visited;
  }

private:
  static constexpr T eps_ = static_cast<T>(1e-6);
};

} // namespace viennacs
