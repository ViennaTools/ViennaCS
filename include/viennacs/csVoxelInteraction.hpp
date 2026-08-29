#pragma once

#include "csVoxelEmbreeTraversal.hpp"
#include "csVoxelSurface.hpp"

#include <vcRNG.hpp>

#include <cmath>
#include <memory>

namespace viennacs {

using namespace viennacore;

/// How a ray decides that it has met the surface, and what normal it sees
/// there.
///
/// This is the physics of a voxel method, and it is where a voxel method
/// differs from a level-set one. A level set carries the surface as a
/// sub-grid quantity and hands back a smooth normal. A voxel grid carries
/// filling fractions, and neither the surface position nor the normal is
/// written down anywhere -- both have to be decided.
///
/// THE INTERACTION RULE. A cell that is 40% solid presents 40% of its
/// cross-section, so a ray crossing it interacts with probability 0.4 and
/// otherwise passes through. Over many rays that places the effective surface
/// between the cell faces, which is how a voxel method recovers sub-grid
/// surface position without reconstructing a surface. Chord length enters
/// through 1 - (1-f)^(L/delta), so a ray clipping a corner is less likely to
/// interact than one crossing the whole cell; for an axis-aligned crossing,
/// L = delta and the probability is exactly f.
///
/// THE NORMAL. Two answers, and choosing between them is the experiment:
///
///   Face          the outward normal of the face the ray entered through.
///                 Quantised to 2D directions, so an ion's angle of incidence
///                 is quantised with it. This is the staircase in its
///                 undisguised form.
///
///   FillGradient  -grad(f) by three-point central differences. Cheap, and
///                 exact where the surface lies along a lattice axis or its
///                 diagonal -- and worst between them, by around ten degrees.
///                 The stencil is anisotropic, so this is not what a voxel
///                 method should be judged on.
///
///   FillGradientYoungs
///                 -grad(f) over the whole 3^D neighbourhood, weighted
///                 (1,2,1) across each axis, as in Youngs' method for volume
///                 of fluid. Isotropic enough that the direction of the
///                 surface stops mattering.
///
/// The last two are still surface reconstruction, implicit rather than
/// explicit: a smoothed gradient of a volume fraction is an interface normal
/// by another name. The claim that voxels avoid reconstructing a surface does
/// not survive contact with an ion yield that depends on cos(theta). What is
/// true is that the reconstruction is local and never stored.
///
/// A comparison against a level set that reports only one of these has
/// answered only part of the question, so all three are here and switchable.
enum class NormalEstimator { Face, FillGradient, FillGradientYoungs };

/// How a ray finds the cell it interacts with. The physics -- the acceptance
/// probability, the normal, the arming distance -- is identical either way;
/// what differs is the cost of the question "which cells does this ray cross".
///
///   GridDDA     Amanatides-Woo cell walking. O(cells crossed) per ray, which
///               is O(depth of gas) -- the reference implementation.
///   EmbreeBVH   the cells as embree user primitives, one BVH descent per
///               ray segment -- the same acceleration library the level-set
///               arm traces against.
///   Hybrid      each engine where it is strong. A primary ray flies the
///               whole gas column, which is the BVH's regime; a re-emitted
///               segment usually lands a few cells away, and marching three
///               cells beats any tree query. The two cases are told apart by
///               the arming distance: only a re-emitted segment carries one.
enum class TraversalEngine { GridDDA, EmbreeBVH, Hybrid };

template <class T, int D> struct VoxelHit {
  int cellId = -1;
  std::array<int, D> index{};
  T distance = 0;
  Vec3D<T> point{0, 0, 0};
  Vec3D<T> normal{0, 0, 0};
  int enteredAxis = -1; ///< the axis whose face the ray crossed to get in
  int enteredSign = 0;  ///< -1 if it entered through the low face, +1 the high

  bool hit() const { return cellId >= 0; }
};

/// Finds where a ray meets a voxel geometry described by filling fractions.
template <class T, int D> class VoxelInteraction {
  const LatticeMap<T, D> *lattice_ = nullptr;
  const std::vector<T> *fill_ = nullptr;
  GridTraversal<T, D> traversal_;
  NormalEstimator estimator_ = NormalEstimator::Face;
  TraversalEngine engine_ = TraversalEngine::GridDDA;
  // A cache of the geometry, not part of the object's value: shared so a
  // copied interaction keeps tracing against the same built scene, mutable so
  // `prepare` can rebuild it on a const object.
  mutable std::shared_ptr<EmbreeCellTraversal<T, D>> bvh_;
  // The gradient normal per cell, computed once per trace instead of once
  // per hit: the fills are frozen while the rays fly, so the 3^D stencil a
  // hit pays -- tens of millions of times per step -- always returns the
  // same vector. `valid` is -1 where nothing is cached (fall through to the
  // direct stencil), 0 where the gradient is degenerate (the hit's own face
  // normal stands in, as gradientNormal itself falls back).
  mutable std::vector<Vec3D<T>> normalCache_;
  mutable std::vector<signed char> normalValid_;

public:
  VoxelInteraction(const LatticeMap<T, D> &lattice, const std::vector<T> &fill,
                   NormalEstimator estimator = NormalEstimator::Face)
      : lattice_(&lattice), fill_(&fill), traversal_(lattice),
        estimator_(estimator) {}

  void setNormalEstimator(NormalEstimator e) {
    estimator_ = e;
    normalValid_.clear(); // the cache belongs to one estimator
    normalCache_.clear();
  }
  NormalEstimator normalEstimator() const { return estimator_; }

  void setTraversalEngine(TraversalEngine e) { engine_ = e; }
  TraversalEngine traversalEngine() const { return engine_; }

  /// Builds the BVH from the CURRENT fills. Must be called before a parallel
  /// tracing region whenever the engine is EmbreeBVH -- the fills change every
  /// step, and rays must never race a build.
  void prepare() const {
    if (engine_ != TraversalEngine::GridDDA) {
      if (!bvh_)
        bvh_ = std::make_shared<EmbreeCellTraversal<T, D>>();
      bvh_->build(*lattice_, *fill_);
    }
    if (estimator_ == NormalEstimator::Face)
      return; // the face normal is per-hit by nature; nothing to cache
    const auto &dims = lattice_->dims();
    size_t sites = 1;
    for (int d = 0; d < D; ++d)
      sites *= static_cast<size_t>(dims[d]);
    normalCache_.assign(fill_->size(), Vec3D<T>{0, 0, 0});
    normalValid_.assign(fill_->size(), 0);
    const bool wide = estimator_ == NormalEstimator::FillGradientYoungs;
#pragma omp parallel for schedule(static)
    for (long long flat = 0; flat < static_cast<long long>(sites); ++flat) {
      std::array<int, D> idx{};
      size_t rem = static_cast<size_t>(flat);
      for (int d = 0; d < D; ++d) {
        idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
        rem /= static_cast<size_t>(dims[d]);
      }
      const int id = lattice_->cellId(idx);
      if (id < 0 || (*fill_)[id] <= T(0))
        continue; // only a cell holding material can be hit
      // A zero sentinel face normal: gradientNormal returns it exactly when
      // the gradient is degenerate, which is the case where the per-hit face
      // normal must stand in.
      const auto n = gradientNormal(idx, Vec3D<T>{0, 0, 0}, wide);
      if (n[0] != T(0) || n[1] != T(0) || n[2] != T(0)) {
        normalCache_[id] = n;
        normalValid_[id] = 1;
      }
    }
  }

  /// The estimator's normal at a hit, from the per-trace cache when one is
  /// built, by the direct stencil when not -- same values either way.
  Vec3D<T> hitNormal(int id, const std::array<int, D> &idx,
                     const Vec3D<T> &faceNormal) const {
    if (estimator_ == NormalEstimator::Face)
      return faceNormal;
    if (normalValid_.size() == fill_->size())
      return normalValid_[id] ? normalCache_[id] : faceNormal;
    return gradientNormal(idx, faceNormal,
                          estimator_ == NormalEstimator::FillGradientYoungs);
  }

  /// The filling fraction at a lattice coordinate; zero where there is no
  /// cell, because a ray that leaves the grid must find nothing to hit.
  T fillAt(const std::array<int, D> &idx) const {
    return fillFieldAt(*lattice_, *fill_, idx);
  }

  /// The same, for a DERIVATIVE: the lattice boundary is a cut through the
  /// material, not a surface, so the field continues across it with zero
  /// gradient. Reading zero instead would give every cell on the edge of the
  /// domain a normal pointing out of it.
  T fillClamped(const std::array<int, D> &idx) const {
    return fillFieldClamped(*lattice_, *fill_, idx);
  }

  /// -grad(f), normalised. Falls back to the face normal where the gradient
  /// vanishes, which happens inside a uniformly filled region and at an
  /// isolated cell.
  ///
  /// `wide` sweeps the whole 3^D neighbourhood with a (1,2,1) weight across
  /// each axis other than the one being differenced -- Youngs' stencil. The
  /// narrow form differences only the two face neighbours, which is cheaper
  /// and markedly more anisotropic.
  Vec3D<T> gradientNormal(const std::array<int, D> &idx,
                          const Vec3D<T> &faceNormal, bool wide) const {
    // the stencil itself is shared with the advance, so the normal a rate law
    // sees and the area that carries its velocity cannot drift apart
    auto n = fillFieldGradient(*lattice_, *fill_, idx, wide);

    T norm = 0;
    for (int d = 0; d < D; ++d)
      norm += n[d] * n[d];
    norm = std::sqrt(norm);
    if (norm < T(1e-12))
      return faceNormal;
    for (int d = 0; d < D; ++d)
      n[d] /= norm;
    return n;
  }

  /// Walks the ray until a cell interacts. `rng` is consumed once per
  /// partially filled cell crossed.
  /// `armAfter` suppresses interaction until the ray has travelled that far
  /// from its origin. It is the voxel analogue of embree's tnear, which is what
  /// keeps the level-set arm's ray from re-hitting the primitive it just left.
  ///
  /// A DISTANCE, deliberately, and not a geometric condition. A fractional
  /// interface is two or three cells thick, so a ray re-emitted from inside it
  /// would otherwise interact again immediately and, where the sticking is
  /// small, deposit its full weight over and over. Two attempts to prevent that
  /// were worse:
  ///
  ///   - displacing the origin along the normal pushed the ray one to three
  ///     cells off the surface at every bounce, which in a trench is several
  ///     percent of the width each time and lifts rays out of the feature;
  ///   - waiting for a cell holding no material stranded grazing rays, which
  ///     hug the interface, never meet an empty cell, and are lost.
  ///
  /// A distance gate does neither: the ray stays where it was emitted, and it
  /// is free to strike anything beyond one interface thickness -- including the
  /// same wall further along, which is a different patch of surface and should
  /// be hit.
  template <class RNG>
  VoxelHit<T, D> firstHit(const std::array<T, D> &origin,
                          const std::array<T, D> &direction, RNG &rng,
                          T armAfter = T(0)) const {
    const bool useBVH =
        bvh_ && bvh_->built() &&
        (engine_ == TraversalEngine::EmbreeBVH ||
         (engine_ == TraversalEngine::Hybrid && armAfter <= T(0)));
    if (useBVH) {
      // One seed per ray segment: it decides every acceptance on this
      // segment idempotently, and a re-emitted segment draws a new one.
      const std::uint64_t seed = (static_cast<std::uint64_t>(rng()) << 32) ^
                                 static_cast<std::uint64_t>(rng());
      const auto raw = bvh_->firstHit(origin, direction, seed, armAfter);
      VoxelHit<T, D> result;
      if (!raw.hit())
        return result;
      result.cellId = raw.cellId;
      result.index = raw.index;
      result.distance = raw.tEntry;
      for (int d = 0; d < D; ++d)
        result.point[d] = origin[d] + direction[d] * raw.tEntry;
      Vec3D<T> faceNormal{0, 0, 0};
      if (raw.axis >= 0) {
        result.enteredAxis = raw.axis;
        result.enteredSign = raw.sign;
      } else {
        // The origin lay inside the cell: same fallback as the DDA's first
        // cell -- the face of the axis the ray travels most steeply against.
        int axis = 0;
        T steepest = 0;
        for (int d = 0; d < D; ++d)
          if (std::abs(direction[d]) > steepest) {
            steepest = std::abs(direction[d]);
            axis = d;
          }
        result.enteredAxis = axis;
        result.enteredSign = direction[axis] > 0 ? -1 : 1;
      }
      faceNormal[result.enteredAxis] = static_cast<T>(result.enteredSign);
      result.normal = hitNormal(raw.cellId, raw.index, faceNormal);
      return result;
    }

    const T delta = lattice_->gridDelta();
    std::uniform_real_distribution<T> uniform(T(0), T(1));

    VoxelHit<T, D> result;
    std::array<int, D> previous{};
    bool havePrevious = false;

    traversal_.traverse(origin, direction, [&](GridStep<T, D> step) {
      if (step.tExit <= armAfter) {
        previous = step.index;
        havePrevious = true;
        return true; // still within the interface it was emitted from
      }
      const T f = fillAt(step.index);
      if (f > T(0)) {
        // A partially filled cell transmits; a full one always interacts.
        const T chord = step.tExit - step.tEntry;
        const T probability =
            f >= T(1) ? T(1) : T(1) - std::pow(T(1) - f, chord / delta);
        if (probability >= T(1) || uniform(rng) < probability) {
          result.cellId = lattice_->cellId(step.index);
          result.index = step.index;
          result.distance = step.tEntry;
          for (int d = 0; d < D; ++d)
            result.point[d] = origin[d] + direction[d] * step.tEntry;

          // Which face did it come in through? The axis that changed.
          Vec3D<T> faceNormal{0, 0, 0};
          if (havePrevious) {
            for (int d = 0; d < D; ++d)
              if (previous[d] != step.index[d]) {
                result.enteredAxis = d;
                result.enteredSign = previous[d] < step.index[d] ? -1 : 1;
                faceNormal[d] = static_cast<T>(result.enteredSign);
              }
          } else {
            // It interacted in the first cell it entered, so the face is the
            // one whose plane the ray crossed on the way in: the axis it is
            // travelling most steeply against.
            int axis = 0;
            T steepest = 0;
            for (int d = 0; d < D; ++d)
              if (std::abs(direction[d]) > steepest) {
                steepest = std::abs(direction[d]);
                axis = d;
              }
            result.enteredAxis = axis;
            result.enteredSign = direction[axis] > 0 ? -1 : 1;
            faceNormal[axis] = static_cast<T>(result.enteredSign);
          }

          result.normal = hitNormal(result.cellId, step.index, faceNormal);
          return false; // stop: the ray has met the surface
        }
      }
      previous = step.index;
      havePrevious = true;
      return true;
    });

    return result;
  }
};

/// Sets filling fractions from a signed distance, so a voxel geometry can
/// start with the sub-grid surface position a level set already knows.
///
/// A cell centred a distance phi from the surface (negative inside) is filled
/// to 0.5 - phi/delta, clamped. Without this, `fromLevelSets` gives every cell
/// a fraction of exactly one and the sub-voxel information is discarded before
/// the first step -- which would hand the voxel arm a worse initial condition
/// than the level-set arm, and make the comparison unfair from step zero.
template <class T, int D, class SignedDistance>
void fillFromSignedDistance(const LatticeMap<T, D> &lattice,
                            std::vector<T> &fill, SignedDistance &&phi) {
  const T delta = lattice.gridDelta();
  const auto &dims = lattice.dims();
  const auto &min = lattice.minCorner();

  size_t sites = 1;
  for (int d = 0; d < D; ++d)
    sites *= static_cast<size_t>(dims[d]);

  std::array<int, D> idx{};
  for (size_t flat = 0; flat < sites; ++flat) {
    size_t rem = flat;
    for (int d = 0; d < D; ++d) {
      idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
      rem /= static_cast<size_t>(dims[d]);
    }
    const int id = lattice.cellId(idx);
    if (id < 0)
      continue;
    Vec3D<T> centre{0, 0, 0};
    for (int d = 0; d < D; ++d)
      centre[d] = min[d] + delta * (static_cast<T>(idx[d]) + T(0.5));
    const T f = T(0.5) - phi(centre) / delta;
    fill[id] = std::min(T(1), std::max(T(0), f));
  }
}

} // namespace viennacs
