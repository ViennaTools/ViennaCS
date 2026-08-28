#pragma once

#include "csVoxelAdvance.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <rayReflection.hpp>
#include <rayUtil.hpp>

namespace viennacs {

using namespace viennacore;

/// Monte Carlo flux onto a voxel geometry.
///
/// Rays are launched from ViennaRay's own source and reflected by ViennaRay's
/// own laws, neither of which touches embree. That is deliberate: a comparison
/// between a level-set and a voxel simulation is only about the geometry if
/// everything else is shared, and the source distribution and the reflection
/// are everything else. What differs is how a ray finds the surface -- a walk
/// through cells here, a BVH descent over triangles there -- and what normal
/// it meets when it does.
///
/// NORMALISATION is the part to get right, and the part that decides whether
/// the two arms can be compared at all. A ray carries a rate, source flux
/// times source area over the ray count. A cell accumulates rates and divides
/// by the interface area it presents, which turns them back into the flux
/// density a rate law expects. On a flat surface facing an unobstructed
/// source, that has to give back the source flux exactly; if it does not, no
/// profile computed from it means anything.
template <class T, int D> class VoxelFlux {
  const LatticeMap<T, D> *lattice_ = nullptr;
  const std::vector<T> *fill_ = nullptr;
  VoxelInteraction<T, D> interaction_;
  VoxelAdvance<T, D> areas_;
  mutable std::vector<T> areaCache_; ///< per-trace, fills frozen while rays fly
  GridTraversal<T, D> traversal_;

public:
  VoxelFlux(const LatticeMap<T, D> &lattice, const std::vector<T> &fill,
            NormalEstimator estimator = NormalEstimator::FillGradientYoungs,
            AreaEstimator areaEstimator = AreaEstimator::Gradient)
      : lattice_(&lattice), fill_(&fill),
        interaction_(lattice, fill, estimator), areas_(lattice, areaEstimator),
        traversal_(lattice) {}

  void setTraversalEngine(TraversalEngine e) {
    interaction_.setTraversalEngine(e);
  }
  TraversalEngine traversalEngine() const {
    return interaction_.traversalEngine();
  }

  /// Rebuilds the acceleration structure AND the per-cell caches from the
  /// current fills. `trace` calls it itself; a caller driving
  /// `traceToSurface` directly -- the ion walker does -- must call it before
  /// its parallel region.
  ///
  /// The area cache exists because deposit() runs per ENCOUNTER -- tens of
  /// millions of times per step -- and each call asked two 3^D
  /// neighbourhoods for a 27-point area stencil that cannot change while the
  /// rays fly: the fills are frozen for the whole trace. Computing every
  /// cell's area once turns those stencils into reads of the same values.
  void prepareTransport() const {
    interaction_.prepare();
    const auto &dims = lattice_->dims();
    size_t sites = 1;
    for (int d = 0; d < D; ++d)
      sites *= static_cast<size_t>(dims[d]);
    areaCache_.assign(fill_->size(), T(0));
#pragma omp parallel for schedule(static)
    for (long long flat = 0; flat < static_cast<long long>(sites); ++flat) {
      std::array<int, D> idx{};
      size_t rem = static_cast<size_t>(flat);
      for (int d = 0; d < D; ++d) {
        idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
        rem /= static_cast<size_t>(dims[d]);
      }
      const int id = lattice_->cellId(idx);
      if (id >= 0)
        areaCache_[id] = areas_.interfaceArea(*fill_, idx);
    }
  }

  /// The cached interface area at a lattice coordinate; zero off the grid.
  /// Falls back to the direct stencil when no cache has been built, so a
  /// caller that never traced still gets the right answer.
  T areaAt(const std::array<int, D> &idx) const {
    const int id = lattice_->cellId(idx);
    if (id < 0)
      return T(0);
    if (areaCache_.size() == fill_->size())
      return areaCache_[id];
    return areas_.interfaceArea(*fill_, idx);
  }

  /// Credits `amount` to the interface at `idx`, spread over the local
  /// interface in proportion to the area each cell carries.
  ///
  /// Not to the single cell the ray interacted in, which is what a binary
  /// voxel code does and what this did at first. Across a smeared interface
  /// that is wrong: a ray meeting a cell of fill 0.3 passes through it seven
  /// times in ten, so the deeper cell collects more while carrying less area.
  /// On a blanket, where the answer is known, that gives 90 and 160 on one
  /// physical surface whose true density is 100, with the outermost interface
  /// cell collecting nothing at all -- a cell that would then refuse to move
  /// while its neighbour moved twice as fast.
  ///
  /// Giving each cell a share of the flux equal to its share of the area makes
  /// the density uniform across the smear by construction. The premise is the
  /// same as Youngs' stencil: the interface is a neighbourhood, not a cell.
  ///
  /// The cost is at a corner, where the neighbourhood straddles two surfaces
  /// and flux leaks between them over one cell. That is real, and it is the
  /// resolution limit rather than a choice: a one-cell neighbourhood cannot
  /// tell a corner from a smear.
  void deposit(std::vector<T> &collected, const std::array<int, D> &idx,
               T amount) const {
    int span = 1;
    for (int d = 0; d < D; ++d)
      span *= 3;

    T totalArea = 0;
    for (int s = 0; s < span; ++s) {
      std::array<int, D> probe = idx;
      int rem = s;
      for (int d = 0; d < D; ++d) {
        probe[d] += rem % 3 - 1;
        rem /= 3;
      }
      if (lattice_->cellId(probe) < 0)
        continue;
      totalArea += areaAt(probe);
    }
    if (totalArea <= T(0)) {
      const int id = lattice_->cellId(idx);
      if (id >= 0)
        collected[id] += amount;
      return;
    }
    for (int s = 0; s < span; ++s) {
      std::array<int, D> probe = idx;
      int rem = s;
      for (int d = 0; d < D; ++d) {
        probe[d] += rem % 3 - 1;
        rem /= 3;
      }
      const int id = lattice_->cellId(probe);
      if (id < 0)
        continue;
      const T area = areaAt(probe);
      if (area > T(0))
        collected[id] += amount * area / totalArea;
    }
  }

  /// Follows a ray until it meets the surface or leaves through the top.
  ///
  /// The lateral boundaries are reflective, as they are in a ViennaPS domain:
  /// a feature is one period of something repeating, not an island. Without
  /// this a ray launched obliquely near the edge simply walks out of the
  /// domain, and the surface measures less flux than was sent -- by a third in
  /// 2D and more than half in 3D, which looks exactly like a normalisation
  /// error and is in fact a missing boundary. ViennaRay's own boundary is
  /// embree geometry, so it is the one piece of its transport that cannot be
  /// borrowed here.
  template <class RNG>
  VoxelHit<T, D> traceToSurface(std::array<T, D> origin,
                                std::array<T, D> direction, RNG &rng,
                                T armAfter = T(0)) const {
    const T delta = lattice_->gridDelta();
    const auto &dims = lattice_->dims();
    const auto &minCorner = lattice_->minCorner();

    for (int crossing = 0; crossing < 64; ++crossing) {
      const auto hit =
          interaction_.firstHit(origin, direction, rng, armAfter);
      if (hit.hit())
        return hit;

      // No hit: find where it left, and mirror it if that was a side.
      T tNear = 0, tFar = 0;
      if (!traversal_.clip(origin, direction, tNear, tFar))
        return hit; // it never entered the grid at all

      std::array<T, D> exit{};
      for (int d = 0; d < D; ++d)
        exit[d] = origin[d] + direction[d] * tFar;

      int axis = -1;
      T closest = std::numeric_limits<T>::max();
      for (int d = 0; d < D - 1; ++d) { // lateral axes only
        const T lo = minCorner[d];
        const T hi = minCorner[d] + delta * static_cast<T>(dims[d]);
        const T distance = std::min(std::abs(exit[d] - lo), std::abs(exit[d] - hi));
        if (distance < closest) {
          closest = distance;
          axis = d;
        }
      }
      if (axis < 0 || closest > T(1e-6) * delta)
        return hit; // it left through the top or the bottom: gone

      direction[axis] = -direction[axis];
      for (int d = 0; d < D; ++d)
        origin[d] = exit[d];
      origin[axis] += direction[axis] * delta * T(1e-6);
    }
    return VoxelHit<T, D>{};
  }

  struct Result {
    std::vector<T> flux;    ///< per cell, in the units of `sourceFlux`
    size_t raysTraced = 0;
    size_t raysAbsorbed = 0; ///< rays that met the surface at least once
    size_t reemissions = 0;
  };

  /// Traces `numRays` and returns the flux density on every cell.
  ///
  /// `sticking` is the fraction of a ray's remaining rate deposited at each
  /// encounter; the rest is re-emitted diffusely about the local normal. A
  /// sticking of one gives line of sight, which is the limit where shadowing
  /// is strongest.
  /// `smoothingNeighbors` averages the flux over that many rings of the
  /// interface once it has been normalised, matching what ViennaPS does for
  /// the level-set arm, where RayTracingParameters::smoothingNeighbors is 1 by
  /// default and every engine applies it. Leaving it out here would not merely
  /// be noisier -- it would be a different transport convention, and the
  /// comparison is only about the geometry if everything else matches.
  ///
  /// It matters more than tidiness. Surface velocity is a nonlinear function
  /// of flux, and a growing front is unstable to noise in it: a cell that
  /// samples high grows faster, protrudes, and then genuinely collects more.
  Result trace(size_t numRays, T sourceFlux, T sticking,
               T cosinePower = T(1), unsigned seed = 1,
               T weightCutoff = T(1e-4), int smoothingNeighbors = 1) const {
    return trace(numRays, sourceFlux, std::vector<T>(fill_->size(), sticking),
                 cosinePower, seed, weightCutoff, smoothingNeighbors);
  }

  /// The same, with the sticking resolved PER CELL.
  ///
  /// A selective mechanism adsorbs differently on different materials, and the
  /// re-emission has to follow: a species that sticks on the substrate and
  /// reflects off the mask reaches deep into a feature, and one that sticks
  /// everywhere does not. A single value for the whole surface gets transport
  /// inside a feature wrong even where the surface solve is right, and a
  /// blanket of one material cannot show it.
  Result trace(size_t numRays, T sourceFlux, const std::vector<T> &sticking,
               T cosinePower = T(1), unsigned seed = 1,
               T weightCutoff = T(1e-4), int smoothingNeighbors = 1) const {
    prepareTransport(); // (re)build the BVH before any ray flies
    const auto &dims = lattice_->dims();
    const auto &minCorner = lattice_->minCorner();
    const T delta = lattice_->gridDelta();

    // The source spans the top of the lattice, looking down.
    //
    // Sampled here rather than through ViennaRay's SourceRandom, whose
    // interface differs between ViennaRay generations: one offers
    // getOriginAndDirection, another getOrigin and getDirection separately,
    // and a header compiling against one fails against the other. The law is
    // ViennaRay's own and unchanged -- a position uniform over the top face,
    // and cos(theta) drawn as r^(1/(power+1)) about the inward axis -- so both
    // arms still sample the same distribution. Only the coupling to an
    // interface that moves is gone.
    std::array<T, D> sourceMin{}, sourceMax{};
    for (int d = 0; d < D; ++d) {
      sourceMin[d] = minCorner[d];
      sourceMax[d] = minCorner[d] + delta * static_cast<T>(dims[d]);
    }
    T sourceArea = 1;
    for (int d = 0; d < D - 1; ++d)
      sourceArea *= (sourceMax[d] - sourceMin[d]);
    const T exponent = T(1) / (cosinePower + T(1));

    Result result;
    result.flux.assign(fill_->size(), T(0));
    result.raysTraced = numRays;

    const T rayRate = sourceFlux * sourceArea / static_cast<T>(numRays);

    std::vector<T> collected(fill_->size(), T(0));

    // Rays are independent, so they are traced in parallel -- as ViennaRay
    // traces the level-set arm's. Left serial, the voxel arm is not merely
    // slow, it is being compared against a parallel implementation, and any
    // statement about cost means nothing.
    //
    // Each thread accumulates into its own buffer and they are summed at the
    // end: contention on a shared buffer would cost more than the tracing.
    // Each also carries its own stream, so the result does not depend on how
    // many threads ran, only on the seed.
#pragma omp parallel
    {
      std::vector<T> mine(fill_->size(), T(0));
      const unsigned thread =
#ifdef _OPENMP
          static_cast<unsigned>(omp_get_thread_num());
#else
          0u;
#endif
      RNG rng(seed * 7919u + thread);

#pragma omp for schedule(static)
      for (long long r = 0; r < static_cast<long long>(numRays); ++r) {
      std::uniform_real_distribution<T> uni(T(0), T(1));
      std::array<T, D> origin{}, direction{};
      for (int d = 0; d < D - 1; ++d)
        origin[d] = sourceMin[d] + (sourceMax[d] - sourceMin[d]) * uni(rng);
      origin[D - 1] = sourceMax[D - 1];

      // Sample the direction in THREE dimensions and, in 2D, project it onto
      // the plane. That is what ViennaRay does, and it is not a detail: the
      // polar angle of the 3D cosine law is not the angle of the 2D one.
      //
      // Applying cos(theta) = u^(1/(n+1)) directly as a 2D angle gives a
      // cumulative of sin^2(theta) where 2D requires sin(theta) -- far too
      // collimated about the vertical. Down a trench of aspect ratio 1.5 that
      // put only sin^2(18.4 deg) = 0.10 of the field flux on the floor against
      // the 0.32 the geometry allows, so the floor received a third of what
      // line of sight alone would deliver. Projecting a 3D sample recovers the
      // 2D law exactly, for any cosine power, without a special case.
      const T cosTheta = std::pow(uni(rng), exponent);
      const T sinTheta = std::sqrt(std::max(T(0), T(1) - cosTheta * cosTheta));
      const T phi = T(2) * T(M_PI) * uni(rng);
      if constexpr (D == 2) {
        const T dx = std::cos(phi) * sinTheta;
        const T dy = -cosTheta;
        const T norm = std::sqrt(dx * dx + dy * dy);
        direction[0] = dx / norm;
        direction[1] = dy / norm;
      } else {
        direction[0] = std::cos(phi) * sinTheta;
        direction[1] = std::sin(phi) * sinTheta;
        direction[2] = -cosTheta;
      }

      T weight = rayRate;
      bool absorbed = false;
      for (int bounce = 0; bounce < 1000; ++bounce) {
        const auto hit = traceToSurface(origin, direction, rng);
        if (!hit.hit())
          break;
        absorbed = true;
        // The FULL weight, not weight * sticking. What a surface receives is
        // what arrives at it; whether it sticks is the rate law's business,
        // and it applies the sticking itself. Depositing the absorbed part
        // here instead would apply it twice. The sticking still governs how
        // much of the ray survives to carry on.
        deposit(mine, hit.index, weight);
        weight *= (T(1) - sticking[hit.cellId]);
        if (weight <= rayRate * weightCutoff)
          break;

        // Re-emit about the local normal, from OUTSIDE the interface.
        //
        // A fractional interface is two or three cells thick, so a ray nudged
        // out by a fraction of a cell is still inside it and can interact
        // again immediately -- and where the sticking is small it keeps its
        // full weight and deposits that full weight every time. On a binary
        // geometry this cannot happen: the ray leaves the one solid cell and
        // is gone. On a fractional one it inflated the measured flux from an
        // incident 1000 to 1600 over a few steps, and the growth rate with it.
        //
        // So the ray is placed past the last cell holding any material along
        // its outward normal. Having left the surface, it restarts outside it.
        Vec3D<T> normal3{0, 0, 0};
        for (int d = 0; d < D; ++d)
          normal3[d] = hit.normal[d];
        const auto newDir = viennaray::ReflectionDiffuse<T, D>(normal3, rng);

        // Restart just outside the interface it was emitted from. Three
        // schemes were measured -- this, a gate on reaching an empty cell, and
        // a distance epsilon -- and all three land within a few percent of one
        // another, so the restart is NOT what separates this arm from the
        // level-set one. This one is kept because it costs the blanket least.
        int axis = 0;
        T steepest = 0;
        for (int d = 0; d < D; ++d)
          if (std::abs(normal3[d]) > steepest) {
            steepest = std::abs(normal3[d]);
            axis = d;
          }
        const int outward = normal3[axis] > 0 ? 1 : -1;
        int clear = 1;
        auto probe = hit.index;
        for (int stepOut = 0; stepOut < 8; ++stepOut) {
          probe[axis] += outward;
          const int nid = lattice_->cellId(probe);
          if (nid < 0 || (*fill_)[nid] <= T(1e-9))
            break;
          ++clear;
        }
        for (int d = 0; d < D; ++d) {
          origin[d] = hit.point[d] +
                      normal3[d] * delta * (static_cast<T>(clear) + T(1e-3));
          direction[d] = newDir[d];
        }
      }
      if (absorbed)
#pragma omp atomic
        ++result.raysAbsorbed;
      }

#pragma omp critical
      for (size_t c = 0; c < mine.size(); ++c)
        collected[c] += mine[c];
    }

    // Rates become flux densities by dividing out the area each cell shows.
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
      const int id = lattice_->cellId(idx);
      if (id < 0 || collected[id] == T(0))
        continue;
      // Same guard as the chemistry applies: a sliver of interface would
      // divide a small rate by a smaller area and report an impossible flux.
      T faceArea = 1;
      for (int d = 0; d < D - 1; ++d)
        faceArea *= delta;
      const T area = areaAt(idx);
      if (area > T(1e-2) * faceArea)
        result.flux[id] = collected[id] / area;
    }

    smooth(result.flux, smoothingNeighbors);
    return result;
  }

  /// Area-weighted average of the flux over the interface neighbourhood.
  void smooth(std::vector<T> &flux, int rings) const {
    if (rings <= 0)
      return;
    const T delta = lattice_->gridDelta();
    T faceArea = 1;
    for (int d = 0; d < D - 1; ++d)
      faceArea *= delta;
    const T minArea = T(1e-2) * faceArea;

    const auto &dims = lattice_->dims();
    size_t sites = 1;
    for (int d = 0; d < D; ++d)
      sites *= static_cast<size_t>(dims[d]);
    const int width = 2 * rings + 1;
    int span = 1;
    for (int d = 0; d < D; ++d)
      span *= width;

    std::vector<T> smoothed = flux;
    std::array<int, D> idx{};
    for (size_t flat = 0; flat < sites; ++flat) {
      size_t rem = flat;
      for (int d = 0; d < D; ++d) {
        idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
        rem /= static_cast<size_t>(dims[d]);
      }
      const int id = lattice_->cellId(idx);
      if (id < 0)
        continue;
      if (areaAt(idx) <= minArea)
        continue;

      T sum = 0, weight = 0;
      for (int s = 0; s < span; ++s) {
        std::array<int, D> probe = idx;
        int r = s;
        for (int d = 0; d < D; ++d) {
          probe[d] += r % width - rings;
          r /= width;
        }
        const int nid = lattice_->cellId(probe);
        if (nid < 0)
          continue;
        const T area = areaAt(probe);
        if (area <= minArea)
          continue;
        sum += flux[nid] * area;
        weight += area;
      }
      if (weight > T(0))
        smoothed[id] = sum / weight;
    }
    flux.swap(smoothed);
  }
};

} // namespace viennacs
