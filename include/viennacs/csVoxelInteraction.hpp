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
enum class NormalEstimator {
  Face,                ///< the face the ray entered through: 4 values in 2D
  FillGradient,        ///< -grad(f), two-point stencil
  FillGradientYoungs,  ///< -grad(f), Youngs' 3^D stencil
  InterfaceAverage,    ///< mean of the exposed cell faces within a radius
  InterfaceFit         ///< least-squares plane through the interface points
};

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

  /// The estimator's normal at a cell, outside a ray trace. Zero when the
  /// estimator has nothing to say (a Face normal is per-hit by nature).
  Vec3D<T> normalAt(const std::array<int, D> &idx) const {
    if (estimator_ == NormalEstimator::Face)
      return Vec3D<T>{0, 0, 0};
    const int id = lattice_->cellId(idx);
    if (id < 0)
      return Vec3D<T>{0, 0, 0};
    if (normalValid_.size() == fill_->size())
      return normalValid_[id] ? normalCache_[id] : Vec3D<T>{0, 0, 0};
    return (estimator_ == NormalEstimator::InterfaceAverage ||
            estimator_ == NormalEstimator::InterfaceFit)
               ? (estimator_ == NormalEstimator::InterfaceFit ? interfaceFitNormal(idx) : interfaceNormal(idx))
               : gradientNormal(idx, Vec3D<T>{0, 0, 0},
                                estimator_ == NormalEstimator::FillGradientYoungs);
  }

  /// Normal AND centroid of the least-squares plane at a cell, in cell
  /// units: the reconstructed continuous surface there.
  bool fitPlaneAt(const std::array<int, D> &idx, Vec3D<T> &normal,
                  std::array<T, D> &centroid) const {
    normal = interfaceFitNormal(idx);
    centroid = lastFitCentroid_;
    return normal[0] != T(0) || normal[1] != T(0) || normal[2] != T(0);
  }

  void setInterfaceRadius(int r) const { interfaceRadius_ = r > 0 ? r : 1; }
  /// Disable to get MCFPM's plain fixed-radius fit.
  void setCurvatureCap(bool on) { curvatureCap_ = on; }
  /// Off restores the plain fit over the whole neighbourhood.
  void setSegmentedFit(bool on) { segment_ = on; }
  void setSegmentTolerance(T t) { segmentTol_ = t > T(0) ? t : T(1.5); }
  bool curvatureCap() const { return curvatureCap_; }
  void setMinInterfaceRadius(int r) { minInterfaceRadius_ = r > 1 ? r : 2; }
  void setCurvatureAlpha(T a) { curvAlpha_ = a > T(0) ? a : T(0.18); }
  /// fits, of those capped, and the mean radius the cap chose
  void capStats(long &fits, long &fired, double &meanR) const {
    fits = capFits_; fired = capFired_;
    meanR = capFired_ ? double(capRSum_) / double(capFired_) : 0.0; }
  T lastFitRms() const { return lastFitRms_; }
  /// how many surface cells the last plane fit was built from
  int lastFitPoints() const { return lastFitPts_; }
  int interfaceRadius() const { return interfaceRadius_; }

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
      const auto n = (estimator_ == NormalEstimator::InterfaceAverage ||
                          estimator_ == NormalEstimator::InterfaceFit)
                         ? (estimator_ == NormalEstimator::InterfaceFit ? interfaceFitNormal(idx) : interfaceNormal(idx))
                         : gradientNormal(idx, Vec3D<T>{0, 0, 0}, wide);
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
    if (estimator_ == NormalEstimator::InterfaceAverage ||
        estimator_ == NormalEstimator::InterfaceFit) {
      const auto n = interfaceNormal(idx);
      return (n[0] != T(0) || n[1] != T(0) || n[2] != T(0)) ? n : faceNormal;
    }
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

  /// Radius, in cells, of the InterfaceAverage stencil. A binary field
  /// carries its orientation in the ARRANGEMENT of cells, not in any one
  /// cell, so the stencil has to be wider than the 3^D a gradient uses:
  /// measured against planes voxelised at known tilts, the mean error is
  /// 45 deg for a face normal, 11 deg for Youngs' 3^D, 3.1 deg at R = 2 and
  /// 2.5 deg at R = 3.
  mutable int interfaceRadius_ = 3;
  mutable T lastFitRms_ = 0;     ///< RMS plane residual of the last fit, cells
  mutable int lastFitPts_ = 0;   ///< surface cells the last fit used
  T curvAlpha_ = T(0.18);        ///< allowed rms per unit radius
  mutable long capFits_ = 0, capFired_ = 0, capRSum_ = 0;
  int minInterfaceRadius_ = 2;   ///< the cap never shrinks below this
  bool curvatureCap_ = true;     ///< shrink the stencil on curved surfaces
  bool segment_ = true;          ///< fit one face only, not across a corner
  int seedRadius_ = 2;           ///< stencil for the segmentation seed
  T segmentTol_ = T(1.5);        ///< cells off the seed plane, kept

  /// Mean of the outward unit normals of every solid/gas cell face within
  /// `interfaceRadius_` cells. This is the staircase's own surface, averaged
  /// -- the orientation a facetted geometry actually presents, rather than
  /// the one facet a ray happened to enter through. Zero when the
  /// neighbourhood holds no interface at all, so the caller can fall back.
  /// Least-squares plane through the solid/gas face midpoints in the
  /// neighbourhood -- the estimator the Monte Carlo feature profile models
  /// use, developed there for specular ion scattering off a cubic mesh.
  ///
  /// Measured against planes voxelised at known tilts, mean error over
  /// 0-90 deg: face normal 45 deg, Youngs' 3^D stencil 11 deg, mean of the
  /// exposed faces 2.5 deg, this fit 1.4 deg. The eigenvector of the smallest
  /// eigenvalue of the point covariance is the plane normal; the mean face
  /// normal only fixes its sign.
  /// The fitted plane's centroid, in CELL units -- the point the plane
  /// passes through. Paired with interfaceFitNormal it is the reconstructed
  /// surface MCFPM intersects the ray with to get the exact impact point.
  mutable std::array<T, D> lastFitCentroid_{};

  /// The same fit at a chosen radius, without disturbing the cached one.
  /// Used for the SPECULAR direction, which needs a smoother normal than the
  /// yield angle does.
  Vec3D<T> fitNormalAt(const std::array<int, D> &idx, int radius) const {
    const int keep = interfaceRadius_;
    interfaceRadius_ = radius > 1 ? radius : 2;
    const auto n = interfaceFitNormal(idx);
    interfaceRadius_ = keep;
    return n;
  }

  Vec3D<T> fitNormalCore(const std::array<int, D> &idx,
                         const Vec3D<T> *seed = nullptr) const {
    // MCFPM's algorithm (Huard thesis 2.5.3; Guo & Sawin 2009 review): fit a
    // plane Ax+By+Cz=D by least squares to the CENTRES of the surface sites
    // within a search distance, typically 4*dx. Cell centres, not cell faces:
    // "a continuous surface must be generated from the surface cells rather
    // than using the faces of the surface cells. The use of the cell surface
    // can give rise to numerical artefacts." D is their centre of mass, so
    // the plane through the centroid with the least-variance direction as its
    // normal -- the smallest-eigenvalue eigenvector of the covariance.
    const int R = interfaceRadius_;
    std::vector<std::array<T, D>> pts;
    std::array<T, D> mean{};
    std::array<int, D> at{}, lo{}, hi{};
    for (int d = 0; d < D; ++d) { lo[d] = idx[d] - R; hi[d] = idx[d] + R; at[d] = lo[d]; }
    while (true) {
      if (fillAt(at) >= T(0.5)) {
        bool exposed = false;              // "one or more faces exposed"
        for (int d = 0; d < D && !exposed; ++d)
          for (int sgn = -1; sgn <= 1; sgn += 2) {
            auto nb = at; nb[d] += sgn;
            if (fillAt(nb) < T(0.5)) { exposed = true; break; }
          }
        if (exposed) {
          std::array<T, D> p{};
          for (int d = 0; d < D; ++d) p[d] = static_cast<T>(at[d]);
          // SEGMENTATION: keep only the cells that lie on the SAME face as
          // the centre cell. A plane fit is meaningful over one face; at a
          // floor/wall corner the neighbourhood holds both, and the wall wins
          // on cell count -- a tall mask contributes a stacked column of
          // exposed cells within the radius while the floor contributes a
          // single row -- so the fitted normal swings toward horizontal and
          // an ion reads grazing incidence on a flat floor. Measured on a
          // W=40 trench, the share of floor impacts binned at 80-90 deg runs
          // 1.1 % at R=4 and 9.6 % at R=12, and the yield in that bin falls
          // from 42.7 to 0.7 atoms per ion.
          bool keep = true;
          if (seed) {
            T off = 0;
            for (int d = 0; d < D; ++d)
              off += (*seed)[d] * (p[d] - static_cast<T>(idx[d]));
            keep = std::abs(off) <= segmentTol_;
          }
          if (keep) {
            for (int d = 0; d < D; ++d) mean[d] += p[d];
            pts.push_back(p);
          }
        }
      }
      int d = 0;
      for (; d < D; ++d) { if (++at[d] <= hi[d]) break; at[d] = lo[d]; }
      if (d == D) break;
    }
    const int m = static_cast<int>(pts.size());
    lastFitPts_ = m;
    if (m < D + 1)
      return Vec3D<T>{0, 0, 0};
    for (int d = 0; d < D; ++d) mean[d] /= static_cast<T>(m);
    lastFitCentroid_ = mean;
    T C[D][D] = {};
    for (const auto &p : pts)
      for (int a = 0; a < D; ++a)
        for (int b = 0; b < D; ++b)
          C[a][b] += (p[a] - mean[a]) * (p[b] - mean[b]);
    T V[D][D] = {};
    for (int d = 0; d < D; ++d) V[d][d] = 1;
    for (int sweep = 0; sweep < 32; ++sweep) {
      int pI = 0, qI = 1; T big = 0;
      for (int a = 0; a < D; ++a)
        for (int b = a + 1; b < D; ++b)
          if (std::abs(C[a][b]) > big) { big = std::abs(C[a][b]); pI = a; qI = b; }
      if (big < T(1e-12)) break;
      const T th = T(0.5) * std::atan2(T(2) * C[pI][qI], C[pI][pI] - C[qI][qI]);
      const T c = std::cos(th), sn = std::sin(th);
      for (int k = 0; k < D; ++k) {
        const T cp = C[pI][k], cq = C[qI][k];
        C[pI][k] = c * cp + sn * cq; C[qI][k] = -sn * cp + c * cq;
      }
      for (int k = 0; k < D; ++k) {
        const T cp = C[k][pI], cq = C[k][qI];
        C[k][pI] = c * cp + sn * cq; C[k][qI] = -sn * cp + c * cq;
        const T vp = V[k][pI], vq = V[k][qI];
        V[k][pI] = c * vp + sn * vq; V[k][qI] = -sn * vp + c * vq;
      }
    }
    int least = 0;
    for (int d = 1; d < D; ++d) if (C[d][d] < C[least][least]) least = d;
    // RMS distance of the points from the fitted plane, in CELL units. On a
    // flat wall this is just the staircase, ~0.29; on a curved one it grows
    // as the patch departs from planarity, which is what caps the stencil.
    lastFitRms_ = std::sqrt(std::max(T(0), C[least][least]) / static_cast<T>(m));
    Vec3D<T> n{0, 0, 0};
    for (int d = 0; d < D; ++d) n[d] = V[d][least];
    // Orientation by polling, as MCFPM does: three points along +n and -n;
    // a solid cell votes against that direction, a gas cell votes for it.
    int votes = 0;
    for (int step = 1; step <= 3; ++step)
      for (int sgn = -1; sgn <= 1; sgn += 2) {
        std::array<int, D> probe{};
        for (int d = 0; d < D; ++d)
          probe[d] = idx[d] + static_cast<int>(std::lround(sgn * step * n[d]));
        const int v = (fillAt(probe) < T(0.5)) ? 1 : -1;
        votes += sgn * v;
      }
    if (votes < 0) for (int d = 0; d < D; ++d) n[d] = -n[d];
    T len = 0;
    for (int d = 0; d < D; ++d) len += n[d] * n[d];
    if (len < T(1e-12)) return Vec3D<T>{0, 0, 0};
    len = std::sqrt(len);
    for (int d = 0; d < D; ++d) n[d] /= len;
    return n;
  }


  /// The plane fit, with the stencil capped by the LOCAL CURVATURE.
  ///
  /// A least-squares plane is only meaningful over a patch that is actually
  /// planar. On a trench sidewall that is any radius, so the widest stencil
  /// wins and the staircase is smoothed away. On a hole of bore radius r it
  /// is not: once the patch spans a large angle the fit cuts across the
  /// bore, the normal swings inward, and the etch collapses.
  ///
  /// Measured on masked holes, the floor error tracks R/r rather than R --
  /// r=8/R=4 and r=20/R=12 agree at -16.9 % and -16.1 %, while r=8/R=8
  /// (R/r = 1) falls to -33 %. So the bound is angular. For a patch of
  /// half-width R on radius r the out-of-plane RMS is ~0.298 R^2 / r, which
  /// turns the empirical R/r <= 0.6 into the scale-free test
  ///
  ///     rms <= curvAlpha_ * R,    kCurvAlpha = 0.18
  ///
  /// costing nothing to evaluate: rms is the smallest eigenvalue of the fit
  /// we already do. A flat staircase gives rms ~0.29, which passes at every
  /// R >= 2, so this never fires on a trench. When it does fire, inverting
  /// the same relation gives the radius that would have passed,
  /// R' = curvAlpha_ * R^2 / rms, and one refit there is enough.
  Vec3D<T> interfaceFitNormal(const std::array<int, D> &idx) const {
    const int R = interfaceRadius_;
    // Seed the segmentation with the mean of the exposed faces over a small
    // stencil: crude (2.5 deg) but it cannot straddle a corner, which is all
    // that is asked of it. The full fit then runs on one face only.
    Vec3D<T> n;
    if (segment_ && R > seedRadius_) {
      const int keepR = interfaceRadius_;
      interfaceRadius_ = seedRadius_;
      const Vec3D<T> n0 = interfaceNormal(idx);
      interfaceRadius_ = keepR;
      T l = 0;
      for (int d = 0; d < D; ++d) l += n0[d] * n0[d];
      n = (l > T(1e-12)) ? fitNormalCore(idx, &n0) : fitNormalCore(idx);
    } else {
      n = fitNormalCore(idx);
    }
    ++capFits_;
    if (!curvatureCap_ || R <= minInterfaceRadius_ || lastFitRms_ <= T(0))
      return n;
    if (lastFitRms_ <= curvAlpha_ * static_cast<T>(R))
      return n;
    int rNew = static_cast<int>(curvAlpha_ * static_cast<T>(R) *
                                static_cast<T>(R) / lastFitRms_);
    if (rNew < minInterfaceRadius_) rNew = minInterfaceRadius_;
    if (rNew >= R) return n;
    ++capFired_; capRSum_ += rNew;
    interfaceRadius_ = rNew;
    const auto nn = fitNormalCore(idx);
    interfaceRadius_ = R;
    // A degenerate neighbourhood at the smaller radius is worse than a
    // curved fit at the larger one, so keep the wide answer in that case.
    T len = 0;
    for (int d = 0; d < D; ++d) len += nn[d] * nn[d];
    return len > T(0.5) ? nn : n;
  }

  Vec3D<T> interfaceNormal(const std::array<int, D> &idx) const {
    const int R = interfaceRadius_;
    Vec3D<T> sum{0, 0, 0};
    std::array<int, D> at{}, nb{};
    std::array<int, D> lo{}, hi{};
    for (int d = 0; d < D; ++d) {
      lo[d] = idx[d] - R;
      hi[d] = idx[d] + R;
      at[d] = lo[d];
    }
    while (true) {
      if (fillAt(at) >= T(0.5)) {
        for (int d = 0; d < D; ++d)
          for (int sgn = -1; sgn <= 1; sgn += 2) {
            nb = at;
            nb[d] += sgn;
            if (fillAt(nb) < T(0.5))
              sum[d] += static_cast<T>(sgn);   // face points solid -> gas
          }
      }
      int d = 0;
      for (; d < D; ++d) {
        if (++at[d] <= hi[d])
          break;
        at[d] = lo[d];
      }
      if (d == D)
        break;
    }
    const T len = std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
    if (len < T(1e-12))
      return Vec3D<T>{0, 0, 0};
    return Vec3D<T>{sum[0] / len, sum[1] / len, sum[2] / len};
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
