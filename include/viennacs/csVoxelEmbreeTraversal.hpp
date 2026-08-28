#pragma once

#include "csGridTraversal.hpp"

#include <rayUtil.hpp> // brings <embreeN/rtcore.h> behind the version guard

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#if VIENNARAY_EMBREE_VERSION < 4
#error "EmbreeCellTraversal is written against the embree4 query API"
#endif

namespace viennacs {

using namespace viennacore;

/// The cells of a voxel geometry as embree primitives.
///
/// The DDA in GridTraversal walks a ray through EVERY cell it crosses, which
/// makes the cost of a ray proportional to the depth of gas above the surface.
/// The level-set arm pays log(triangles) for the same question, so a wall-time
/// comparison between the two arms would measure the acceleration structure,
/// not the method. This closes that gap with the same library the level-set
/// arm uses: every cell holding material becomes one user primitive -- a box
/// -- in an embree BVH, and a ray finds its first interaction in one descent.
///
/// THE PHYSICS DOES NOT MOVE. The interaction rule stays exactly the one in
/// VoxelInteraction: a cell of fill f met over a chord L interacts with
/// probability 1-(1-f)^(L/delta). It lives in the intersect callback, where
/// the ray-box slabs hand over the chord for free; declining to report a hit
/// lets embree continue to the next candidate on its own, which is what the
/// DDA's "pass through" is.
///
/// ORDER INDEPENDENCE. Embree visits candidates in approximate, not strict,
/// near-to-far order, and may test one primitive twice. Both are harmless
/// under two provisions made here. First, each cell's acceptance is a HASH of
/// (ray seed, primitive), not a draw from a stream: retesting a cell gives
/// the same answer. Second, the nearest ACCEPTED candidate wins: with rolls
/// independent per cell, P(nearest accepted = k) = p_k * prod_{j<k}(1-p_j),
/// which is exactly the sequential transmission the DDA performs. A test with
/// binary fills, where no probability is rolled at all, must and does agree
/// with the DDA hit for hit.
template <class T, int D> class EmbreeCellTraversal {
  RTCScene scene_ = nullptr;
  const LatticeMap<T, D> *lattice_ = nullptr;
  const std::vector<T> *fill_ = nullptr;
  std::vector<std::array<int, D>> primIndex_; ///< primID -> lattice coordinate
  T delta_ = T(0);
  std::array<T, D> min_{};

  /// One device for the process. Creating a device spawns worker threads, and
  /// a scene is rebuilt every step; the device must not be.
  static RTCDevice device() {
    static RTCDevice dev = rtcNewDevice(nullptr);
    return dev;
  }

  static std::uint64_t mix(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
  }
  static double toUnit(std::uint64_t h) {
    return static_cast<double>(h >> 11) * 0x1.0p-53;
  }

  /// The extended query context: embree hands the base pointer back to the
  /// intersect callback, which casts it to this. The context is the whole
  /// per-ray state -- the double-precision ray, the gate, the seed, and the
  /// best accepted candidate so far, kept in double because embree's float
  /// tfar is used only to cull.
  struct Query {
    RTCRayQueryContext base; // must stay the first member
    const EmbreeCellTraversal *self = nullptr;
    std::array<double, D> org{}, dir{};
    double armAfter = 0;
    std::uint64_t seed = 0;
    double bestT = std::numeric_limits<double>::max();
    int bestPrim = -1;
    int bestAxis = -1;
    int bestSign = 0;
    bool bestInside = false;
  };

  static void boundsFunc(const RTCBoundsFunctionArguments *args) {
    const auto *self =
        static_cast<const EmbreeCellTraversal *>(args->geometryUserPtr);
    const auto &ii = self->primIndex_[args->primID];
    // Pad by a hair: the query ray embree culls with is float, the slab test
    // that decides is double, and a tangential ray must reach the test.
    const float pad = static_cast<float>(self->delta_) * 1e-4f;
    float lo[3], hi[3];
    for (int d = 0; d < D; ++d) {
      lo[d] = static_cast<float>(self->min_[d] +
                                 self->delta_ * static_cast<T>(ii[d])) -
              pad;
      hi[d] = static_cast<float>(self->min_[d] +
                                 self->delta_ * static_cast<T>(ii[d] + 1)) +
              pad;
    }
    if constexpr (D == 2) { // a 2D lattice is a one-cell-thick 3D one
      lo[2] = -static_cast<float>(self->delta_);
      hi[2] = static_cast<float>(self->delta_);
    }
    args->bounds_o->lower_x = lo[0];
    args->bounds_o->lower_y = lo[1];
    args->bounds_o->lower_z = lo[2];
    args->bounds_o->upper_x = hi[0];
    args->bounds_o->upper_y = hi[1];
    args->bounds_o->upper_z = hi[2];
  }

  static void intersectFunc(const RTCIntersectFunctionNArguments *args) {
    if (args->valid[0] == 0)
      return;
    auto *q = reinterpret_cast<Query *>(args->context);
    const auto *self = q->self;
    const unsigned prim = args->primID;
    const auto &ii = self->primIndex_[prim];

    // Ray-box slabs, in double and in D dimensions. The axis that closes the
    // entry interval is the face the ray came in through.
    double tEnter = -std::numeric_limits<double>::max();
    double tExit = std::numeric_limits<double>::max();
    int axis = -1;
    for (int d = 0; d < D; ++d) {
      const double lo =
          static_cast<double>(self->min_[d]) +
          static_cast<double>(self->delta_) * static_cast<double>(ii[d]);
      const double hi = lo + static_cast<double>(self->delta_);
      const double o = q->org[d], v = q->dir[d];
      if (std::abs(v) < 1e-15) {
        if (o < lo || o > hi)
          return;
        continue;
      }
      double t1 = (lo - o) / v, t2 = (hi - o) / v;
      if (t1 > t2)
        std::swap(t1, t2);
      if (t1 > tEnter) {
        tEnter = t1;
        axis = d;
      }
      if (t2 < tExit)
        tExit = t2;
    }
    if (tExit <= std::max(tEnter, 0.0))
      return; // missed, or entirely behind the origin
    if (tExit <= q->armAfter)
      return; // still within the interface it was emitted from

    const double tEff = std::max(tEnter, 0.0);
    if (tEff >= q->bestT)
      return; // a nearer candidate is already accepted

    const int cellId = self->lattice_->cellId(ii);
    const double f = static_cast<double>((*self->fill_)[cellId]);
    if (f < 1.0) {
      // The chord is measured from the origin for the cell that contains it,
      // exactly as the DDA measures it.
      const double chord = tExit - tEff;
      const double p =
          1.0 - std::pow(1.0 - f, chord / static_cast<double>(self->delta_));
      const double u = toUnit(mix(q->seed ^ mix(prim)));
      if (u >= p)
        return; // transmitted
    }

    q->bestT = tEff;
    q->bestPrim = static_cast<int>(prim);
    q->bestInside = tEnter < 0.0;
    if (!q->bestInside) {
      q->bestAxis = axis;
      q->bestSign = q->dir[axis] > 0 ? -1 : 1;
    } else {
      q->bestAxis = -1;
      q->bestSign = 0;
    }

    // Shrink the ray so embree culls beyond the accepted candidate -- to the
    // float JUST ABOVE tEff, so nothing at the same depth is culled by the
    // float rounding; ties and near-ties are settled by bestT in double.
    auto *rh = reinterpret_cast<RTCRayHit *>(args->rayhit);
    rh->ray.tfar = std::nextafterf(static_cast<float>(tEff),
                                   std::numeric_limits<float>::max());
    rh->hit.geomID = args->geomID;
    rh->hit.primID = prim;
  }

public:
  struct RawHit {
    int cellId = -1;
    std::array<int, D> index{};
    T tEntry = T(0);
    int axis = -1; ///< entry face axis, or -1 when the origin was inside
    int sign = 0;
    bool hit() const { return cellId >= 0; }
  };

  EmbreeCellTraversal() = default;
  EmbreeCellTraversal(const EmbreeCellTraversal &) = delete;
  EmbreeCellTraversal &operator=(const EmbreeCellTraversal &) = delete;
  ~EmbreeCellTraversal() { release(); }

  void release() {
    if (scene_ != nullptr) {
      rtcReleaseScene(scene_);
      scene_ = nullptr;
    }
  }

  bool built() const { return scene_ != nullptr; }

  /// Builds the BVH over every cell currently holding material. Fills change
  /// every step, so this is called once per trace, before the rays fly; the
  /// level-set arm remeshes its surface on the same cadence.
  void build(const LatticeMap<T, D> &lattice, const std::vector<T> &fill) {
    release();
    lattice_ = &lattice;
    fill_ = &fill;
    delta_ = lattice.gridDelta();
    min_ = lattice.minCorner();

    primIndex_.clear();
    const auto &dims = lattice.dims();
    size_t total = 1;
    for (int d = 0; d < D; ++d)
      total *= static_cast<size_t>(dims[d]);
    // Only the SURFACE BAND: a cell every ray must cross before any other
    // solid cell. A cell buried behind fully solid neighbours can never be a
    // first interaction, and carrying the bulk (hundreds of thousands of
    // cells here) deepens the tree that every ray descends. The band is any
    // solid cell with an opening -- fill below one, or a lattice hole --
    // anywhere in its 3^D neighbourhood: VERTEX neighbours included, so a ray
    // passing exactly through a staircase corner still meets a primitive and
    // cannot slip between face-adjacent band cells into the bulk.
    int span = 1;
    for (int d = 0; d < D; ++d)
      span *= 3;
    for (size_t s = 0; s < total; ++s) {
      size_t rem = s;
      std::array<int, D> ii{};
      for (int d = 0; d < D; ++d) {
        ii[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
        rem /= static_cast<size_t>(dims[d]);
      }
      const int id = lattice.cellId(ii);
      if (id < 0 || fill[id] <= T(0))
        continue;
      bool open = false;
      for (int n = 0; n < span && !open; ++n) {
        if (n == span / 2)
          continue; // the cell itself
        auto probe = ii;
        int rem2 = n;
        bool inside = true;
        for (int d = 0; d < D; ++d) {
          probe[d] += rem2 % 3 - 1;
          rem2 /= 3;
          if (probe[d] < 0 || probe[d] >= dims[d])
            inside = false;
        }
        if (!inside)
          continue; // beyond the lattice the material continues, as in
                    // fillClamped: an edge is a cut, not a surface
        const int nid = lattice.cellId(probe);
        if (nid < 0 || fill[nid] < T(1))
          open = true;
      }
      if (open)
        primIndex_.push_back(ii);
    }

    scene_ = rtcNewScene(device());
    rtcSetSceneBuildQuality(scene_, RTC_BUILD_QUALITY_LOW);
    RTCGeometry geom = rtcNewGeometry(device(), RTC_GEOMETRY_TYPE_USER);
    rtcSetGeometryUserPrimitiveCount(geom, primIndex_.size());
    rtcSetGeometryUserData(geom, this);
    rtcSetGeometryBoundsFunction(geom, &boundsFunc, nullptr);
    rtcSetGeometryIntersectFunction(geom, &intersectFunc);
    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene_, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene_);
  }

  /// The first cell the ray interacts with, by one BVH query. `seed` decides
  /// every acceptance on this ray segment; a re-emitted segment must bring a
  /// fresh seed, so the same cell is rolled anew after a bounce.
  RawHit firstHit(const std::array<T, D> &origin,
                  const std::array<T, D> &direction, std::uint64_t seed,
                  T armAfter = T(0)) const {
    Query q;
    rtcInitRayQueryContext(&q.base);
    q.self = this;
    for (int d = 0; d < D; ++d) {
      q.org[d] = static_cast<double>(origin[d]);
      q.dir[d] = static_cast<double>(direction[d]);
    }
    q.armAfter = static_cast<double>(armAfter);
    q.seed = seed;

    RTCRayHit rh{};
    rh.ray.org_x = static_cast<float>(origin[0]);
    rh.ray.org_y = static_cast<float>(origin[1]);
    rh.ray.org_z = D == 3 ? static_cast<float>(origin[2]) : 0.0f;
    rh.ray.dir_x = static_cast<float>(direction[0]);
    rh.ray.dir_y = static_cast<float>(direction[1]);
    rh.ray.dir_z = D == 3 ? static_cast<float>(direction[2]) : 0.0f;
    rh.ray.tnear = 0.0f;
    rh.ray.tfar = std::numeric_limits<float>::max();
    rh.ray.mask = 0xFFFFFFFFu;
    rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    RTCIntersectArguments args;
    rtcInitIntersectArguments(&args);
    args.context = &q.base;
    rtcIntersect1(scene_, &rh, &args);

    RawHit out;
    if (q.bestPrim >= 0) {
      out.index = primIndex_[static_cast<size_t>(q.bestPrim)];
      out.cellId = lattice_->cellId(out.index);
      out.tEntry = static_cast<T>(q.bestT);
      out.axis = q.bestAxis;
      out.sign = q.bestSign;
    }
    return out;
  }
};

} // namespace viennacs
