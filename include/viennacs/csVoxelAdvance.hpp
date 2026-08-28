#pragma once

#include "csVoxelInteraction.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <utility>

namespace viennacs {

using namespace viennacore;

/// How much interface a cell presents. It enters the surface advance twice --
/// once turning a delivered particle rate into the flux density the chemistry
/// needs, and once turning a surface velocity back into a volume -- and the
/// two cancel only where the velocity is linear in the flux. With coverages
/// they do not, so the estimate matters.
///
///   StaircaseFaces  the exposed cell faces, (D-1)-dimensional, counted. It is
///                   the area the geometry literally has, and it is wrong: a
///                   surface at 45 degrees is a flight of steps whose total
///                   length exceeds the true interface by sqrt(2).
///
///   Gradient        |grad(f)| times the cell volume, which is the interfacial
///                   area per unit volume of a volume-fraction field. Uses the
///                   same Youngs stencil as the normal, so a mechanism sees an
///                   area and a normal that agree with each other.
enum class AreaEstimator { StaircaseFaces, Gradient };

/// Moves a voxel surface by changing filling fractions.
///
/// A cell whose interface moves at v_n gains v_n * A * dt of solid, so its
/// fraction changes by that over the cell volume. A fraction driven past one
/// or below zero is not clamped away -- the excess is handed to the neighbour
/// the surface is moving towards, because the material is real and a voxel
/// method's one structural advantage over a level set is that it need never
/// lose any. Clamping would throw exactly that away.
template <class T, int D> class VoxelAdvance {
  const LatticeMap<T, D> *lattice_ = nullptr;
  AreaEstimator areaEstimator_ = AreaEstimator::Gradient;
  bool wideStencil_ = true;

public:
  explicit VoxelAdvance(const LatticeMap<T, D> &lattice,
                        AreaEstimator estimator = AreaEstimator::Gradient)
      : lattice_(&lattice), areaEstimator_(estimator) {}

  void setAreaEstimator(AreaEstimator e) { areaEstimator_ = e; }

  T fillAt(const std::vector<T> &fill, const std::array<int, D> &idx) const {
    const int id = lattice_->cellId(idx);
    return id < 0 ? T(0) : fill[id];
  }

  /// The filling fraction for the purpose of a DERIVATIVE, with the lattice
  /// boundary treated as a cut through the material rather than as void.
  ///
  /// Reading zero off the lattice would put a phantom interface along every
  /// edge of the domain: the bottom row of a solid block would see void
  /// beneath it, find a gradient, and grow downwards out of the grid. A domain
  /// boundary is where the simulation stops, not where the material does, so
  /// the field continues with zero gradient across it.
  ///
  /// Note that `fillAt` above must NOT do this: a ray leaving the grid has to
  /// find nothing there, which is the opposite convention for the opposite
  /// question.
  T fillClamped(const std::vector<T> &fill,
                const std::array<int, D> &idx) const {
    auto inside = idx;
    for (int d = 0; d < D; ++d)
      inside[d] = std::min(std::max(inside[d], 0), lattice_->dims()[d] - 1);
    const int id = lattice_->cellId(inside);
    return id < 0 ? T(0) : fill[id];
  }

  /// -grad(f) in physical units, unnormalised: its magnitude is the
  /// interfacial area per unit volume, its direction the outward normal.
  Vec3D<T> fillGradient(const std::vector<T> &fill,
                        const std::array<int, D> &idx) const {
    const T delta = lattice_->gridDelta();
    Vec3D<T> g{0, 0, 0};

    int span = 1;
    for (int d = 0; d < D; ++d)
      span *= 3;
    // Youngs' (1,2,1) weighting across the axes not being differenced. The
    // normalisation is the total weight on one side, 2^(D-1), times the
    // two-cell span of the central difference.
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
      const T f = fillClamped(fill, probe);
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

  /// The interface area of one cell.
  T interfaceArea(const std::vector<T> &fill,
                  const std::array<int, D> &idx) const {
    const T delta = lattice_->gridDelta();
    T cellVolume = 1;
    for (int d = 0; d < D; ++d)
      cellVolume *= delta;

    if (areaEstimator_ == AreaEstimator::Gradient) {
      const auto g = fillGradient(fill, idx);
      T norm = 0;
      for (int d = 0; d < D; ++d)
        norm += g[d] * g[d];
      return std::sqrt(norm) * cellVolume;
    }

    // Staircase: count faces with solid on this side and not on the other.
    // Both halves of that matter. A cell that is not itself solid has no
    // interface, however empty its neighbours are -- without that test every
    // void cell in the grid claims 2D faces and grows.
    if (fillClamped(fill, idx) < T(0.5))
      return T(0);
    const T faceArea = cellVolume / delta;
    T area = 0;
    for (int d = 0; d < D; ++d)
      for (int sign = -1; sign <= 1; sign += 2) {
        auto neighbour = idx;
        neighbour[d] += sign;
        if (fillClamped(fill, neighbour) < T(0.5))
          area += faceArea;
      }
    return area;
  }

  /// Result of one step, for checking that nothing was lost.
  struct Step {
    T volumeRequested = 0; ///< sum of v_n * A * dt over all cells
    T volumeApplied = 0;   ///< what actually landed in the grid
    T volumeLost = 0;      ///< what had nowhere to go: off the lattice
  };

  /// Advances the surface by `dt`. `velocity` is indexed by cell id, positive
  /// for growth. Cells with no interface contribute nothing.
  /// `solid`, when given, makes the advance respect material boundaries:
  /// matter may not be moved between cells of DIFFERENT solids. Without it,
  /// fill is one substance, and at a two-material corner one material's
  /// over-etched deficit walks into the other and drains it -- a masked etch
  /// ate its mask foot at a thousand times the mask's own sputter rate, cell
  /// by cell, through exactly this. Cells labelled `wildcard` (gas: matter
  /// below the labelling threshold) match anything. Blocked dose is counted
  /// in volumeLost, never vanished.
  /// `velMat`, when given, names PER CELL the material the velocity was
  /// computed for. A gas skin cell at a two-material corner legitimately
  /// carries the floor's silicon etch -- but its inward placement walk can
  /// point laterally into the mask, and a wildcard-only guard lets the
  /// silicon-computed share drain mask matter. The sharper rule: a share may
  /// only land in the material it was computed for; a mismatch drops the
  /// share from volumeRequested, reported, never quietly delivered.
  /// The cells this step can possibly act on, as a mask over cell ids.
  ///
  /// Marked: every cell whose filling differs from a face neighbour (a
  /// lattice hole reads as empty, as fillAt reads it), dilated by TWO rings
  /// of the full 3^D neighbourhood. That is a provable superset of every
  /// cell any stencil here can see as interface: a non-uniform 3^D
  /// neighbourhood contains a face-differing pair (the block is
  /// face-connected), whose members the first pass marks, and one dilation
  /// reaches the centre; the second ring covers the cells a step's placement
  /// and spill can newly disturb. Everything OUTSIDE the mask is provably
  /// inert -- zero gradient, zero area, nothing to anchor -- so sweeping
  /// only the band changes no result, only the bill: the sweeps stop paying
  /// for the lattice and start paying for the surface.
  std::vector<unsigned char> interfaceBand(const std::vector<T> &fill) const {
    const auto &dims = lattice_->dims();
    size_t sites = 1;
    for (int d = 0; d < D; ++d)
      sites *= static_cast<size_t>(dims[d]);

    std::vector<unsigned char> band(fill.size(), 0);
    std::vector<std::array<int, D>> ring;

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
      const T f = fill[id];
      bool differs = false;
      for (int d = 0; d < D && !differs; ++d)
        for (int sgn = -1; sgn <= 1 && !differs; sgn += 2) {
          auto nb = idx;
          nb[d] += sgn;
          bool inside = nb[d] >= 0 && nb[d] < dims[d];
          const int nid = inside ? lattice_->cellId(nb) : -1;
          const T nf = nid < 0 ? (inside ? T(0) : f) : fill[nid];
          if (nf != f)
            differs = true;
        }
      if (differs && !band[id]) {
        band[id] = 1;
        ring.push_back(idx);
      }
    }

    int span = 1;
    for (int d = 0; d < D; ++d)
      span *= 3;
    for (int round = 0; round < 2; ++round) {
      std::vector<std::array<int, D>> next;
      for (const auto &c : ring)
        for (int n = 0; n < span; ++n) {
          auto nb = c;
          int rem = n;
          for (int d = 0; d < D; ++d) {
            nb[d] += rem % 3 - 1;
            rem /= 3;
          }
          const int nid = lattice_->cellId(nb);
          if (nid >= 0 && !band[nid]) {
            band[nid] = 1;
            next.push_back(nb);
          }
        }
      ring = std::move(next);
    }
    return band;
  }

  Step apply(std::vector<T> &fill, const std::vector<T> &velocity, T dt,
             const std::vector<int> *solid = nullptr,
             int wildcard = -1,
             const std::vector<int> *velMat = nullptr) const {
    auto sameSolid = [&](int a, int b) {
      if (!solid)
        return true;
      const int ma = (*solid)[a], mb = (*solid)[b];
      return ma == mb || ma == wildcard || mb == wildcard;
    };
    const T delta = lattice_->gridDelta();
    const auto &dims = lattice_->dims();
    T cellVolume = 1;
    for (int d = 0; d < D; ++d)
      cellVolume *= delta;

    size_t sites = 1;
    for (int d = 0; d < D; ++d)
      sites *= static_cast<size_t>(dims[d]);

    Step step;
    std::vector<T> change(fill.size(), T(0));
    const auto band = interfaceBand(fill);

    std::array<int, D> idx{};
    auto unflatten = [&](size_t flat) {
      size_t rem = flat;
      for (int d = 0; d < D; ++d) {
        idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
        rem /= static_cast<size_t>(dims[d]);
      }
    };
    auto flatten = [&](const std::array<int, D> &q) {
      size_t flat = 0;
      for (int d = D - 1; d >= 0; --d)
        flat = flat * static_cast<size_t>(dims[d]) + static_cast<size_t>(q[d]);
      return flat;
    };

    // WHERE the volume is computed and WHERE it is placed are separate
    // questions, and conflating them is what makes an algebraic volume-of-fluid
    // advance diffuse.
    //
    // The amount is right: |grad(f)| over the wide stencil recovers the
    // interface area exactly. But that stencil is nonzero over a band wider
    // than the interface, so placing each cell's share in the cell that
    // computed it deposits material AHEAD of the front, into cells that are
    // still empty. That tail then feeds on itself: one partial cell becomes
    // six in twenty steps, the apparent area grows with it, and a silane
    // blanket runs forty percent fast.
    //
    // So the volume is computed over the whole band and placed only in cells
    // that already hold material. A cell ahead of the front hands its share
    // back along the normal; the front reaches it when the cell behind fills
    // up and spills, which is what the redistribution below is for. The
    // interface then stays one cell wide however far it travels.
    // PARALLEL, with a deterministic reduction. Each thread walks one
    // contiguous block of sites (schedule(static) assigns blocks in thread
    // order) and records its placements and requested volume locally; the
    // blocks are merged in thread order afterwards. Contributions to any one
    // target cell therefore arrive in ascending site order -- the serial
    // order -- so the FILLS are bit-identical to the serial loop for any
    // thread count. Only the requested-volume scalar is summed in a different
    // grouping, which can differ in the last ulp; it is a diagnostic.
#ifdef _OPENMP
    const int numThreads = omp_get_max_threads();
#else
    const int numThreads = 1;
#endif
    std::vector<std::vector<std::pair<int, T>>> placedPer(numThreads);
    std::vector<T> requestedPer(numThreads, T(0));
#pragma omp parallel
    {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      auto &placed = placedPer[tid];
      T &requested = requestedPer[tid];
      std::array<int, D> idx{}; // shadows the shared scratch, deliberately
#pragma omp for schedule(static)
      for (long long sflat = 0; sflat < static_cast<long long>(sites);
           ++sflat) {
        size_t rem = static_cast<size_t>(sflat);
        for (int d = 0; d < D; ++d) {
          idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
          rem /= static_cast<size_t>(dims[d]);
        }
        const int id = lattice_->cellId(idx);
        if (id < 0)
          continue;
        if (velocity[id] == T(0))
          continue; // a zero velocity moves no volume wherever it would land
        const T area = interfaceArea(fill, idx);
        if (area <= T(0))
          continue;
        const T volume = velocity[id] * area * dt;
        requested += volume;

      // Placement is symmetric about the front, and the symmetry is not
      // cosmetic. A share computed at an EMPTY cell (deposition reaching ahead
      // of the front) walks inward to a cell holding material; a share
      // computed at a FULL cell (etching reaching behind the front) walks
      // outward to the partial front cell. Only the first half of this
      // existed at first, so deposition stayed one cell sharp while an etch
      // eroded the cell UNDER its front, widened the interface every step --
      // 560 surface cells became 1144 over twenty -- and the measured motion
      // decayed toward zero while the assigned velocities stayed correct. An
      // interface must be drained from its front for the same reason it must
      // be grown there.
      int target = id;
      if (fill[id] <= T(0)) {
        // Empty, and ahead of the front: walk back along the normal until a
        // cell with material is found.
        const auto g = fillGradient(fill, idx);
        int axis = 0;
        T best = 0;
        for (int d = 0; d < D; ++d)
          if (std::abs(g[d]) > best) {
            best = std::abs(g[d]);
            axis = d;
          }
        target = -1;
        if (best > T(1e-12)) {
          const int back = g[axis] > 0 ? -1 : 1; // opposite the outward normal
          auto probe = idx;
          for (int stepBack = 0; stepBack < 3 && target < 0; ++stepBack) {
            probe[axis] += back;
            const int nid = lattice_->cellId(probe);
            if (nid < 0)
              break;
            const bool matches =
                velMat ? ((*velMat)[id] == wildcard ||
                          (*solid)[nid] == wildcard ||
                          (*velMat)[id] == (*solid)[nid])
                       : sameSolid(id, nid);
            if (fill[nid] > T(0) && matches)
              target = nid;
          }
        }
        if (target < 0) {
          requested -= volume; // nowhere to put it: do not claim it
          continue;
        }
      } else if (fill[id] >= T(1)) {
        // Full, behind the front: hand the share to the partial front cell,
        // outward along the normal. If the interface is perfectly sharp there
        // is no partial cell and this cell IS the front: keep it.
        const auto g = fillGradient(fill, idx);
        int axis = 0;
        T best = 0;
        for (int d = 0; d < D; ++d)
          if (std::abs(g[d]) > best) {
            best = std::abs(g[d]);
            axis = d;
          }
        if (best > T(1e-12)) {
          const int outward = g[axis] > 0 ? 1 : -1;
          auto probe = idx;
          for (int stepOut = 0; stepOut < 3; ++stepOut) {
            probe[axis] += outward;
            const int nid = lattice_->cellId(probe);
            if (nid < 0 || fill[nid] <= T(0))
              break; // sharp interface: this cell is the front
            if (!sameSolid(id, nid))
              break; // a material boundary is a wall, not a conduit
            if (fill[nid] < T(1)) {
              target = nid;
              break;
            }
          }
        }
      }
        placed.emplace_back(target, volume / cellVolume);
      }
    }
    for (int t = 0; t < numThreads; ++t) {
      for (const auto &[cid, amount] : placedPer[t])
        change[cid] += amount;
      step.volumeRequested += requestedPer[t];
    }

    // The direction surplus travels must be read from the surface BEFORE it
    // is disturbed. Reading it afterwards asks the gradient of a field that
    // may hold values above one and that is degenerate at a cell just filled
    // to capacity; the axis then comes back lateral as often as not, and
    // surplus walks sideways into the neighbouring column instead of forward
    // into the void. On a flat blanket that showed up as columns differing by
    // thirteen cells of material while every fraction was a legal one.
    std::vector<int> pushAxis(fill.size(), -1);
    std::vector<int> pushDir(fill.size(), 0);
    // Independent per-cell writes: safe to parallelise, and deterministic.
#pragma omp parallel for schedule(dynamic, 256)
    for (long long pflat = 0; pflat < static_cast<long long>(sites); ++pflat) {
      std::array<int, D> pidx{};
      size_t rem = static_cast<size_t>(pflat);
      for (int d = 0; d < D; ++d) {
        pidx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
        rem /= static_cast<size_t>(dims[d]);
      }
      const int id = lattice_->cellId(pidx);
      if (id < 0)
        continue;
      if (!band[id])
        continue; // outside the band the gradient is zero and -1 stands
      const auto g = fillGradient(fill, pidx);
      int axis = 0;
      T best = 0;
      for (int d = 0; d < D; ++d)
        if (std::abs(g[d]) > best) {
          best = std::abs(g[d]);
          axis = d;
        }
      if (best > T(1e-12)) {
        pushAxis[id] = axis;
        pushDir[id] = g[axis] > 0 ? 1 : -1; // outward, towards the void
      }
    }

    // Applying the change also names the only cells that can now be out of
    // range: the redistribution below works through this list and what it
    // spills into, instead of rescanning the lattice once per sweep.
    std::vector<size_t> work;
    std::vector<unsigned char> queued(fill.size(), 0);
    for (size_t flat = 0; flat < sites; ++flat) {
      unflatten(flat);
      const int id = lattice_->cellId(idx);
      if (id < 0)
        continue;
      if (change[id] == T(0))
        continue;
      fill[id] += change[id];
      work.push_back(flat);
      queued[id] = 1;
    }

    // Redistribute what will not fit. A cell over one pushes the excess to the
    // neighbour the surface is advancing into; a cell under zero takes the
    // deficit from the neighbour behind it. Repeat, because that neighbour may
    // overflow in turn.
    //
    // The sweep count has to be generous, not a small constant. Surplus moves
    // one cell per sweep, so a step that deposits more than a couple of cells'
    // worth into one place -- which happens as soon as the velocity is not
    // uniform -- needs more passes than a fixed handful. Stopping early leaves
    // fractions above one, and a fraction above one is not a surface at all:
    // it is volume that the geometry cannot represent and that every later
    // measurement counts anyway. Silane reached 5.67 that way, and reported
    // twice the growth it had.
    const int maxSweeps = static_cast<int>(dims[0]) * 4;
    std::vector<size_t> next;
    for (int sweep = 0; sweep < maxSweeps && !work.empty(); ++sweep) {
      next.clear();
      for (const size_t flat : work) {
        unflatten(flat);
        const int id = lattice_->cellId(idx);
        queued[id] = 0;
        const T f = fill[id];
        if (f <= T(1) + T(1e-15) && f >= -T(1e-15))
          continue;

        const bool overflowing = f > T(1);
        const T surplus = overflowing ? f - T(1) : f; // signed
        int axis = pushAxis[id];
        int direction = pushDir[id];
        if (axis < 0) {
          // No interface here before the step: fall back to the axis the
          // surface as a whole is moving along, which for a grid built from a
          // level set is the last one.
          axis = D - 1;
          direction = 1;
        }
        if (!overflowing)
          direction = -direction; // a deficit eats backwards, into the solid

        auto neighbour = idx;
        neighbour[axis] += direction;
        const int nid = lattice_->cellId(neighbour);
        if (nid < 0 || !sameSolid(id, nid)) {
          step.volumeLost += surplus * cellVolume;
          fill[id] = overflowing ? T(1) : T(0);
          continue;
        }
        fill[id] = overflowing ? T(1) : T(0);
        fill[nid] += surplus;
        if (!queued[nid]) { // it may overflow in turn: next sweep's problem
          queued[nid] = 1;
          next.push_back(flatten(neighbour));
        }
      }
      std::swap(work, next);
    }

    // No partial cell may be fuller than every one of its neighbours. A cell
    // that is has come loose from the front: nothing anchors it, and it is
    // the droplet debris transport noise leaves a cell or two off the
    // surface. A sound interface is monotone -- each partial cell has a
    // neighbour at least as full, chaining back to a full cell -- so a
    // detached cell's content is merged into its fullest neighbour, toward
    // the front, and only matter with NOTHING around it is removed, counted
    // in volumeLost rather than vanished.
    for (int sweep = 0; sweep < 4; ++sweep) {
      bool moved = false;
      for (size_t flat = 0; flat < sites; ++flat) {
        unflatten(flat);
        const int id = lattice_->cellId(idx);
        if (id < 0)
          continue;
        // Everything a step writes lands inside the band, except a spill
        // chain's last recipient -- which its now-full donor anchors, so it
        // needs no visit either.
        if (!band[id])
          continue;
        const T f = fill[id];
        if (f <= T(1e-12) || f >= T(1) - T(1e-12))
          continue;
        T maxNbr = -1;
        int best = -1;
        for (int d = 0; d < D; ++d)
          for (int sgn = -1; sgn <= 1; sgn += 2) {
            auto nb = idx;
            nb[d] += sgn;
            const int nid = lattice_->cellId(nb);
            if (nid < 0) { // the boundary continues the field: anchored
              maxNbr = std::max(maxNbr, f);
              continue;
            }
            if (fill[nid] > maxNbr && sameSolid(id, nid)) {
              maxNbr = fill[nid];
              best = nid;
            }
          }
        if (f <= maxNbr + T(1e-12))
          continue; // anchored
        moved = true;
        if (maxNbr <= T(0) || best < 0) {
          step.volumeLost += f * cellVolume; // isolated debris, reported
          fill[id] = T(0);
          continue;
        }
        fill[best] += f;
        if (fill[best] > T(1)) {
          fill[id] = fill[best] - T(1); // remainder stays, now anchored
          fill[best] = T(1);
        } else
          fill[id] = T(0);
      }
      if (!moved)
        break;
    }

    // Nothing may be left out of range. If it is, the surplus had nowhere to
    // go and the geometry is not representable -- report it as lost rather
    // than let it stand as a fraction above one.
    for (size_t i = 0; i < fill.size(); ++i) {
      if (fill[i] > T(1)) {
        step.volumeLost += (fill[i] - T(1)) * cellVolume;
        fill[i] = T(1);
      } else if (fill[i] < T(0)) {
        step.volumeLost += fill[i] * cellVolume;
        fill[i] = T(0);
      }
      step.volumeApplied += fill[i] * cellVolume;
    }
    return step;
  }
};

} // namespace viennacs
