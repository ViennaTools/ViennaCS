#pragma once

#include "csVoxelInteraction.hpp"

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
  Step apply(std::vector<T> &fill, const std::vector<T> &velocity, T dt) const {
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

    std::array<int, D> idx{};
    auto unflatten = [&](size_t flat) {
      size_t rem = flat;
      for (int d = 0; d < D; ++d) {
        idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
        rem /= static_cast<size_t>(dims[d]);
      }
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
    for (size_t flat = 0; flat < sites; ++flat) {
      unflatten(flat);
      const int id = lattice_->cellId(idx);
      if (id < 0)
        continue;
      const T area = interfaceArea(fill, idx);
      if (area <= T(0))
        continue;
      const T volume = velocity[id] * area * dt;
      step.volumeRequested += volume;

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
            if (fill[nid] > T(0))
              target = nid;
          }
        }
        if (target < 0) {
          step.volumeRequested -= volume; // nowhere to put it: do not claim it
          continue;
        }
      }
      change[target] += volume / cellVolume;
    }

    for (size_t i = 0; i < fill.size(); ++i)
      fill[i] += change[i];

    // Redistribute what will not fit. A cell over one pushes the excess to the
    // neighbour the surface is advancing into; a cell under zero takes the
    // deficit from the neighbour behind it. Repeat, because that neighbour may
    // overflow in turn.
    for (int sweep = 0; sweep < 2 * D + 2; ++sweep) {
      bool moved = false;
      for (size_t flat = 0; flat < sites; ++flat) {
        unflatten(flat);
        const int id = lattice_->cellId(idx);
        if (id < 0)
          continue;
        const T f = fill[id];
        if (f <= T(1) + T(1e-15) && f >= -T(1e-15))
          continue;

        const bool overflowing = f > T(1);
        const T surplus = overflowing ? f - T(1) : f; // signed
        const auto g = fillGradient(fill, idx);
        // Outward, i.e. towards the void, is +grad since g = -grad(f).
        int axis = 0;
        T best = 0;
        for (int d = 0; d < D; ++d)
          if (std::abs(g[d]) > best) {
            best = std::abs(g[d]);
            axis = d;
          }
        int direction = 0;
        if (best > T(1e-12))
          direction = g[axis] > 0 ? 1 : -1;
        else
          direction = 1; // no interface to read: pick an axis and be consistent
        if (!overflowing)
          direction = -direction; // a deficit eats backwards, into the solid

        auto neighbour = idx;
        neighbour[axis] += direction;
        const int nid = lattice_->cellId(neighbour);
        if (nid < 0) {
          step.volumeLost += surplus * cellVolume;
          fill[id] = overflowing ? T(1) : T(0);
          moved = true;
          continue;
        }
        fill[id] = overflowing ? T(1) : T(0);
        fill[nid] += surplus;
        moved = true;
      }
      if (!moved)
        break;
    }

    for (size_t i = 0; i < fill.size(); ++i)
      step.volumeApplied += fill[i] * cellVolume;
    return step;
  }
};

} // namespace viennacs
