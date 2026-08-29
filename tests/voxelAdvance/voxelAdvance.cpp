#include <csVoxelAdvance.hpp>

#include <lsMakeGeometry.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>

// Moving a surface by changing filling fractions.
//
// Three things have to be true, and the third is the one that decides whether
// the area estimate matters.
//
//   1. Nothing is lost. Mass conservation is the voxel method's structural
//      advantage over a level set, whose advection conserves volume only to
//      the accuracy of the scheme. If material vanishes at a clamp, that
//      advantage is gone and there is little reason to prefer voxels.
//
//   2. A flat surface under a uniform velocity moves at that velocity.
//
//   3. So does a TILTED one. This is where the staircase area is wrong: the
//      steps of a 45 degree surface are longer than the surface, so a cell
//      believes it has more interface than it has and advances too fast. The
//      error is geometric, not statistical, and it does not average out.

namespace ls = viennals;
namespace cs = viennacs;

using T = double;

int failures = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
  if (!detail.empty())
    std::cout << "   " << detail;
  std::cout << "\n";
  if (!ok)
    ++failures;
}

template <int D> cs::DenseCellSet<T, D> makeBlock(T gridDelta, T lo, T hi) {
  ls::BoundaryConditionEnum bc[D];
  for (int i = 0; i < D - 1; ++i)
    bc[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  bc[D - 1] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;
  T bounds[2 * D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = lo;
    bounds[2 * i + 1] = hi;
  }
  T origin[D] = {}, normal[D] = {};
  origin[D - 1] = hi;
  normal[D - 1] = 1.;
  auto plane = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bc, gridDelta);
  ls::MakeGeometry<T, D>(plane,
                         ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
      .apply();
  std::vector<ls::SmartPointer<ls::Domain<T, D>>> lss{plane};
  cs::DenseCellSet<T, D> cellSet;
  cellSet.setCellSetPosition(false);
  cellSet.setCoverMaterial(0);
  cellSet.fromLevelSets(lss, nullptr, lo);
  return cellSet;
}

/// Total solid in the grid.
template <int D> T totalVolume(const std::vector<T> &fill, T gridDelta) {
  T cellVolume = 1;
  for (int d = 0; d < D; ++d)
    cellVolume *= gridDelta;
  T sum = 0;
  for (const auto f : fill)
    sum += f;
  return sum * cellVolume;
}

void checkConservation() {
  constexpr int D = 2;
  const T gridDelta = 1.0;
  auto cellSet = makeBlock<D>(gridDelta, -15., 15.);
  cs::LatticeMap<T, D> lattice(cellSet);
  std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
  cs::fillFromSignedDistance<T, D>(
      lattice, fill, [](const viennacore::Vec3D<T> &p) { return p[1] - T(0); });

  cs::VoxelAdvance<T, D> advance(lattice);
  std::vector<T> velocity(fill.size(), T(0.2)); // uniform growth

  T before = totalVolume<D>(fill, gridDelta);
  T requested = 0, lost = 0;
  for (int s = 0; s < 10; ++s) {
    const auto step = advance.apply(fill, velocity, 0.5);
    requested += step.volumeRequested;
    lost += step.volumeLost;
  }
  const T after = totalVolume<D>(fill, gridDelta);

  check("every requested volume lands in the grid",
        std::abs((after - before) - requested) < 1e-9,
        "gained " + std::to_string(after - before) + " of " +
            std::to_string(requested) + " requested");
  check("nothing runs off the lattice", std::abs(lost) < 1e-12,
        "lost " + std::to_string(lost));
}

/// The position of a level-set-free surface: where the filling fractions sum
/// to, which for a monotone profile is the interface height.
T surfaceHeight(const std::vector<T> &fill, const cs::LatticeMap<T, 2> &lattice,
                int column) {
  const T delta = lattice.gridDelta();
  const auto &dims = lattice.dims();
  T height = lattice.minCorner()[1];
  for (int j = 0; j < dims[1]; ++j) {
    const int id = lattice.cellId({column, j});
    if (id >= 0)
      height += fill[id] * delta;
  }
  return height;
}

void checkFlatAdvance() {
  constexpr int D = 2;
  const T gridDelta = 1.0;
  auto cellSet = makeBlock<D>(gridDelta, -15., 15.);
  cs::LatticeMap<T, D> lattice(cellSet);
  std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
  cs::fillFromSignedDistance<T, D>(
      lattice, fill, [](const viennacore::Vec3D<T> &p) { return p[1] - T(0); });

  cs::VoxelAdvance<T, D> advance(lattice);
  std::vector<T> velocity(fill.size(), T(0.3));

  const int column = lattice.dims()[0] / 2;
  const T start = surfaceHeight(fill, lattice, column);
  const T dt = 0.5;
  const int steps = 12;
  for (int s = 0; s < steps; ++s)
    advance.apply(fill, velocity, dt);
  const T moved = surfaceHeight(fill, lattice, column) - start;
  const T expected = 0.3 * dt * steps;

  // Not exact: |grad(f)| over a discrete stencil recovers the interface area
  // to a fraction of a percent, not to machine precision, and the surface
  // velocity inherits that. A percent is the honest tolerance.
  check("a flat surface moves at its velocity",
        std::abs(moved - expected) < 0.01 * expected,
        "moved " + std::to_string(moved) + ", expected " +
            std::to_string(expected) + "  (" +
            std::to_string(100 * (moved - expected) / expected) + "%)");
}

void checkTiltedAdvance() {
  constexpr int D = 2;
  const T gridDelta = 1.0;

  std::cout << "     tilt     gradient area     staircase area   (expected "
               "3.000)\n";
  T worstGradient = 0, worstStaircase = 0;
  for (int deg = 0; deg <= 45; deg += 15) {
    const T angle = deg * M_PI / 180.;
    const viennacore::Vec3D<T> n{std::sin(angle), std::cos(angle), 0};

    auto run = [&](cs::AreaEstimator estimator) {
      auto cellSet = makeBlock<D>(gridDelta, -25., 25.);
      cs::LatticeMap<T, D> lattice(cellSet);
      std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
      cs::fillFromSignedDistance<T, D>(lattice, fill,
                                       [&](const viennacore::Vec3D<T> &p) {
                                         return p[0] * n[0] + p[1] * n[1];
                                       });
      cs::VoxelAdvance<T, D> advance(lattice, estimator);
      std::vector<T> velocity(fill.size(), T(0.25));

      const T before = totalVolume<D>(fill, gridDelta);
      const T dt = 0.5, steps = 24;
      for (int s = 0; s < steps; ++s)
        advance.apply(fill, velocity, dt);
      const T after = totalVolume<D>(fill, gridDelta);

      // A surface of length L advancing a distance x sweeps L*x of area. The
      // interface spans the domain width W at a tilt, so L = W / cos(angle),
      // and the distance moved follows from the volume gained.
      const T width = 50.0;
      const T length = width / std::cos(angle);
      return (after - before) / length;
    };

    const T grad = run(cs::AreaEstimator::Gradient);
    const T stair = run(cs::AreaEstimator::StaircaseFaces);
    const T expected = 0.25 * 0.5 * 24;
    worstGradient = std::max(worstGradient, std::abs(grad - expected));
    worstStaircase = std::max(worstStaircase, std::abs(stair - expected));
    std::cout << "     " << std::setw(4) << deg << " deg" << std::fixed
              << std::setprecision(3) << std::setw(14) << grad << std::setw(19)
              << stair << "\n";
  }

  check("a tilted surface advances at its velocity, with the gradient area",
        worstGradient < 0.25,
        "worst error " + std::to_string(worstGradient) + " of 3.0");
  check("the staircase area advances it too fast, as the geometry predicts",
        worstStaircase > worstGradient,
        "staircase off by " + std::to_string(worstStaircase) + " vs gradient " +
            std::to_string(worstGradient));
}

void checkReversibility() {
  constexpr int D = 2;
  const T gridDelta = 1.0;
  auto cellSet = makeBlock<D>(gridDelta, -15., 15.);
  cs::LatticeMap<T, D> lattice(cellSet);
  std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
  cs::fillFromSignedDistance<T, D>(
      lattice, fill, [](const viennacore::Vec3D<T> &p) { return p[1] - T(0); });
  const T before = totalVolume<D>(fill, gridDelta);

  cs::VoxelAdvance<T, D> advance(lattice);
  std::vector<T> grow(fill.size(), T(0.2)), etch(fill.size(), T(-0.2));
  T moved = 0;
  for (int s = 0; s < 6; ++s)
    moved += std::abs(advance.apply(fill, grow, 0.5).volumeRequested);
  for (int s = 0; s < 6; ++s)
    moved += std::abs(advance.apply(fill, etch, 0.5).volumeRequested);

  // Not exact, and it should not be expected to be: the interface area a cell
  // presents changes as the surface moves through it, so the volume removed on
  // the way back is not identical to the volume added on the way out. What
  // must not happen is material disappearing at a clamp, which would show up
  // far above the area estimate's own discretisation error.
  const T drift = totalVolume<D>(fill, gridDelta) - before;
  check("deposit then etch loses nothing at the clamps",
        std::abs(drift) < 0.001 * moved,
        "drift " + std::to_string(drift) + " over " + std::to_string(moved) +
            " moved (" + std::to_string(100 * drift / moved) + "%)");
}

/// Support must be a GEOMETRIC question, not a material one. A sub-cell-thick
/// film of one material resting on a full substrate of another is attached
/// matter: no cell of the film ever reaches fill 1, so the film cannot seed
/// the support flood itself, and if the flood refuses to cross the material
/// boundary the whole film is classified unsupported and poured away. That is
/// a silent deletion of legal geometry, and it is exactly the case a selective
/// or conformal deposition produces on its first steps.
void checkThinFilmOnForeignSubstrate() {
  constexpr int D = 2;
  const T gridDelta = 1.0;
  auto cellSet = makeBlock<D>(gridDelta, -10., 10.);
  cs::LatticeMap<T, D> lattice(cellSet);
  const auto &dims = lattice.dims();

  const int substrateId = 1, filmId = 2, gasId = 0;
  std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
  std::vector<int> material(cellSet.getNumberOfCells(), gasId);

  // full substrate below y0, one row of a 0.4-filled film above it
  const int y0 = dims[1] / 2;
  for (int i = 0; i < dims[0]; ++i) {
    for (int j = 0; j < y0; ++j) {
      const int id = lattice.cellId({i, j});
      if (id < 0)
        continue;
      fill[id] = T(1);
      material[id] = substrateId;
    }
    const int id = lattice.cellId({i, y0});
    if (id < 0)
      continue;
    fill[id] = T(0.4);
    material[id] = filmId;
  }

  T filmBefore = 0;
  for (size_t c = 0; c < fill.size(); ++c)
    if (material[c] == filmId)
      filmBefore += fill[c];

  cs::VoxelAdvance<T, D> advance(lattice);
  std::vector<T> velocity(fill.size(), T(0)); // no motion: only the repairs run
  const auto step = advance.apply(fill, velocity, T(0), &material, gasId);

  T filmAfter = 0;
  for (size_t c = 0; c < fill.size(); ++c)
    if (material[c] == filmId)
      filmAfter += fill[c];

  check("a thin film on a foreign substrate survives",
        std::abs(filmAfter - filmBefore) < T(1e-9),
        "film volume " + std::to_string(filmBefore) + " -> " +
            std::to_string(filmAfter) + ", lost " +
            std::to_string(step.volumeLost));
}

/// The other half of the support rule: matter with nothing under it goes,
/// including on the lattice edge, where the field-continuation convention
/// used for gradients would otherwise let a speck anchor itself forever.
void checkDetachedSpecksRemoved() {
  constexpr int D = 2;
  const T gridDelta = 1.0;
  auto cellSet = makeBlock<D>(gridDelta, -10., 10.);
  cs::LatticeMap<T, D> lattice(cellSet);
  const auto &dims = lattice.dims();

  std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
  const int y0 = dims[1] / 2;
  for (int i = 0; i < dims[0]; ++i)
    for (int j = 0; j < y0; ++j) {
      const int id = lattice.cellId({i, j});
      if (id >= 0)
        fill[id] = T(1);
    }
  // three specks in the gas: interior, on the lateral edge, and a flat pair
  const std::array<std::array<int, 2>, 4> specks{
      {{dims[0] / 2, y0 + 3}, {0, y0 + 3}, {2, y0 + 5}, {3, y0 + 5}}};
  T strayBefore = 0;
  for (const auto &sp : specks) {
    const int id = lattice.cellId(sp);
    if (id >= 0) {
      fill[id] = T(0.3);
      strayBefore += T(0.3);
    }
  }

  cs::VoxelAdvance<T, D> advance(lattice);
  std::vector<T> velocity(fill.size(), T(0));
  const auto step = advance.apply(fill, velocity, T(0));

  T strayAfter = 0;
  for (const auto &sp : specks) {
    const int id = lattice.cellId(sp);
    if (id >= 0)
      strayAfter += fill[id];
  }
  check("detached specks are removed, edge and flat pair included",
        strayAfter < T(1e-9),
        "stray " + std::to_string(strayBefore) + " -> " +
            std::to_string(strayAfter) + ", lost " +
            std::to_string(step.volumeLost));
}

int main() {
  std::cout << "Voxel advance: moving a surface through filling fractions\n\n";

  std::cout << "1) mass\n";
  checkConservation();
  checkReversibility();

  std::cout << "\n2) a flat surface\n";
  checkFlatAdvance();

  std::cout << "\n3) a tilted surface: distance moved after 3.0 of travel\n";
  checkTiltedAdvance();

  std::cout << "\n4) support is geometric, not material\n";
  checkThinFilmOnForeignSubstrate();
  checkDetachedSpecksRemoved();

  std::cout << "\n";
  if (failures) {
    std::cout << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
