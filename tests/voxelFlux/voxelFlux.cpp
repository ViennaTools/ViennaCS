#include <csVoxelFlux.hpp>

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>

// Flux onto a voxel geometry, and whether it is the flux that was sent.
//
// The first check is the one the whole comparison hangs on, and it is not
// about geometry at all: an unobstructed flat surface must measure back the
// source flux. Rays carry a rate, cells divide by the area they show, and if
// that round trip is off by a factor then every profile computed from it is
// wrong by that factor -- silently, and in a way that looks like a difference
// between the two methods rather than a bug in one of them.
//
// Only once that holds does shadowing mean anything.

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

template <int D>
cs::DenseCellSet<T, D> makeDomain(T gridDelta, T lo, T hi) {
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

/// The flux density of the surface, weighted by the area each cell carries.
///
/// NOT a mean over cells that received something. A smeared interface spreads
/// its area over two or three cells, and the outermost of them may collect
/// nothing at all, so an unweighted mean over the cells with flux reports a
/// density that no part of the surface has. What is physically meaningful is
/// the total rate delivered over the total area presented, which is what a
/// rate law means by a flux.
template <int D>
T meanSurfaceFlux(const cs::LatticeMap<T, D> &lattice,
                  const std::vector<T> &flux, int margin,
                  const std::vector<T> &fill) {
  const cs::VoxelAdvance<T, D> areas(lattice);
  const auto &dims = lattice.dims();
  T collected = 0, area = 0;
  size_t sites = 1;
  for (int d = 0; d < D; ++d)
    sites *= static_cast<size_t>(dims[d]);
  std::array<int, D> idx{};
  for (size_t flat = 0; flat < sites; ++flat) {
    size_t rem = flat;
    bool edge = false;
    for (int d = 0; d < D; ++d) {
      idx[d] = static_cast<int>(rem % static_cast<size_t>(dims[d]));
      rem /= static_cast<size_t>(dims[d]);
      if (d < D - 1 && (idx[d] < margin || idx[d] >= dims[d] - margin))
        edge = true;
    }
    if (edge)
      continue;
    const int id = lattice.cellId(idx);
    if (id < 0)
      continue;
    const T a = areas.interfaceArea(fill, idx);
    if (a <= 0)
      continue;
    collected += flux[id] * a;
    area += a;
  }
  return area > 0 ? collected / area : 0;
}

/// A flat surface under an unobstructed source must measure the source flux.
template <int D> void checkBlanketNormalisation(const std::string &label) {
  const T gridDelta = 1.0;
  auto cellSet = makeDomain<D>(gridDelta, -12., 12.);
  cs::LatticeMap<T, D> lattice(cellSet);
  std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
  cs::fillFromSignedDistance<T, D>(lattice, fill,
                                   [](const viennacore::Vec3D<T> &p) {
                                     return p[D - 1] - T(0);
                                   });

  cs::VoxelFlux<T, D> flux(lattice, fill);
  const T sourceFlux = 100.0;
  const auto result = flux.trace(D == 2 ? 400000 : 800000, sourceFlux, 1.0);

  // The same, at a sticking small enough that a ray keeps its weight. What a
  // surface receives is what ARRIVES at it: the rate law applies the sticking
  // itself, so depositing the absorbed part here would apply it twice. Silane
  // sticks at 4e-4 at 900 K, where that mistake is a rate three orders of
  // magnitude too small -- and a check run only at sticking one cannot see it.
  for (const T sticking : {T(0.3), T(4.4e-4)}) {
    const auto low = flux.trace(D == 2 ? 200000 : 400000, sourceFlux, sticking);
    const T lowMean = meanSurfaceFlux<D>(lattice, low.flux, 4, fill);
    check(label + ": the flux is the incident one, not the absorbed one",
          std::abs(lowMean / sourceFlux - 1.0) < 0.05,
          "sticking " + std::to_string(sticking) + " gives " +
              std::to_string(lowMean / sourceFlux) + " of the source");
  }

  const T mean = meanSurfaceFlux<D>(lattice, result.flux, 4, fill);
  const T ratio = mean / sourceFlux;

  std::cout << "     " << label << " mean surface flux " << std::fixed
            << std::setprecision(2) << mean << " of " << sourceFlux
            << " sent  (ratio " << std::setprecision(3) << ratio << ")\n";
  check(label + ": a blanket measures the flux that was sent",
        std::abs(ratio - 1.0) < 0.05,
        "ratio " + std::to_string(ratio));
  // A handful still leave: a ray can reflect off a lateral boundary many
  // times over and run out of crossings, and one re-emitted upwards is gone
  // for good. It must be a handful, not a fraction.
  const T escaped = T(result.raysTraced - result.raysAbsorbed) /
                    T(result.raysTraced);
  check(label + ": essentially every ray finds the surface", escaped < 1e-3,
        std::to_string(100 * escaped) + "% escaped");
}

/// A trench, and what the sticking does to the flux that reaches its floor.
void checkShadowing() {
  constexpr int D = 2;
  const T gridDelta = 1.0;
  auto cellSet = makeDomain<D>(gridDelta, -30., 30.);
  cs::LatticeMap<T, D> lattice(cellSet);
  std::vector<T> fill(cellSet.getNumberOfCells(), T(0));

  // A trench 10 wide and 30 deep, in a substrate whose surface is at y = 0.
  cs::fillFromSignedDistance<T, D>(
      lattice, fill, [&](const viennacore::Vec3D<T> &p) {
        const T toSurface = p[1];
        const T inTrench =
            std::max(std::abs(p[0]) - 5.0, -(p[1] + 30.0));
        return std::max(toSurface, -inTrench);
      });

  std::cout << "     sticking     field flux      floor flux    floor/field\n";
  T ratioHigh = 0, ratioLow = 0;
  for (const T sticking : {1.0, 0.1}) {
    cs::VoxelFlux<T, D> flux(lattice, fill);
    const auto result = flux.trace(600000, 100.0, sticking);

    // The field is the flat top away from the trench; the floor is the bottom
    // of the trench.
    auto meanIn = [&](int i0, int i1, int j0, int j1) {
      T sum = 0;
      int n = 0;
      for (int i = i0; i < i1; ++i)
        for (int j = j0; j < j1; ++j) {
          const int id = lattice.cellId({i, j});
          if (id >= 0 && result.flux[id] > 0) {
            sum += result.flux[id];
            ++n;
          }
        }
      return n ? sum / n : 0.0;
    };
    const auto &dims = lattice.dims();
    const T field = meanIn(4, 18, 0, dims[1]);
    const T floor = meanIn(dims[0] / 2 - 3, dims[0] / 2 + 3, 0, 6);
    const T ratio = field > 0 ? floor / field : 0;
    if (sticking > 0.5)
      ratioHigh = ratio;
    else
      ratioLow = ratio;
    std::cout << "     " << std::fixed << std::setprecision(2) << std::setw(6)
              << sticking << std::setw(15) << field << std::setw(16) << floor
              << std::setprecision(4) << std::setw(15) << ratio << "\n";
  }

  check("a deep trench shadows its floor", ratioHigh < 0.5,
        "floor is " + std::to_string(ratioHigh) + " of the field at s = 1");
  check("a low sticking re-emits flux down the trench", ratioLow > ratioHigh,
        "s = 0.1 gives " + std::to_string(ratioLow) + " against " +
            std::to_string(ratioHigh) + " at s = 1");
}

int main() {
  std::cout << "Voxel flux: what the surface receives\n\n";

  std::cout << "1) normalisation on a blanket\n";
  checkBlanketNormalisation<2>("2D");
  checkBlanketNormalisation<3>("3D");

  std::cout << "\n2) shadowing in a trench\n";
  checkShadowing();

  std::cout << "\n";
  if (failures) {
    std::cout << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
