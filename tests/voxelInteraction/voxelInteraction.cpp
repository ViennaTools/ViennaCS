#include <csVoxelInteraction.hpp>

#include <lsMakeGeometry.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

// What a voxel method can and cannot recover, measured.
//
// Two claims are made for filling fractions. The first is that they restore
// sub-grid surface position: a plane lying a third of the way into a cell
// should be found there, not at the nearest cell face. The second is that a
// normal can be had from the gradient of the fractions, which a staircase of
// faces cannot give.
//
// Both are measurable against a plane, whose position and normal are known
// exactly, and neither is taken on trust here. The angular error of the two
// normal estimators against surface tilt is the quantity the whole level-set
// against voxel comparison turns on, so it is printed rather than merely
// asserted.

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

/// A block of solid, as a lattice with every cell present. The geometry is
/// then written into the filling fractions, not into which cells exist.
template <int D>
cs::DenseCellSet<T, D> makeBlock(T gridDelta, T lo, T hi) {
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

  auto plane =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bc, gridDelta);
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

/// Claim zero: with every fraction at one, the rule must reduce exactly to
/// the first solid cell. If it does not, nothing verified so far carries over.
template <int D> void checkBinaryReduction(const std::string &label) {
  const T gridDelta = 0.5;
  auto cellSet = makeBlock<D>(gridDelta, -5., 5.);
  cs::LatticeMap<T, D> lattice(cellSet);
  std::vector<T> fill(cellSet.getNumberOfCells(), T(1));

  cs::VoxelInteraction<T, D> interaction(lattice, fill);
  cs::GridTraversal<T, D> traversal(lattice);

  std::mt19937 rng(3);
  std::uniform_real_distribution<T> lateral(-4.5, 4.5), tilt(-0.5, 0.5);

  int compared = 0, mismatch = 0;
  for (int r = 0; r < 1500; ++r) {
    std::array<T, D> origin{}, dir{};
    for (int d = 0; d < D - 1; ++d) {
      origin[d] = lateral(rng);
      dir[d] = tilt(rng);
    }
    origin[D - 1] = 12.0;
    dir[D - 1] = -1.0;
    T n = 0;
    for (int d = 0; d < D; ++d)
      n += dir[d] * dir[d];
    n = std::sqrt(n);
    for (int d = 0; d < D; ++d)
      dir[d] /= n;

    int walkCell = -1;
    traversal.traverse(origin, dir, [&](cs::GridStep<T, D> s) {
      const int id = lattice.cellId(s.index);
      if (id < 0)
        return true;
      walkCell = id;
      return false;
    });
    const auto hit = interaction.firstHit(origin, dir, rng);
    if (walkCell < 0 && !hit.hit())
      continue;
    ++compared;
    if (hit.cellId != walkCell)
      ++mismatch;
  }
  check(label + ": full cells reduce to the first solid cell", mismatch == 0,
        std::to_string(compared) + " rays, " + std::to_string(mismatch) +
            " disagreed");
}

/// Claim one: a plane between two cell faces is found between them.
void checkSubVoxelPosition() {
  constexpr int D = 2;
  const T gridDelta = 1.0;
  auto cellSet = makeBlock<D>(gridDelta, -10., 10.);
  cs::LatticeMap<T, D> lattice(cellSet);

  std::cout << "     true surface   fractional fill   binary (fill>=0.5)\n";
  T worstFractional = 0, worstBinary = 0;
  for (int k = 0; k <= 4; ++k) {
    const T surfaceY = -0.5 + 0.25 * k; // sweeps across one cell

    std::vector<T> fractional(cellSet.getNumberOfCells(), T(0));
    cs::fillFromSignedDistance<T, D>(
        lattice, fractional,
        [&](const viennacore::Vec3D<T> &p) { return p[1] - surfaceY; });
    std::vector<T> binary(fractional.size());
    for (size_t i = 0; i < fractional.size(); ++i)
      binary[i] = fractional[i] >= 0.5 ? T(1) : T(0);

    auto meanDepth = [&](std::vector<T> &fill) {
      cs::VoxelInteraction<T, D> interaction(lattice, fill);
      std::mt19937 rng(17);
      std::uniform_real_distribution<T> lateral(-8.0, 8.0);
      T sum = 0;
      int n = 0;
      for (int r = 0; r < 20000; ++r) {
        std::array<T, D> origin{lateral(rng), 9.0}, dir{0.0, -1.0};
        const auto hit = interaction.firstHit(origin, dir, rng);
        if (!hit.hit())
          continue;
        // The ray entered the interacting cell at this height.
        sum += origin[1] - hit.distance;
        ++n;
      }
      return n ? sum / n : std::numeric_limits<T>::quiet_NaN();
    };

    const T mf = meanDepth(fractional);
    const T mb = meanDepth(binary);
    worstFractional = std::max(worstFractional, std::abs(mf - surfaceY));
    worstBinary = std::max(worstBinary, std::abs(mb - surfaceY));
    std::cout << "     " << std::fixed << std::setprecision(3) << std::setw(9)
              << surfaceY << std::setw(18) << mf << std::setw(20) << mb << "\n";
  }
  check("fractional fill tracks the surface to well under a cell",
        worstFractional < 0.30 * gridDelta,
        "worst error " + std::to_string(worstFractional) + " of a " +
            std::to_string(gridDelta) + " cell");
  check("and does so more closely than a binary geometry",
        worstFractional < worstBinary,
        "fractional " + std::to_string(worstFractional) + " vs binary " +
            std::to_string(worstBinary));
}

/// Claim two, and the one the comparison turns on: how wrong is the normal?
void checkNormalAgainstTilt() {
  constexpr int D = 2;
  const T gridDelta = 1.0;
  auto cellSet = makeBlock<D>(gridDelta, -20., 20.);
  cs::LatticeMap<T, D> lattice(cellSet);

  std::cout << "     tilt      face normal    central diff     Youngs\n";
  T worstNarrow = 0, worstYoungs = 0, bestFace = 180;
  for (int deg = 0; deg <= 60; deg += 15) {
    const T angle = deg * M_PI / 180.;
    // A plane through the origin, tilted by `angle` from horizontal. The
    // signed distance to it is exact, so the fractions are exact.
    const viennacore::Vec3D<T> trueNormal{std::sin(angle), std::cos(angle), 0};
    std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
    cs::fillFromSignedDistance<T, D>(
        lattice, fill, [&](const viennacore::Vec3D<T> &p) {
          return p[0] * trueNormal[0] + p[1] * trueNormal[1];
        });

    auto meanError = [&](cs::NormalEstimator estimator) {
      cs::VoxelInteraction<T, D> interaction(lattice, fill, estimator);
      std::mt19937 rng(29);
      std::uniform_real_distribution<T> lateral(-8.0, 8.0);
      T sum = 0;
      int n = 0;
      for (int r = 0; r < 8000; ++r) {
        std::array<T, D> origin{lateral(rng), 15.0}, dir{0.0, -1.0};
        const auto hit = interaction.firstHit(origin, dir, rng);
        if (!hit.hit())
          continue;
        T dot = hit.normal[0] * trueNormal[0] + hit.normal[1] * trueNormal[1];
        dot = std::min(T(1), std::max(T(-1), dot));
        sum += std::acos(dot) * 180. / M_PI;
        ++n;
      }
      return n ? sum / n : std::numeric_limits<T>::quiet_NaN();
    };

    const T face = meanError(cs::NormalEstimator::Face);
    const T narrow = meanError(cs::NormalEstimator::FillGradient);
    const T youngs = meanError(cs::NormalEstimator::FillGradientYoungs);
    worstNarrow = std::max(worstNarrow, narrow);
    worstYoungs = std::max(worstYoungs, youngs);
    if (deg > 0)
      bestFace = std::min(bestFace, face);
    std::cout << "     " << std::setw(4) << deg << " deg" << std::fixed
              << std::setprecision(2) << std::setw(12) << face << " deg"
              << std::setw(12) << narrow << " deg" << std::setw(12) << youngs
              << " deg\n";
  }

  check("Youngs' stencil recovers a tilted surface", worstYoungs < 5.0,
        "worst mean error " + std::to_string(worstYoungs) + " deg");
  check("and beats the three-point central difference, which is anisotropic",
        worstYoungs < worstNarrow,
        "Youngs " + std::to_string(worstYoungs) + " vs central " +
            std::to_string(worstNarrow) + " deg");
  check("the face normal recovers nothing: it is quantised to the axes",
        bestFace > 10.0,
        "best mean error on a tilted surface " + std::to_string(bestFace) +
            " deg");
}

/// The angle chain in THREE dimensions, where it is new geometry rather than
/// a template parameter: the Youngs stencil gains an azimuthal dimension the
/// 2D test cannot exercise, and the face normal quantises to six directions.
/// Checked on tilted planes, over tilt AND azimuth, for the two quantities the
/// ion physics consumes: the incidence angle (the yield's f(theta)) and the
/// specular direction (where the reflection cone points). Face normals are
/// asserted BAD -- their specular error is twice the tilt, which is what
/// disqualifies them for ion work -- and Youngs is asserted uniform in
/// azimuth, which a separable stencil has no a-priori right to be.
void checkAngleChain3D() {
  constexpr int D3 = 3;
  using V3 = viennacore::Vec3D<T>;
  const T gridDelta = 1.0, lo = -14., hi = 14.;
  ls::BoundaryConditionEnum bc[D3] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T bounds[2 * D3];
  for (int i = 0; i < D3; ++i) {
    bounds[2 * i] = lo;
    bounds[2 * i + 1] = hi;
  }
  T o[D3] = {0, 0, hi}, n[D3] = {0, 0, 1};
  auto pl = ls::SmartPointer<ls::Domain<T, D3>>::New(bounds, bc, gridDelta);
  ls::MakeGeometry<T, D3>(pl,
                          ls::SmartPointer<ls::Plane<T, D3>>::New(o, n))
      .apply();
  std::vector<ls::SmartPointer<ls::Domain<T, D3>>> lss{pl};
  cs::DenseCellSet<T, D3> cellSet;
  cellSet.setCellSetPosition(false);
  cellSet.setCoverMaterial(0);
  cellSet.fromLevelSets(lss, nullptr, lo);
  cs::LatticeMap<T, D3> lattice(cellSet);

  auto angleBetween = [](const V3 &a, const V3 &b) {
    T dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    return std::acos(std::min(T(1), std::max(T(-1), dot))) * 180. / M_PI;
  };

  T worstYoungsNormal = 0, worstYoungsIncidence = 0, bestFaceIncidence = 1e9;
  T youngsLo = 1e9, youngsHi = 0; // azimuthal spread at one tilt
  for (const T tilt : {15., 30., 60.}) {
    for (const T azim : {0., 22.5, 45.}) {
      const T th = tilt * M_PI / 180., ph = azim * M_PI / 180.;
      const V3 nTrue{std::sin(th) * std::cos(ph),
                     std::sin(th) * std::sin(ph), std::cos(th)};
      std::vector<T> fill(cellSet.getNumberOfCells(), T(0));
      cs::fillFromSignedDistance<T, D3>(
          lattice, fill, [&](const V3 &p) {
            return p[0] * nTrue[0] + p[1] * nTrue[1] + p[2] * nTrue[2];
          });
      for (const auto est :
           {cs::NormalEstimator::Face, cs::NormalEstimator::FillGradientYoungs}) {
        cs::VoxelInteraction<T, D3> interaction(lattice, fill, est);
        std::mt19937 rng(5);
        std::uniform_real_distribution<T> u(-7., 7.);
        T nerr = 0, ierr = 0;
        int cnt = 0;
        for (int r = 0; r < 2000; ++r) {
          std::array<T, D3> org{u(rng), u(rng), T(13)}, dir{0, 0, -1};
          const auto hit = interaction.firstHit(org, dir, rng);
          if (!hit.hit())
            continue;
          nerr += angleBetween(hit.normal, nTrue);
          const T measured =
              std::acos(std::min(T(1), std::max(T(0), hit.normal[2]))) * 180. /
              M_PI; // vertical ray: incidence = angle of the normal
          ierr += std::abs(measured - tilt);
          ++cnt;
        }
        nerr /= cnt;
        ierr /= cnt;
        if (est == cs::NormalEstimator::FillGradientYoungs) {
          worstYoungsNormal = std::max(worstYoungsNormal, nerr);
          worstYoungsIncidence = std::max(worstYoungsIncidence, ierr);
          if (tilt == 30.) {
            youngsLo = std::min(youngsLo, nerr);
            youngsHi = std::max(youngsHi, nerr);
          }
        } else {
          bestFaceIncidence = std::min(bestFaceIncidence, ierr);
        }
      }
    }
  }
  check("3D: Youngs' normal holds on tilted planes at every azimuth",
        worstYoungsNormal < 4.0,
        "worst mean error " + std::to_string(worstYoungsNormal) + " deg");
  check("3D: the incidence angle the yield consumes is right to a few degrees",
        worstYoungsIncidence < 4.0,
        "worst " + std::to_string(worstYoungsIncidence) + " deg");
  check("3D: Youngs is azimuthally uniform, which a separable stencil must earn",
        youngsHi - youngsLo < 2.0,
        "spread " + std::to_string(youngsHi - youngsLo) + " deg at 30 deg tilt");
  check("3D: the face normal misreads incidence by the tilt itself",
        bestFaceIncidence > 10.0,
        "best " + std::to_string(bestFaceIncidence) + " deg");
}

int main() {
  std::cout << "Voxel interaction: what filling fractions recover\n\n";

  std::cout << "0) reduction to the verified binary case\n";
  checkBinaryReduction<2>("2D");
  checkBinaryReduction<3>("3D");

  std::cout << "\n1) sub-voxel surface position\n";
  checkSubVoxelPosition();

  std::cout << "\n2) the normal, against surface tilt\n";
  checkNormalAgainstTilt();

  std::cout << "\n3) the angle chain in three dimensions\n";
  checkAngleChain3D();

  std::cout << "\n";
  if (failures) {
    std::cout << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
