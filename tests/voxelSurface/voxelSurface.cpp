#include <csVoxelSurface.hpp>

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>

#include <rayGeometryTriangle.hpp>

#include <cmath>
#include <iostream>
#include <random>

// The claim this file has to establish is the one everything downstream rests
// on: walking the grid and intersecting the surface find the SAME cell.
//
// They are two descriptions of one geometry. A binary voxel set's boundary is
// exactly its exposed faces, so a DDA that stops at the first solid cell and a
// ray tracer that stops at the first face must agree -- on which cell, and on
// how far away it is. Where they disagree, one of them is wrong, and it is far
// cheaper to find that out here than in a profile.
//
// Ground truth is a brute-force intersection against every face, which is slow
// and obviously right. embree is then checked against the same truth, so the
// ViennaRay plumbing the level-set arm uses is verified on this geometry too.

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
cs::DenseCellSet<T, D> makeTrenchCellSet(T gridDelta, T depth) {
  ls::BoundaryConditionEnum bc[D];
  for (int i = 0; i < D - 1; ++i)
    bc[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  bc[D - 1] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;

  T bounds[2 * D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -5.;
    bounds[2 * i + 1] = 5.;
  }
  T origin[D] = {}, normal[D] = {};
  normal[D - 1] = 1.;

  auto substrate =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bc, gridDelta);
  ls::MakeGeometry<T, D>(substrate,
                         ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
      .apply();

  T lo[D], hi[D];
  for (int i = 0; i < D - 1; ++i) {
    lo[i] = -2.;
    hi[i] = 2.;
  }
  lo[D - 1] = -3.;
  hi[D - 1] = 1.;
  auto trench = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bc, gridDelta);
  ls::MakeGeometry<T, D>(trench, ls::SmartPointer<ls::Box<T, D>>::New(lo, hi))
      .apply();
  ls::BooleanOperation<T, D>(substrate, trench,
                             ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  std::vector<ls::SmartPointer<ls::Domain<T, D>>> lss{substrate};
  cs::DenseCellSet<T, D> cellSet;
  cellSet.setCellSetPosition(false);
  cellSet.setCoverMaterial(0);
  cellSet.fromLevelSets(lss, nullptr, depth);
  return cellSet;
}

/// Nearest face hit by the ray, by trying every one of them.
template <int D>
int bruteForceHit(const cs::VoxelSurface<T, D> &surface,
                  const std::array<T, D> &origin, const std::array<T, D> &dir,
                  T &tHit) {
  int best = -1;
  tHit = std::numeric_limits<T>::max();
  for (size_t f = 0; f < surface.numFaces(); ++f) {
    const auto &n = surface.normal[f];
    T denom = 0;
    for (int d = 0; d < D; ++d)
      denom += n[d] * dir[d];
    if (std::abs(denom) < 1e-12)
      continue;
    const auto &p0 = surface.nodes[surface.faces[f][0]];
    T num = 0;
    for (int d = 0; d < D; ++d)
      num += n[d] * (p0[d] - origin[d]);
    const T t = num / denom;
    if (t < 1e-9 || t >= tHit)
      continue;

    // Inside the face? Both a segment and an axis-aligned quad are a box in
    // the plane, so a bounding-box test on the face's own nodes is exact.
    T pt[D];
    for (int d = 0; d < D; ++d)
      pt[d] = origin[d] + dir[d] * t;
    bool inside = true;
    for (int d = 0; d < D; ++d) {
      T lo = std::numeric_limits<T>::max(),
        hi = std::numeric_limits<T>::lowest();
      for (int k = 0; k < D; ++k) {
        lo = std::min(lo, surface.nodes[surface.faces[f][k]][d]);
        hi = std::max(hi, surface.nodes[surface.faces[f][k]][d]);
      }
      if (pt[d] < lo - 1e-9 || pt[d] > hi + 1e-9)
        inside = false;
    }
    if (!inside)
      continue;
    tHit = t;
    best = static_cast<int>(f);
  }
  return best;
}

template <int D> void checkNormals(const std::string &label) {
  const T gridDelta = 0.5;
  auto cellSet = makeTrenchCellSet<D>(gridDelta, -4.);
  cs::LatticeMap<T, D> lattice(cellSet);
  auto surface = cs::extractVoxelSurface<T, D>(lattice);

  check(label + ": the geometry has exposed faces", surface.numFaces() > 0,
        std::to_string(surface.numFaces()) + " faces");

  // A face normal must point away from the cell that owns it.
  int inward = 0;
  for (size_t f = 0; f < surface.numFaces(); ++f) {
    const auto centre = cellSet.getCellCenter(surface.cellId[f]);
    T dot = 0;
    for (int d = 0; d < D; ++d) {
      T faceMid = 0;
      for (int k = 0; k < D; ++k)
        faceMid += surface.nodes[surface.faces[f][k]][d];
      faceMid /= static_cast<T>(D);
      dot += surface.normal[f][d] * (faceMid - centre[d]);
    }
    if (dot <= 0)
      ++inward;
  }
  check(label + ": every face normal points out of its own cell", inward == 0,
        std::to_string(inward) + " pointing inward");

  // The stored normal is assigned, so the check above tests where the face was
  // placed, not how it was wound. ViennaRay derives its normals from the
  // winding -- cross(p1-p0, p2-p0) for a triangle, (-dy, dx) for a segment --
  // so a face wound backwards would trace with an inverted normal and reflect
  // rays into the solid. Recompute it the way ViennaRay does and compare.
  int misWound = 0;
  for (size_t f = 0; f < surface.numFaces(); ++f) {
    const auto &p0 = surface.nodes[surface.faces[f][0]];
    const auto &p1 = surface.nodes[surface.faces[f][1]];
    viennacore::Vec3D<T> wound{0, 0, 0};
    if constexpr (D == 2) {
      wound = {-(p1[1] - p0[1]), p1[0] - p0[0], 0};
    } else {
      const auto &p2 = surface.nodes[surface.faces[f][2]];
      const viennacore::Vec3D<T> a{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
      const viennacore::Vec3D<T> b{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
      wound = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
               a[0] * b[1] - a[1] * b[0]};
    }
    T len = std::sqrt(wound[0] * wound[0] + wound[1] * wound[1] +
                      wound[2] * wound[2]);
    T dot = 0;
    for (int d = 0; d < D; ++d)
      dot += (wound[d] / len) * surface.normal[f][d];
    if (dot < 0.999)
      ++misWound;
  }
  check(label + ": the winding gives the same normal ViennaRay would compute",
        misWound == 0, std::to_string(misWound) + " wound backwards");

  // A face is exposed by definition: solid on the inside, not solid on the
  // outside. Checking that directly catches buried faces, which no ray test
  // can see -- a ray from outside meets the outer face first either way -- but
  // which would inflate the geometry and give every surface cell more area
  // than it has.
  int buried = 0;
  for (size_t f = 0; f < surface.numFaces(); ++f) {
    const auto centre = cellSet.getCellCenter(surface.cellId[f]);
    viennacore::Vec3D<T> beyond{0, 0, 0};
    for (int d = 0; d < D; ++d)
      beyond[d] = centre[d] + surface.normal[f][d] * gridDelta;
    if (cellSet.getIndex(beyond) >= 0)
      ++buried;
  }
  check(label + ": no face has solid on both sides", buried == 0,
        std::to_string(buried) + " buried faces");
}

template <int D>
void checkAgainstTraversal(const std::string &label, int rays) {
  const T gridDelta = 0.5;
  auto cellSet = makeTrenchCellSet<D>(gridDelta, -4.);
  cs::LatticeMap<T, D> lattice(cellSet);
  cs::GridTraversal<T, D> traversal(lattice);
  auto surface = cs::extractVoxelSurface<T, D>(lattice);

  std::mt19937 rng(11);
  std::uniform_real_distribution<T> lateral(-4.9, 4.9);
  std::uniform_real_distribution<T> tilt(-0.6, 0.6);

  int compared = 0, cellMismatch = 0, distMismatch = 0;
  T worstDist = 0;
  for (int r = 0; r < rays; ++r) {
    std::array<T, D> origin{}, dir{};
    for (int d = 0; d < D - 1; ++d) {
      origin[d] = lateral(rng);
      dir[d] = tilt(rng);
    }
    origin[D - 1] = 6.0;
    dir[D - 1] = -1.0;
    T norm = 0;
    for (int d = 0; d < D; ++d)
      norm += dir[d] * dir[d];
    norm = std::sqrt(norm);
    for (int d = 0; d < D; ++d)
      dir[d] /= norm;

    // The walk: the first cell that is solid, and where the ray entered it.
    int walkCell = -1;
    T walkT = 0;
    traversal.traverse(origin, dir, [&](cs::GridStep<T, D> s) {
      const int id = lattice.cellId(s.index);
      if (id < 0)
        return true;
      walkCell = id;
      walkT = s.tEntry;
      return false; // an interaction: stop
    });

    // The surface: the nearest face, and the cell behind it.
    T faceT = 0;
    const int face = bruteForceHit<D>(surface, origin, dir, faceT);
    const int faceCell = face >= 0 ? surface.cellId[face] : -1;

    if (walkCell < 0 && faceCell < 0)
      continue; // the ray missed the geometry entirely
    ++compared;
    if (walkCell != faceCell)
      ++cellMismatch;
    if (std::abs(walkT - faceT) > 1e-7)
      ++distMismatch;
    worstDist = std::max(worstDist, std::abs(walkT - faceT));
  }

  check(label + ": the walk and the surface agree on which cell is hit",
        cellMismatch == 0,
        std::to_string(compared) + " rays, " + std::to_string(cellMismatch) +
            " disagreed");
  check(label + ": and on how far away it is", distMismatch == 0,
        "worst " + std::to_string(worstDist));
}

/// The same surface through embree, which is the tracer the level-set arm
/// uses. If this disagrees, the two arms are not seeing the same geometry.
void checkEmbree(int rays) {
  constexpr int D = 3;
  const T gridDelta = 0.5;
  auto cellSet = makeTrenchCellSet<D>(gridDelta, -4.);
  cs::LatticeMap<T, D> lattice(cellSet);
  cs::GridTraversal<T, D> traversal(lattice);
  auto surface = cs::extractVoxelSurface<T, D>(lattice);

  std::vector<viennacore::Vec3D<float>> points;
  points.reserve(surface.nodes.size());
  for (const auto &n : surface.nodes)
    points.push_back({(float)n[0], (float)n[1], (float)n[2]});
  std::vector<viennacore::Vec3D<unsigned>> tris;
  tris.reserve(surface.numFaces());
  for (const auto &f : surface.faces)
    tris.push_back({f[0], f[1], f[2]});

  auto device = rtcNewDevice("");
  viennaray::GeometryTriangle<float, 3> geometry;
  geometry.initGeometry(device, points, tris);

  // The geometry does not own a scene -- ViennaRay's kernel builds one and
  // attaches the geometry and the boundary to it. Here there is no boundary,
  // so a scene with this geometry alone is what a first-hit query needs.
  auto scene = rtcNewScene(device);
  rtcAttachGeometry(scene, geometry.getRTCGeometry());
  rtcCommitScene(scene);

  check("embree: the geometry keeps every face",
        geometry.getNumPrimitives() == surface.numFaces(),
        std::to_string(geometry.getNumPrimitives()) + " of " +
            std::to_string(surface.numFaces()));

  std::mt19937 rng(23);
  std::uniform_real_distribution<T> lateral(-4.9, 4.9);
  std::uniform_real_distribution<T> tilt(-0.6, 0.6);

  int compared = 0, mismatch = 0;
  T worstT = 0;
  for (int r = 0; r < rays; ++r) {
    std::array<T, D> origin{}, dir{};
    for (int d = 0; d < D - 1; ++d) {
      origin[d] = lateral(rng);
      dir[d] = tilt(rng);
    }
    origin[D - 1] = 6.0;
    dir[D - 1] = -1.0;
    T norm = 0;
    for (int d = 0; d < D; ++d)
      norm += dir[d] * dir[d];
    norm = std::sqrt(norm);
    for (int d = 0; d < D; ++d)
      dir[d] /= norm;

    int walkCell = -1;
    T walkT = 0;
    traversal.traverse(origin, dir, [&](cs::GridStep<T, D> s) {
      const int id = lattice.cellId(s.index);
      if (id < 0)
        return true;
      walkCell = id;
      walkT = s.tEntry;
      return false;
    });

    RTCRayHit rayHit;
    rayHit.ray.org_x = (float)origin[0];
    rayHit.ray.org_y = (float)origin[1];
    rayHit.ray.org_z = (float)origin[2];
    rayHit.ray.dir_x = (float)dir[0];
    rayHit.ray.dir_y = (float)dir[1];
    rayHit.ray.dir_z = (float)dir[2];
    rayHit.ray.tnear = 1e-4f;
    rayHit.ray.tfar = std::numeric_limits<float>::max();
    rayHit.ray.mask = -1;
    rayHit.ray.flags = 0;
    rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayHit.hit.primID = RTC_INVALID_GEOMETRY_ID;
    rtcIntersect1(scene, &rayHit);

    const int embreeCell = rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID
                               ? -1
                               : surface.cellId[rayHit.hit.primID];
    if (walkCell < 0 && embreeCell < 0)
      continue;
    ++compared;
    if (walkCell != embreeCell)
      ++mismatch;
    else
      worstT = std::max(worstT, std::abs(walkT - (T)rayHit.ray.tfar));
  }

  check("embree: the traced surface and the grid walk hit the same cell",
        mismatch == 0,
        std::to_string(compared) + " rays, " + std::to_string(mismatch) +
            " disagreed");
  check("embree: and at the same distance, to float precision", worstT < 1e-4,
        "worst " + std::to_string(worstT));

  rtcReleaseScene(scene);
  geometry.releaseGeometry();
  rtcReleaseDevice(device);
}

int main() {
  std::cout << "Voxel surface: the staircase a ray tracer sees\n\n";

  std::cout << "1) the extracted faces\n";
  checkNormals<2>("2D");
  checkNormals<3>("3D");

  std::cout << "\n2) the grid walk against the surface it implies\n";
  checkAgainstTraversal<2>("2D", 2000);
  checkAgainstTraversal<3>("3D", 800);

  std::cout << "\n3) the same surface through embree\n";
  checkEmbree(800);

  std::cout << "\n";
  if (failures) {
    std::cout << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
