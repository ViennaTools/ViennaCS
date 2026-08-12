#include <csDenseCellSet.hpp>

#include <lsMakeGeometry.hpp>
#include <lsMaterialMap.hpp>

#include <vcTestAsserts.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace ls = viennals;
namespace cs = viennacs;

using T = double;
constexpr int D = 2;

using LevelSet = ls::SmartPointer<ls::Domain<T, D>>;
using MaterialMap = ls::SmartPointer<ls::MaterialMap>;

namespace {

LevelSet makePlane(const T *bounds, ls::BoundaryConditionEnum *boundaryConds,
                   T gridDelta, const std::array<T, D> &origin,
                   const std::array<T, D> &normal) {
  auto plane = LevelSet::New(bounds, boundaryConds, gridDelta);
  ls::MakeGeometry<T, D>(plane, ls::Plane<T, D>::New(origin, normal)).apply();
  return plane;
}

cs::DenseCellSet<T, D> makeCellSet(const std::vector<LevelSet> &levelSets,
                                   T depth, MaterialMap materialMap = nullptr) {
  cs::DenseCellSet<T, D> cellSet;
  cellSet.setCellSetPosition(true);
  cellSet.setCoverMaterial(0);
  cellSet.enableEmbeddedBoundaries();
  cellSet.fromLevelSets(levelSets, materialMap, depth);
  return cellSet;
}

void testDisabledBoundaryPoints() {
  ls::BoundaryConditionEnum boundaryConds[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T bounds[2 * D] = {-1., 1., -1., 1.};
  constexpr T gridDelta = 0.25;

  auto plane = makePlane(bounds, boundaryConds, gridDelta, {0., 0.}, {0., 1.});
  cs::DenseCellSet<T, D> cellSet;
  cellSet.setCellSetPosition(true);
  cellSet.setCoverMaterial(0);
  cellSet.fromLevelSets({plane}, nullptr, 1.);

  VC_TEST_ASSERT(!cellSet.embeddedBoundariesEnabled())
  VC_TEST_ASSERT(!cellSet.hasEmbeddedBoundaries())
  VC_TEST_ASSERT(cellSet.numEmbeddedBoundaryPoints() == 0)
  for (unsigned cellId = 0; cellId < cellSet.getNumberOfCells(); ++cellId)
    VC_TEST_ASSERT(cellSet.getEmbeddedBoundaryPointIds(cellId).empty())
}

void testFlatBoundaryPoints() {
  ls::BoundaryConditionEnum boundaryConds[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T bounds[2 * D] = {-1., 1., -1., 1.};
  constexpr T gridDelta = 0.25;

  auto lower = makePlane(bounds, boundaryConds, gridDelta, {0., 0.}, {0., 1.});
  auto upper = makePlane(bounds, boundaryConds, gridDelta, {0., 0.5}, {0., 1.});
  auto cellSet = makeCellSet({lower, upper}, 1.);

  VC_TEST_ASSERT(cellSet.hasEmbeddedBoundaries())

  bool foundLower = false;
  bool foundUpper = false;
  for (const auto &point : cellSet.getEmbeddedBoundaryPoints()) {
    if (point.levelSetIndex == 0 && std::abs(point.coordinate[1]) < 1.e-12) {
      foundLower = true;
      VC_TEST_ASSERT(std::abs(std::abs(point.normal[1]) - 1.) < 1.e-12)
    }
    if (point.levelSetIndex == 1 &&
        std::abs(point.coordinate[1] - 0.5) < 1.e-12) {
      foundUpper = true;
      VC_TEST_ASSERT(std::abs(std::abs(point.normal[1]) - 1.) < 1.e-12)
    }
  }

  VC_TEST_ASSERT(foundLower)
  VC_TEST_ASSERT(foundUpper)

  bool foundCellAttachment = false;
  for (unsigned cellId = 0; cellId < cellSet.getNumberOfCells(); ++cellId) {
    if (!cellSet.getEmbeddedBoundaryPointIds(cellId).empty()) {
      foundCellAttachment = true;
      break;
    }
  }
  VC_TEST_ASSERT(foundCellAttachment)
}

void testBoundaryPointAttachmentAndMaterials() {
  ls::BoundaryConditionEnum boundaryConds[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T bounds[2 * D] = {-1., 1., -1., 1.};
  constexpr T gridDelta = 0.25;

  auto plane = makePlane(bounds, boundaryConds, gridDelta, {0., 0.}, {0., 1.});
  auto materialMap = MaterialMap::New();
  materialMap->insertNextMaterial(7);
  materialMap->insertNextMaterial(11);
  auto cellSet = makeCellSet({plane}, 1., materialMap);

  std::vector<bool> seen(cellSet.numEmbeddedBoundaryPoints(), false);
  for (unsigned cellId = 0; cellId < cellSet.getNumberOfCells(); ++cellId) {
    for (const auto pointId : cellSet.getEmbeddedBoundaryPointIds(cellId)) {
      VC_TEST_ASSERT(pointId < cellSet.numEmbeddedBoundaryPoints())
      const auto &point = cellSet.getEmbeddedBoundaryPoints()[pointId];
      VC_TEST_ASSERT(point.adjacentCell == static_cast<int>(cellId))
      VC_TEST_ASSERT(point.axis >= 0 && point.axis < D)
      VC_TEST_ASSERT(point.edgeFraction >= 0. && point.edgeFraction <= 1.)
      VC_TEST_ASSERT(std::abs(point.signedDistance) < 1.e-14)
      seen[pointId] = true;

      if (point.levelSetIndex == 0) {
        VC_TEST_ASSERT(point.negativeMaterial == 7)
        VC_TEST_ASSERT(point.positiveMaterial == 11)
      }
    }
  }

  for (const auto wasSeen : seen)
    VC_TEST_ASSERT(wasSeen)
}

void testRebuildDoesNotAppendGeometryOrBoundaryPoints() {
  ls::BoundaryConditionEnum boundaryConds[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T bounds[2 * D] = {-1., 1., -1., 1.};
  constexpr T gridDelta = 0.25;

  auto lower = makePlane(bounds, boundaryConds, gridDelta, {0., 0.}, {0., 1.});
  auto upper = makePlane(bounds, boundaryConds, gridDelta, {0., 0.5}, {0., 1.});
  auto cellSet = makeCellSet({lower, upper}, 1.);

  const auto firstCells = cellSet.getNumberOfCells();
  const auto firstNodes = cellSet.getNodes().size();
  const auto firstPoints = cellSet.numEmbeddedBoundaryPoints();

  cellSet.buildNeighborhood();
  cellSet.fromLevelSets({lower, upper}, nullptr, 1.);

  VC_TEST_ASSERT(cellSet.getNumberOfCells() == firstCells)
  VC_TEST_ASSERT(cellSet.getNodes().size() == firstNodes)
  VC_TEST_ASSERT(cellSet.numEmbeddedBoundaryPoints() == firstPoints)
  VC_TEST_ASSERT(cellSet.getScalarData("Material")->size() == firstCells)
  VC_TEST_ASSERT(cellSet.getFillingFractions()->size() == firstCells)

  cellSet.buildNeighborhood();
}

void testFaceBoundaryCache() {
  ls::BoundaryConditionEnum boundaryConds[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T bounds[2 * D] = {-1., 1., -1., 1.};
  constexpr T gridDelta = 0.25;
  constexpr T interfaceY = 0.10;

  auto plane =
      makePlane(bounds, boundaryConds, gridDelta, {0., interfaceY}, {0., 1.});
  auto cellSet = makeCellSet({plane}, 1.);

  // Face index for the y- side: axis*2 + 0 = (D-1)*2
  constexpr unsigned yMinusFace = (D - 1) * 2;

  // Find a cell just above the interface that has a boundary point on its y-
  // face.
  bool foundCell = false;
  for (unsigned cellId = 0; cellId < cellSet.getNumberOfCells(); ++cellId) {
    const auto center = cellSet.getCellCenter(cellId);
    if (center[D - 1] <= interfaceY || center[D - 1] >= interfaceY + gridDelta)
      continue;
    const int pointId = cellSet.getFaceBoundaryPointId(cellId, yMinusFace);
    if (pointId < 0)
      continue;
    foundCell = true;
    const T expectedDist = center[D - 1] - interfaceY;
    VC_TEST_ASSERT(expectedDist > 0.)
    VC_TEST_ASSERT(expectedDist < 0.5 * gridDelta)
    VC_TEST_ASSERT_ISCLOSE(cellSet.getFaceBoundaryDistance(cellId, yMinusFace),
                           expectedDist, 1e-10)
    break;
  }
  VC_TEST_ASSERT(foundCell)

  // getMinFaceBoundaryDistance must be strictly between 0 and gridDelta/2.
  const T minDist = cellSet.getMinFaceBoundaryDistance();
  VC_TEST_ASSERT(minDist > 0.)
  VC_TEST_ASSERT(minDist < 0.5 * gridDelta)
}

void testTiltedBoundaryPoints() {
  ls::BoundaryConditionEnum boundaryConds[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T bounds[2 * D] = {-1., 1., -1., 1.};
  constexpr T gridDelta = 0.25;
  constexpr T planeOffset = 0.1;

  auto tilted =
      makePlane(bounds, boundaryConds, gridDelta, {0., planeOffset}, {1., 1.});
  auto cellSet = makeCellSet({tilted}, 1.);

  const T invSqrt2 = 1. / std::sqrt(2.);
  bool foundTilted = false;
  bool foundOnPlane = false;
  bool foundTiltedNormal = false;
  unsigned onPlaneCount = 0;
  for (const auto &point : cellSet.getEmbeddedBoundaryPoints()) {
    if (point.levelSetIndex != 0)
      continue;
    foundTilted = true;
    if (std::abs(point.coordinate[0] + point.coordinate[1] - planeOffset) >
        1.e-10)
      continue;
    foundOnPlane = true;
    ++onPlaneCount;
    const T normalNorm = std::sqrt(point.normal[0] * point.normal[0] +
                                   point.normal[1] * point.normal[1]);
    if (normalNorm < 1.e-12)
      continue;
    const T alignment =
        std::abs(point.normal[0] * invSqrt2 + point.normal[1] * invSqrt2);
    if (alignment > 0.99)
      foundTiltedNormal = true;
  }

  VC_TEST_ASSERT(foundTilted)
  VC_TEST_ASSERT(foundOnPlane)
  VC_TEST_ASSERT(foundTiltedNormal)
  VC_TEST_ASSERT(onPlaneCount == cellSet.getEmbeddedBoundaryPoints().size())
}

} // namespace

int main() {
  testDisabledBoundaryPoints();
  testFlatBoundaryPoints();
  testBoundaryPointAttachmentAndMaterials();
  testRebuildDoesNotAppendGeometryOrBoundaryPoints();
  testFaceBoundaryCache();
  testTiltedBoundaryPoints();
  return 0;
}
