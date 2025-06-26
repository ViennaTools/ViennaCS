#include <csDenseCellSet.hpp>

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>

#include <vcTestAsserts.hpp>

namespace ls = viennals;
namespace cs = viennacs;

int main() {
  using T = double;
  constexpr int D = 2;

  // two plane geometries
  ls::BoundaryConditionEnum boundaryConds[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  T bounds[2 * D] = {-1., 1., -1., 1.};
  T gridDelta = 0.02;

  T origin[D] = {0., 0.};
  T normal[D] = {0., 0.};
  normal[D - 1] = 1.;

  auto plane1 = ls::Domain<T, D>::New(bounds, boundaryConds, gridDelta);
  ls::MakeGeometry<T, D>(plane1, ls::Plane<T, D>::New(origin, normal)).apply();

  auto box = ls::Domain<T, D>::New(bounds, boundaryConds, gridDelta);
  T minPoint[D] = {-0.5, -0.5};
  T maxPoint[D] = {0.5, 0.0};
  ls::MakeGeometry<T, D>(box, ls::Box<T, D>::New(minPoint, maxPoint)).apply();

  ls::BooleanOperation<T, D>(plane1, box,
                             ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  auto levelSets = std::vector<ls::SmartPointer<ls::Domain<T, D>>>{};
  levelSets.push_back(plane1);

  cs::DenseCellSet<T, D> cellSet;
  int coverMaterial = 0;
  bool isAboveSurface = true;
  T depth = 1.;
  cellSet.setCellSetPosition(isAboveSurface);
  cellSet.setCoverMaterial(coverMaterial);
  cellSet.fromLevelSets(levelSets, nullptr, depth);

  //   VC_TEST_ASSERT(cellSet.getDepth() == depth);
  //   VC_TEST_ASSERT(cellSet.getNumberOfCells() == 160);

  cellSet.writeVTU("trench.vtu");
}
