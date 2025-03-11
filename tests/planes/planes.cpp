#include <csDenseCellSet.hpp>

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
  T gridDelta = 0.2;

  T origin[D] = {};
  T normal[D] = {};
  normal[D - 1] = 1.;

  auto plane1 =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryConds, gridDelta);
  ls::MakeGeometry<T, D>(plane1,
                         ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
      .apply();

  origin[D - 1] = 1.;
  auto plane2 =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryConds, gridDelta);
  ls::MakeGeometry<T, D>(plane2,
                         ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
      .apply();

  auto levelSets = std::vector<ls::SmartPointer<ls::Domain<T, D>>>{};
  levelSets.push_back(plane1);
  levelSets.push_back(plane2);

  cs::DenseCellSet<T, D> cellSet;
  int coverMaterial = 0;
  bool isAboveSurface = true;
  T depth = 3.;
  cellSet.setCellSetPosition(isAboveSurface);
  cellSet.setCoverMaterial(coverMaterial);
  cellSet.fromLevelSets(levelSets, nullptr, depth);

  VC_TEST_ASSERT(cellSet.getDepth() == depth);
  VC_TEST_ASSERT(cellSet.getNumberOfCells() == 160);
}
