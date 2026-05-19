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
  T gridDelta = 0.05;

  auto sphere = ls::Domain<T, D>::New(bounds, boundaryConds, gridDelta);
  T center[D] = {0., 0.};
  T radius = 0.5;
  ls::MakeGeometry<T, D>(sphere, ls::Sphere<T, D>::New(center, radius)).apply();

  auto mesh = ls::Mesh<T>::New();
  ls::ToMesh<T, D>(sphere, mesh).apply();
  ls::VTKWriter<T>(mesh, "sphere_2.vtp").apply();

  ls::Expand<T, D>(sphere, 3).apply();
  ls::ToMesh<T, D>(sphere, mesh).apply();
  ls::VTKWriter<T>(mesh, "sphere_3.vtp").apply();

  auto levelSets = std::vector<ls::SmartPointer<ls::Domain<T, D>>>{};
  levelSets.push_back(sphere);

  cs::DenseCellSet<T, D> cellSet;
  int coverMaterial = 0;
  bool isAboveSurface = true;
  T depth = 0.6;
  cellSet.setCellSetPosition(isAboveSurface);
  cellSet.setCoverMaterial(coverMaterial);
  cellSet.fromLevelSets(levelSets, nullptr, depth);

  //   VC_TEST_ASSERT(cellSet.getDepth() == depth);
  //   VC_TEST_ASSERT(cellSet.getNumberOfCells() == 160);

  cellSet.writeVTU("sphere.vtu");
}
