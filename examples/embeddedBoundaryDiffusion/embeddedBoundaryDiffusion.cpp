#include <csBoundaryDiffusionSolver.hpp>
#include <csDenseCellSet.hpp>
#include <csDiffusionSolver.hpp>

#include <lsMakeGeometry.hpp>
#include <lsMaterialMap.hpp>
#include <lsVTKWriter.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

namespace cs = viennacs;
namespace ls = viennals;

using T = double;
constexpr int D = 2;

constexpr int substrateMaterial = 0;
constexpr int coverMaterial = 1;

using CellSet = cs::SmartPointer<cs::DenseCellSet<T, D>>;
using LevelSet = cs::SmartPointer<ls::Domain<T, D>>;

namespace {

LevelSet makeTiltedSurface(T gridDelta, T halfWidth, T yMin, T yMax, T y0,
                           T slope) {
  T bounds[2 * D] = {-halfWidth, halfWidth, yMin, yMax};
  ls::BoundaryConditionEnum bc[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};

  auto surface = LevelSet::New(bounds, bc, gridDelta);
  const std::array<T, D> origin{0., y0};
  const std::array<T, D> normal{-slope, 1.};
  ls::MakeGeometry<T, D>(surface, ls::Plane<T, D>::New(origin, normal)).apply();
  return surface;
}

CellSet makeCellSet(T gridDelta, bool withEmbeddedBoundaries) {
  constexpr T halfWidth = 5.;
  constexpr T yMin = -3.;
  constexpr T yMax = 4.;
  constexpr T surfaceY = 0.25;
  constexpr T surfaceSlope = 0.35;
  constexpr T coverDepth = 2.5;

  auto matMap = cs::SmartPointer<ls::MaterialMap>::New();
  matMap->insertNextMaterial(substrateMaterial);

  auto cellSet = CellSet::New();
  cellSet->setCellSetPosition(true);
  cellSet->setCoverMaterial(coverMaterial);
  cellSet->enableEmbeddedBoundaries(withEmbeddedBoundaries);
  cellSet->fromLevelSets({makeTiltedSurface(gridDelta, halfWidth, yMin, yMax,
                                            surfaceY, surfaceSlope)},
                         matMap, coverDepth);
  cellSet->buildNeighborhood();
  return cellSet;
}

bool isSubstrate(T material) {
  return static_cast<int>(std::round(material)) == substrateMaterial;
}

bool isCover(T material) {
  return static_cast<int>(std::round(material)) == coverMaterial;
}

void initializeDirichletCover(CellSet cellSet, std::vector<T> &field,
                              const T boundaryValue) {
  auto materials = cellSet->getScalarData("Material");
  for (int cellId = 0; cellId < static_cast<int>(cellSet->getNumberOfCells());
       ++cellId) {
    if (isCover((*materials)[cellId]))
      field[cellId] = boundaryValue;
  }
}

void runCartesian(CellSet cellSet, std::vector<T> &field, T gridDelta,
                  T duration, T diffusivity) {
  auto materials = cellSet->getScalarData("Material");
  cs::DiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::Explicit);
  solver.setDiffusiveMaterials({substrateMaterial});
  solver.setSourceMaterials({coverMaterial});

  T dt = T(0.2) * gridDelta * gridDelta / (T(2) * D * diffusivity);
  T time = 0.;
  while (time < duration) {
    if (time + dt > duration)
      dt = duration - time;
    solver.step(field, *materials, gridDelta, dt, diffusivity);
    time += dt;
  }
}

void runEmbedded(CellSet cellSet, std::vector<T> &field, T gridDelta,
                 T duration, T diffusivity, T boundaryValue) {
  const auto materials = cellSet->getScalarData("Material");
  cs::BoundaryDiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::GaussSeidel);
  solver.setMaxIterations(1000);
  solver.setTolerance(1.e-10);
  solver.setDiffusiveMaterials({substrateMaterial});
  solver.setSourceMaterials({coverMaterial});
  solver.setBoundaryCondition(
      0, cs::BoundaryCondition<T>::dirichlet(boundaryValue));

  T dt = duration;
  T time = 0.;
  while (time < duration) {
    if (time + dt > duration)
      dt = duration - time;
    solver.step(field, *materials, gridDelta, dt, diffusivity);
    time += dt;
  }
}

void writeComparison(CellSet cellSet, const std::vector<T> &cartesian,
                     const std::vector<T> &embedded) {
  cellSet->addScalarData("cartesian", 0.);
  cellSet->setScalarData("cartesian", cartesian);
  cellSet->addScalarData("embedded", 0.);
  cellSet->setScalarData("embedded", embedded);

  std::vector<T> difference(cartesian.size(), 0.);
  for (std::size_t i = 0; i < difference.size(); ++i)
    difference[i] = embedded[i] - cartesian[i];
  cellSet->addScalarData("embedded_minus_cartesian", 0.);
  cellSet->setScalarData("embedded_minus_cartesian", difference);

  std::vector<T> boundaryPointCount(cartesian.size(), 0.);
  std::vector<T> physicalCutCell(cartesian.size(), 0.);
  for (int cellId = 0; cellId < static_cast<int>(cellSet->getNumberOfCells());
       ++cellId) {
    boundaryPointCount[cellId] =
        static_cast<T>(cellSet->getEmbeddedBoundaryPointIds(cellId).size());
    for (const auto pointId : cellSet->getEmbeddedBoundaryPointIds(cellId)) {
      if (cellSet->getEmbeddedBoundaryPoints()[pointId].levelSetIndex == 0) {
        physicalCutCell[cellId] = 1.;
        break;
      }
    }
  }
  cellSet->addScalarData("embedded_boundary_point_count", 0.);
  cellSet->setScalarData("embedded_boundary_point_count", boundaryPointCount);
  cellSet->addScalarData("physical_interface_cut_cell", 0.);
  cellSet->setScalarData("physical_interface_cut_cell", physicalCutCell);

  cellSet->writeVTU("embeddedBoundaryDiffusion.vtu");
}

void printComparison(CellSet cellSet, const std::vector<T> &cartesian,
                     const std::vector<T> &embedded) {
  auto materials = cellSet->getScalarData("Material");
  T sumCartesian = 0.;
  T sumEmbedded = 0.;
  T l2 = 0.;
  T maxAbs = 0.;
  int substrateCells = 0;

  for (int i = 0; i < static_cast<int>(cartesian.size()); ++i) {
    if (!isSubstrate((*materials)[i]))
      continue;
    const T diff = embedded[i] - cartesian[i];
    sumCartesian += cartesian[i];
    sumEmbedded += embedded[i];
    l2 += diff * diff;
    maxAbs = std::max(maxAbs, std::abs(diff));
    ++substrateCells;
  }

  l2 = substrateCells > 0 ? std::sqrt(l2 / substrateCells) : 0.;
  std::cout << std::setprecision(6) << std::fixed;
  std::cout << "Substrate cells             : " << substrateCells << "\n";
  std::cout << "Embedded boundary points    : "
            << cellSet->getEmbeddedBoundaryPoints().size() << "\n";
  std::cout << "Mean Cartesian concentration: " << sumCartesian / substrateCells
            << "\n";
  std::cout << "Mean embedded concentration : " << sumEmbedded / substrateCells
            << "\n";
  std::cout << "RMS embedded-Cartesian diff : " << l2 << "\n";
  std::cout << "Max embedded-Cartesian diff : " << maxAbs << "\n";
  std::cout << "Wrote embeddedBoundaryDiffusion.vtu\n";
}

void printBoundaryDiagnostics(CellSet cellSet) {
  std::map<unsigned, unsigned> pointsPerLevelSet;
  T maxPlaneError = 0.;
  unsigned physicalInterfacePoints = 0;
  unsigned physicalInterfaceCutCells = 0;

  std::vector<bool> cutCell(cellSet->getNumberOfCells(), false);
  for (const auto &point : cellSet->getEmbeddedBoundaryPoints()) {
    ++pointsPerLevelSet[point.levelSetIndex];
    if (point.levelSetIndex != 0)
      continue;

    constexpr T surfaceY = 0.25;
    constexpr T surfaceSlope = 0.35;
    maxPlaneError =
        std::max(maxPlaneError,
                 std::abs(point.coordinate[1] -
                          (surfaceY + surfaceSlope * point.coordinate[0])));
    ++physicalInterfacePoints;
    if (point.adjacentCell >= 0 &&
        point.adjacentCell < static_cast<int>(cutCell.size())) {
      cutCell[point.adjacentCell] = true;
    }
  }

  for (const auto isCut : cutCell)
    if (isCut)
      ++physicalInterfaceCutCells;

  std::cout << "Boundary points per level set:";
  for (const auto &[levelSet, count] : pointsPerLevelSet)
    std::cout << " L" << levelSet << "=" << count;
  std::cout << "\n";
  std::cout << "Physical interface points   : " << physicalInterfacePoints
            << "\n";
  std::cout << "Physical cut cells          : " << physicalInterfaceCutCells
            << "\n";
  std::cout << "Max tilted-plane error      : " << maxPlaneError << "\n";
}

} // namespace

int main() {
  cs::Logger::setLogLevel(cs::LogLevel::WARNING);

  constexpr T gridDelta = 0.25;
  constexpr T duration = 0.35;
  constexpr T diffusivity = 0.5;
  constexpr T boundaryValue = 1.;

  auto cartesianCellSet = makeCellSet(gridDelta, false);
  auto embeddedCellSet = makeCellSet(gridDelta, true);

  auto cartesian = std::vector<T>(cartesianCellSet->getNumberOfCells(), 0.);
  auto embedded = std::vector<T>(embeddedCellSet->getNumberOfCells(), 0.);
  initializeDirichletCover(cartesianCellSet, cartesian, boundaryValue);
  initializeDirichletCover(embeddedCellSet, embedded, boundaryValue);

  runCartesian(cartesianCellSet, cartesian, gridDelta, duration, diffusivity);
  runEmbedded(embeddedCellSet, embedded, gridDelta, duration, diffusivity,
              boundaryValue);

  writeComparison(embeddedCellSet, cartesian, embedded);
  printBoundaryDiagnostics(embeddedCellSet);
  printComparison(embeddedCellSet, cartesian, embedded);
  return 0;
}
