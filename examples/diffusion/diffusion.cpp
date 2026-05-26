#include <csDenseCellSet.hpp>
#include <csDiffusionSolver.hpp>

#include "geometry.hpp"

namespace cs = viennacs;
namespace ls = viennals;

using T = double;
constexpr int D = 2;

const int substrateMaterial = 0;
const int maskMaterial = 1;
const int coverMaterial = 2;

int main(int argc, char **argv) {
  cs::Logger::setLogLevel(cs::LogLevel::INTERMEDIATE);

  cs::util::Parameters params;
  if (argc > 1) {
    params.readConfigFile(argv[1]);
  } else {
    params.readConfigFile("config.txt");
    if (params.m.empty()) {
      std::cout << "No configuration file provided!\n"
                << "Usage: " << argv[0] << " <config file>\n";
      return 1;
    }
  }

  omp_set_num_threads(params.get<int>("numThreads"));

  auto matMap = cs::SmartPointer<ls::MaterialMap>::New();
  auto levelSets = geometry::makeStructure<T, D>(
      params, matMap, substrateMaterial, maskMaterial);

  auto cellSet = cs::SmartPointer<cs::DenseCellSet<T, D>>::New();
  const T depth = params.get("substrateHeight") + params.get("coverHeight") + 10.;
  cellSet->setCellSetPosition(true);
  cellSet->setCoverMaterial(coverMaterial);
  cellSet->fromLevelSets(levelSets, matMap, depth);
  cellSet->buildNeighborhood();

  // Dirichlet BC: fixed concentration in cover cells at the surface boundary
  auto concentration = cellSet->addScalarData("dopant", 0.);
  auto materials = cellSet->getScalarData("Material");
  const T boundaryValue = params.get("boundaryValue");
  const T heightLimit = params.get("substrateHeight") + params.get("gridDelta");
#pragma omp parallel for
  for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
    if (static_cast<int>((*materials)[i]) == coverMaterial &&
        cellSet->getCellCenter(i)[D - 1] < heightLimit)
      (*concentration)[i] = boundaryValue;
  }
  cellSet->writeVTU("initial.vtu");

  const T diffCoeff = params.get("diffusionCoefficient");
  const T dx = params.get("gridDelta");
  const T duration = params.get("duration");

  if (params.get("velocity") != 0.) {
    const T stability = 2 * diffCoeff / params.get("velocity");
    std::cout << "Stability: " << stability << "\n";
    if (0.5 * stability <= dx)
      std::cout << "Unstable parameters. Reduce grid spacing!\n";
  }

  T dt = std::min(dx * dx / (diffCoeff * 2 * D) *
                      params.get("timeStabilityFactor"),
                  duration);

  cs::DiffusionSolver<T, D> solver;
  solver.setCellSet(cellSet);
  solver.setMode(cs::DiffusionSolverMode::Explicit);
  solver.setDiffusiveMaterials({substrateMaterial});
  solver.setSourceMaterials({coverMaterial}); // fixed-value Dirichlet neighbours
  solver.setBlockedMaterials({maskMaterial});

  T time = 0.;
  while (time < duration) {
    if (time + dt > duration)
      dt = duration - time;
    solver.step(*concentration, *materials, dx, dt, diffCoeff);
    time += dt;
  }

  cellSet->writeVTU("final.vtu");
}
