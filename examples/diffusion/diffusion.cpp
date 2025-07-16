#include <csDenseCellSet.hpp>

#include "geometry.hpp"

namespace cs = viennacs;
namespace ls = viennals;

using T = double;
constexpr int D = 2;

const int substrateMaterial = 0;
const int maskMaterial = 1;
const int coverMaterial = 2;

template <typename Material> bool isMaterial(Material x, int material) {
  return static_cast<int>(x) == material;
}

void addConcentration(cs::DenseCellSet<T, D> &cellSet,
                      cs::util::Parameters &params) {
  // Add quantity to be diffused (on top of the cell material)
  auto concentration = cellSet.addScalarData("dopant", 0.);
  auto materials = cellSet.getScalarData("Material");

  const T boundaryValue = params.get("boundaryValue");
  const T heightLimit = params.get("substrateHeight") + params.get("gridDelta");

  // Boundary condition: constant concentration at the top (outside)
#pragma omp parallel for
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isMaterial((*materials)[i], coverMaterial) &&
        cellSet.getCellCenter(i)[D - 1] < heightLimit) {
      (*concentration)[i] = boundaryValue;
    }
  }
}

// Explicit diffusion time-step (forward Euler)
void solveDiffusionStep(cs::DenseCellSet<T, D> &cellSet,
                        cs::util::Parameters &params, T dt) {
  auto data = cellSet.getScalarData("dopant");
  auto materials = cellSet.getScalarData("Material");
  std::vector<T> solution(data->size(), 0.);
  const T dx = params.get("gridDelta");
  const T C = dt * params.get("diffusionCoefficient") / (dx * dx);

#pragma omp parallel for
  for (int e = 0; e < data->size(); e++) {
    if (!isMaterial((*materials)[e], substrateMaterial)) {
      solution[e] = (*data)[e];
      continue;
    }

    auto coord = cellSet.getCellCenter(e);
    int numNeighbors = 0;

    const auto &cellNeighbors = cellSet.getNeighbors(e);
    for (auto n : cellNeighbors) {
      if (n == -1 || isMaterial((*materials)[n], maskMaterial))
        continue;

      solution[e] += (*data)[n];
      numNeighbors++;
    }

    // Diffusion
    solution[e] = (*data)[e] + C * (solution[e] - numNeighbors * (*data)[e]);
  }
  *data = std::move(solution);
}

int main(int argc, char **argv) {
  cs::Logger::setLogLevel(cs::LogLevel::INTERMEDIATE);

  cs::util::Parameters params;
  if (argc > 1) {
    params.readConfigFile(argv[1]);
  } else {
    std::cout << "Usage: " << argv[0] << " <config file>" << std::endl;
    return 1;
  }
  omp_set_num_threads(params.get<int>("numThreads"));

  auto matMap = cs::SmartPointer<ls::MaterialMap>::New();
  auto levelSets = geometry::makeStructure<T, D>(
      params, matMap, substrateMaterial, maskMaterial);

  cs::DenseCellSet<T, D> cellSet;
  T depth = params.get("substrateHeight") + params.get("coverHeight") + 10.;
  cellSet.setCellSetPosition(true); // isAboveSurface
  cellSet.setCoverMaterial(coverMaterial);
  cellSet.fromLevelSets(levelSets, matMap, depth);

  // We need neighborhood information for solving the diffusion equation
  cellSet.buildNeighborhood();

  addConcentration(cellSet, params);
  cellSet.writeVTU("initial.vtu");

  if (params.get("velocity") != 0.) {
    const T stability =
        2 * params.get("diffusionCoefficient") / params.get("velocity");
    std::cout << "Stability: " << stability << std::endl;
    if (0.5 * stability <= params.get("gridDelta"))
      std::cout << "Unstable parameters. Reduce grid spacing!" << std::endl;
  }

  T duration = params.get("duration");
  T dx = params.get("gridDelta");
  T dt = std::min(dx * dx / (params.get("diffusionCoefficient") * 2 * D) *
                      params.get("timeStabilityFactor"),
                  duration);
  T time = 0.;
  while (time < duration) {

    if (time + dt > duration)
      dt = duration - time;

    solveDiffusionStep(cellSet, params, dt);

    time += dt;
  }

  // saveVolumeMesh("final", levelSets, matMap);
  cellSet.writeVTU("final.vtu");
}
