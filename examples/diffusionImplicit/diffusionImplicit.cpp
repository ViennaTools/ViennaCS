#include <csDenseCellSet.hpp>

#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseLU>

#include "geometry.hpp"

namespace cs = viennacs;
namespace ls = viennals;

using T = double;
constexpr int D = 2;

const int substrateMaterial = 0;
const int maskMaterial = 1;
const int coverMaterial = 2;

struct SolutionData {
  Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
  std::unordered_map<int, int> cellMapping; // cellSet -> solution index
  unsigned numCells = 0;
  Eigen::SparseMatrix<T> systemMatrix;
  Eigen::Matrix<T, Eigen::Dynamic, 1> rhs;
};

template <typename Material> bool isMaterial(Material x, int material) {
  return static_cast<int>(x) == material;
}

template <typename Material>
bool isDirichletBoundary(std::array<T, 3> center, Material material,
                         cs::util::Parameters &params) {
  return isMaterial(material, coverMaterial) &&
         center[D - 1] <
             params.get("substrateHeight") + params.get("gridDelta");
}

void addConcentration(cs::DenseCellSet<T, D> &cellSet,
                      cs::util::Parameters &params) {
  // Add quantity to be diffused (on top of the cell material)
  auto concentration = cellSet.addScalarData("dopant", 0.);
  auto materials = cellSet.getScalarData("Material");

  const T boundaryValue = params.get("boundaryValue");

  // Boundary condition: constant concentration at the top (outside)
#pragma omp parallel for
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isDirichletBoundary(cellSet.getCellCenter(i), (*materials)[i],
                            params)) {
      (*concentration)[i] = boundaryValue;
    }
  }
}

// Implicit diffusion time-step (backward Euler)
void solveDiffusionStep(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet,
                        cs::util::Parameters &params, T dt) {
  auto data = cellSet.getScalarData("dopant");

  sol.rhs = sol.solver.solve(sol.rhs);

  if (sol.solver.info() != Eigen::Success) {
    cs::Logger::getInstance().addError("Solving failed.").print();
  }

  // write results to cellSet
  for (const auto &i : sol.cellMapping) {
    (*data)[i.first] = sol.rhs[i.second];
  }
}

void initCellMapping(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet,
                     cs::util::Parameters &params) {
  auto materials = cellSet.getScalarData("Material");
  sol.numCells = 0;
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isMaterial((*materials)[i], substrateMaterial) ||
        isDirichletBoundary(cellSet.getCellCenter(i), (*materials)[i],
                            params)) {
      sol.cellMapping[i] = sol.numCells++;
    }
  }
}

void assembleMatrix(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet,
                    cs::util::Parameters &params, T dt) {
  auto materials = cellSet.getScalarData("Material");
  int material = substrateMaterial;
  const T dx = cellSet.getGridDelta();
  const T dtdx2 = dt / (dx * dx);

  std::vector<Eigen::Triplet<T>> triplets; // (i,j,value) matrix entries
  triplets.reserve(2 * D * sol.numCells);

  for (const auto &ids : sol.cellMapping) {
    auto i = ids.first;
    auto cellMapIdx = ids.second;
    bool onBoundary =
        isDirichletBoundary(cellSet.getCellCenter(i), (*materials)[i], params);

    // Boundary conditions
    if (onBoundary) {
      triplets.push_back({cellMapIdx, cellMapIdx, 1.}); // diagonal
      continue;
    }

    const auto &neighbors = cellSet.getNeighbors(i);
    T D_sum = 0;
    for (auto n : neighbors) {
      if (n >= 0 && !isMaterial((*materials)[n], maskMaterial)) {
        T D_loc = params.get("diffusionCoefficient");
        D_sum += D_loc;
        triplets.emplace_back(cellMapIdx, sol.cellMapping[n], -dtdx2 * D_loc);
      }
    }

    triplets.emplace_back(cellMapIdx, cellMapIdx, 1. + dtdx2 * D_sum);
  }

  sol.systemMatrix.resize(sol.numCells, sol.numCells);
  sol.systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
  sol.systemMatrix.makeCompressed();
}

void assembleRHS(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet,
                 cs::util::Parameters &params) {
  auto concentration = cellSet.getScalarData("dopant");
  sol.rhs.resize(sol.numCells);

  // Add boundary conditions
  for (const auto &ids : sol.cellMapping) {
    sol.rhs[ids.second] = (*concentration)[ids.first];
  }
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

  SolutionData solution;

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

  initCellMapping(solution, cellSet, params);
  assembleRHS(solution, cellSet, params);

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

  assembleMatrix(solution, cellSet, params, dt);
  solution.solver.compute(solution.systemMatrix);

  if (solution.solver.info() != Eigen::Success) {
    cs::Logger::getInstance().addError("Decomposition failed.").print();
  }

  T time = 0.;
  while (time < duration) {

    if (time + dt > duration) {
      dt = duration - time;
      assembleMatrix(solution, cellSet, params, dt);
      solution.solver.compute(solution.systemMatrix);
    }

    solveDiffusionStep(solution, cellSet, params, dt);

    time += dt;
  }

  // saveVolumeMesh("final", levelSets, matMap);
  cellSet.writeVTU("final.vtu");
}
