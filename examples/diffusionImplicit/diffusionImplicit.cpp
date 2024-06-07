#include <csDenseCellSet.hpp>
#include <lsBooleanOperation.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseLU>

using T = double;
constexpr int D = 2;

struct Parameters {
  // Domain
  T gridDelta = 1.5; // nm
  T xExtent = 80.0;  // nm
  T yExtent = 80.0;  // nm

  // Geometry
  T substrateHeight = 50.; // nm
  T coverHeight = 30.;
  T maskHeight = 10.;
  T holeRadius = xExtent / 6.;

  // Process
  T duration = 20.;
  T diffusionCoefficient = 1.; // nm²/s
  T velocity = 0.;             // Advection
  T boundaryValue = 1.;
  T timeStabilityFactor = 0.95;

  int substrateMaterial = 0;
  int maskMaterial = 1;
  int coverMaterial = 2;
};

struct SolutionData {
  Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
  std::unordered_map<int, int> cellMapping; // cellSet -> solution index
  unsigned numCells = 0;
  Eigen::SparseMatrix<T> systemMatrix;
  Eigen::Matrix<T, Eigen::Dynamic, 1> rhs;
};

namespace cs = viennacs;
namespace ls = viennals;

using levelSetType = cs::SmartPointer<ls::Domain<T, D>>;
using levelSetsType = std::vector<levelSetType>;
using materialMapType = cs::SmartPointer<ls::MaterialMap>;

void addLevelSet(levelSetsType &levelSets, levelSetType levelSet,
                 materialMapType matMap, int material,
                 bool wrapLowerLevelSet = true) {
  if (!levelSets.empty() && wrapLowerLevelSet) {
    ls::BooleanOperation<T, D>(levelSet, levelSets.back(),
                               ls::BooleanOperationEnum::UNION)
        .apply();
  }

  levelSets.push_back(levelSet);
  matMap->insertNextMaterial(material);
}

void makePlane(levelSetType &domain, const T *origin, const T *normal) {
  ls::MakeGeometry<T, D>(domain,
                         cs::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
      .apply();
}

void makeBox(levelSetType &domain, const T *minPoint, const T *maxPoint) {
  ls::MakeGeometry<T, D>(
      domain, cs::SmartPointer<ls::Box<T, D>>::New(minPoint, maxPoint))
      .apply();
}

auto makeStructure(const Parameters &params, materialMapType matMap) {
  const T gridDelta = params.gridDelta;
  const T substrateHeight = params.substrateHeight;
  const T coverHeight = params.coverHeight;
  const T maskHeight = params.maskHeight;
  const T holeRadius = params.holeRadius;
  ls::BoundaryConditionEnum<D> boundaryConds[D] = {
      ls::BoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY};
  boundaryConds[D - 1] = ls::BoundaryConditionEnum<D>::INFINITE_BOUNDARY;
  T bounds[2 * D] = {-params.xExtent / 2., params.xExtent / 2.,
                     -params.yExtent / 2., params.yExtent / 2.};
  bounds[2 * D - 2] = 0.;
  bounds[2 * D - 1] = substrateHeight + maskHeight + gridDelta;

  T origin[D] = {};
  T normal[D] = {};
  normal[D - 1] = 1.;

  levelSetsType levelSets;

  // Substrate
  origin[D - 1] = 0.;
  auto bottom = levelSetType::New(bounds, boundaryConds, gridDelta);
  makePlane(bottom, origin, normal);
  addLevelSet(levelSets, bottom, matMap, params.substrateMaterial);

  origin[D - 1] = substrateHeight;
  auto substrate = levelSetType::New(bounds, boundaryConds, gridDelta);
  makePlane(substrate, origin, normal);
  addLevelSet(levelSets, substrate, matMap, params.substrateMaterial);

  // Mask
  if (maskHeight > 0.) {
    auto mask = levelSetType::New(bounds, boundaryConds, gridDelta);
    origin[D - 1] = substrateHeight + maskHeight;
    makePlane(mask, origin, normal);

    auto maskAdd = levelSetType::New(bounds, boundaryConds, gridDelta);

    T minPoint[D] = {-holeRadius, -holeRadius};
    T maxPoint[D] = {holeRadius, holeRadius};
    minPoint[D - 1] = substrateHeight - gridDelta;
    maxPoint[D - 1] = substrateHeight + maskHeight + gridDelta;

    makeBox(maskAdd, minPoint, maxPoint);

    ls::BooleanOperation<T, D>(mask, maskAdd,
                               ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    addLevelSet(levelSets, mask, matMap, params.maskMaterial);
  }

  return levelSets;
}

template <typename Material> bool isMaterial(Material x, int material) {
  return static_cast<int>(x) == material;
}

template <typename Material>
bool isDirichletBoundary(std::array<T, 3> center, Material material,
                         const Parameters &params) {
  return isMaterial(material, params.coverMaterial) &&
         center[D - 1] < params.substrateHeight + params.gridDelta;
}

void addConcentration(cs::DenseCellSet<T, D> &cellSet,
                      const Parameters &params) {
  // Add quantity to be diffused (on top of the cell material)
  auto concentration = cellSet.addScalarData("dopant", 0.);
  auto materials = cellSet.getScalarData("Material");

  // Boundary condition: constant concentration at the top (outside)
#pragma omp parallel for
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isDirichletBoundary(cellSet.getCellCenter(i), (*materials)[i],
                            params)) {
      (*concentration)[i] = params.boundaryValue;
    }
  }
}

// Implicit diffusion time-step (backward Euler)
void solveDiffusionStep(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet,
                        const Parameters &params, T dt) {
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
                     const Parameters &params) {
  auto materials = cellSet.getScalarData("Material");
  sol.numCells = 0;
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isMaterial((*materials)[i], params.substrateMaterial) ||
        isDirichletBoundary(cellSet.getCellCenter(i), (*materials)[i],
                            params)) {
      sol.cellMapping[i] = sol.numCells++;
    }
  }
}

void assembleMatrix(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet,
                    const Parameters &params, T dt) {
  auto materials = cellSet.getScalarData("Material");
  int material = params.substrateMaterial;
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
      if (n >= 0 && !isMaterial((*materials)[n], params.maskMaterial)) {
        T D_loc = params.diffusionCoefficient;
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
                 const Parameters &params) {
  auto concentration = cellSet.getScalarData("dopant");
  sol.rhs.resize(sol.numCells);

  // Add boundary conditions
  for (const auto &ids : sol.cellMapping) {
    sol.rhs[ids.second] = (*concentration)[ids.first];
  }
}

void saveVolumeMesh(std::string name, levelSetsType &levelSets,
                    materialMapType matMap) {
  ls::WriteVisualizationMesh<T, D> writer;
  writer.setFileName(name);
  for (const auto &ls : levelSets)
    writer.insertNextLevelSet(ls);
  writer.setMaterialMap(matMap);
  writer.apply();
}

int main(int argc, char **argv) {
  omp_set_num_threads(4);
  cs::Logger::setLogLevel(cs::LogLevel::INTERMEDIATE);

  Parameters params;
  SolutionData solution;

  auto matMap = materialMapType::New();
  auto levelSets = makeStructure(params, matMap);

  cs::DenseCellSet<T, D> cellSet;
  T depth = params.substrateHeight + params.coverHeight + 10.;
  cellSet.setCellSetPosition(true); // isAboveSurface
  cellSet.setCoverMaterial(params.coverMaterial);
  cellSet.fromLevelSets(levelSets, matMap, depth);

  addConcentration(cellSet, params);
  // We need neighborhood information for solving the diffusion equation
  cellSet.buildNeighborhood();
  cellSet.writeVTU("initial.vtu");

  initCellMapping(solution, cellSet, params);
  assembleRHS(solution, cellSet, params);

  if (params.velocity != 0.) {
    const T stability = 2 * params.diffusionCoefficient / params.velocity;
    std::cout << "Stability: " << stability << std::endl;
    if (0.5 * stability <= params.gridDelta)
      std::cout << "Unstable parameters. Reduce grid spacing!" << std::endl;
  }

  T duration = params.duration;
  T dt = std::min(params.gridDelta * params.gridDelta /
                      (params.diffusionCoefficient * 2 * D) *
                      params.timeStabilityFactor,
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
