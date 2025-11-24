/**
  Electric Field Driven Oxidation Model
  
  Adapts the diffusion example to model Si -> SiO2 oxidation.
  - Oxidant diffuses from Gas through Oxide to Silicon.
  - Gas/Oxide interface has a prescribed concentration (Dirichlet).
  - Reaction Si + O -> SiO2 depends on E-field magnitude.
  - Tracks oxide growth via "oxideFraction" scalar.
*/

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream> // Added for std::cout

#include <csDenseCellSet.hpp>

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

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

// -----------------------------------------------------------------------------
// DATA LOADING WITH SPATIAL MAPPING
// -----------------------------------------------------------------------------

void addElectricField(cs::DenseCellSet<T, D> &cellSet,
                      cs::util::Parameters &params) {
  std::string csvFileName = params.get<std::string>("EfieldFile");
  const auto numCells = cellSet.getNumberOfCells();
  if (numCells == 0) return;

  std::ifstream in(csvFileName);
  if (!in.good()) {
    cs::Logger::getInstance()
        .addError("addElectricField: cannot open file " + csvFileName)
        .print();
    return;
  }

  // Initialize vector data
  auto E = cellSet.getVectorData("Efield");
  if (!E) {
    E = cellSet.addVectorData("Efield", {0., 0., 0.});
  } else {
    E->assign(numCells, {0., 0., 0.});
  }

  // 1. Load CSV into a flat buffer
  std::vector<std::array<T, 3>> rawData;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line[0] == '%') continue;
    std::stringstream ss(line);
    std::string segment;
    std::array<T, 3> val = {0., 0., 0.};
    int component = 0;
    while (std::getline(ss, segment, ',') && component < 3) {
        try {
            val[component++] = std::stod(segment);
        } catch (...) {}
    }
    if(component > 0) rawData.push_back(val);
  }

  // 2. Reconstruct Grid Info from Params
  T dx = params.get("gridDelta");
  T xExtent = params.get("xExtent");
  // Assuming E-field was generated for the substrate region height
  T yExtent = params.get("substrateHeight"); 

  int nx = static_cast<int>(xExtent / dx);
  int ny = static_cast<int>(yExtent / dx);

  if (rawData.size() < nx * ny) {
      cs::Logger::getInstance()
          .addWarning("CSV has fewer rows than expected grid size. Mapping may fail.")
          .print();
  }

  // 3. Map Cell Centers to E-field Grid
  int mappedCount = 0;
  for (unsigned int i = 0; i < numCells; ++i) {
      auto center = cellSet.getCellCenter(i);
      
      // Transform Simulation Coordinates to Generator Coordinates
      // Sim X: [-xExtent/2, xExtent/2] -> Gen X: [0, xExtent]
      T gen_x = center[0] + (xExtent / 2.0);
      T gen_y = center[1]; // Y matches (0 at bottom)

      int ix = static_cast<int>(gen_x / dx);
      int iy = static_cast<int>(gen_y / dx);

      // Check bounds
      if (ix >= 0 && ix < nx && iy >= 0 && iy < ny) {
          // Index in flat array (assuming Y-major from generator loops: for y... for x...)
          // Generator: for j, y... for i, x... -> rawData index = j * nx + i
          size_t gridIdx = iy * nx + ix;
          
          if (gridIdx < rawData.size()) {
              (*E)[i] = rawData[gridIdx];
              mappedCount++;
          }
      }
  }

  cs::Logger::getInstance()
      .addInfo("Mapped E-field to " + std::to_string(mappedCount) + " cells.")
      .print();
}

void initOxidationFields(cs::DenseCellSet<T, D> &cellSet,
                         cs::util::Parameters &params) {
  auto oxidant = cellSet.addScalarData("oxidant", 0.);
  auto oxideFraction = cellSet.addScalarData("oxideFraction", 0.);
  auto materials = cellSet.getScalarData("Material");
  const T boundaryValue = params.get("boundaryValue");

  #pragma omp parallel for
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isMaterial((*materials)[i], coverMaterial)) {
      (*oxidant)[i] = boundaryValue;
    }
  }
}

void initCellMapping(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet) {
  auto materials = cellSet.getScalarData("Material");
  sol.numCells = 0;
  sol.cellMapping.clear();
  
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isMaterial((*materials)[i], substrateMaterial)) {
      sol.cellMapping[i] = sol.numCells++;
    }
  }
}

// -----------------------------------------------------------------------------
// PHYSICS KERNELS
// -----------------------------------------------------------------------------

T getReactionRate(const std::array<T, 3>& E, T oxideFrac, cs::util::Parameters &params) {
  const T k_base = params.get("reactionRateConstant"); 
  const T alpha = params.get("eFieldInfluence");
  
  T E_mag = std::sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]);
  T availableSi = std::max(0.0, 1.0 - oxideFrac);
  
  return k_base * (1.0 + alpha * E_mag) * availableSi;
}

T getDiffusivity(int material, T oxideFrac, cs::util::Parameters &params) {
  if (isMaterial(material, maskMaterial)) return 0.0;
  if (isMaterial(material, coverMaterial)) return 0.0; // Not used in solve directly
  
  const T D_ox = params.get("diffusivityOxide");
  return D_ox * oxideFrac;
}

// -----------------------------------------------------------------------------
// SOLVER
// -----------------------------------------------------------------------------

void assembleMatrix(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet,
                    cs::util::Parameters &params, T dt) {
  auto materials = cellSet.getScalarData("Material");
  auto oxideFractions = cellSet.getScalarData("oxideFraction");
  auto EFields = cellSet.getVectorData("Efield");
  
  const T dx = cellSet.getGridDelta();
  const T dtdx2 = dt / (dx * dx);

  std::vector<Eigen::Triplet<T>> triplets;
  triplets.reserve(2 * D * sol.numCells);

  for (const auto &ids : sol.cellMapping) {
    auto i = ids.first;
    auto row = ids.second;
    
    // 1. Reaction Sink
    T reactionSink = 0.0;
    std::array<T,3> E = (EFields) ? (*EFields)[i] : std::array<T,3>{0.,0.,0.};
    T k = getReactionRate(E, (*oxideFractions)[i], params);
    reactionSink = dt * k; 

    // 2. Diffusion
    const auto &neighbors = cellSet.getNeighbors(i);
    T D_center = getDiffusivity((*materials)[i], (*oxideFractions)[i], params);
    T D_sum_neighbors = 0.0;

    for (auto n : neighbors) {
      if (n < 0) continue; // Boundary

      int mat_n = static_cast<int>((*materials)[n]);
      if (isMaterial(mat_n, maskMaterial)) continue; // Impermeable

      T D_eff = 0.0;

      if (sol.cellMapping.count(n)) {
        // Substrate Neighbor
        T D_neighbor = getDiffusivity(mat_n, (*oxideFractions)[n], params);
        D_eff = 0.5 * (D_center + D_neighbor);
        
        triplets.emplace_back(row, sol.cellMapping[n], -dtdx2 * D_eff);
        
      } else if (isMaterial(mat_n, coverMaterial)) {
        // Gas Neighbor (Dirichlet)
        D_eff = D_center; 
        // Flux in handled in RHS
      }
      
      D_sum_neighbors += D_eff;
    }

    // Diagonal
    triplets.emplace_back(row, row, 1. + (dtdx2 * D_sum_neighbors) + reactionSink);
  }

  sol.systemMatrix.resize(sol.numCells, sol.numCells);
  sol.systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
  sol.systemMatrix.makeCompressed();
}

void assembleRHS(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet, 
                 cs::util::Parameters &params, T dt) {
  auto oxidant = cellSet.getScalarData("oxidant");
  auto materials = cellSet.getScalarData("Material");
  auto oxideFractions = cellSet.getScalarData("oxideFraction");
  
  const T boundaryValue = params.get("boundaryValue");
  const T dx = cellSet.getGridDelta();
  const T dtdx2 = dt / (dx * dx);

  sol.rhs.resize(sol.numCells);

  for (const auto &ids : sol.cellMapping) {
    auto i = ids.first;
    auto row = ids.second;
    
    T val = (*oxidant)[i]; // Old value

    // Add Flux from Gas neighbors
    const auto &neighbors = cellSet.getNeighbors(i);
    T D_center = getDiffusivity((*materials)[i], (*oxideFractions)[i], params);

    for (auto n : neighbors) {
      if (n >= 0 && isMaterial((*materials)[n], coverMaterial)) {
         T D_interface = params.get("diffusivityOxide"); // Use D_ox for the interface flux
         val += dtdx2 * D_interface * boundaryValue;
      }
    }
    
    sol.rhs[row] = val;
  }
}

void solveDiffusionStep(SolutionData &sol, cs::DenseCellSet<T, D> &cellSet) {
  auto oxidant = cellSet.getScalarData("oxidant");

  sol.rhs = sol.solver.solve(sol.rhs);

  if (sol.solver.info() != Eigen::Success) {
    cs::Logger::getInstance().addError("Solving failed.").print();
    return;
  }

  for (const auto &i : sol.cellMapping) {
    (*oxidant)[i.first] = sol.rhs[i.second];
  }
}

void updateOxideFraction(cs::DenseCellSet<T, D> &cellSet, 
                         cs::util::Parameters &params, T dt) {
    auto oxidant = cellSet.getScalarData("oxidant");
    auto oxideFractions = cellSet.getScalarData("oxideFraction");
    auto materials = cellSet.getScalarData("Material");
    auto EFields = cellSet.getVectorData("Efield");

    T conversionFactor = params.get("oxideConversionRate");

    #pragma omp parallel for
    for(int i=0; i<cellSet.getNumberOfCells(); ++i) {
        if(isMaterial((*materials)[i], substrateMaterial)) {
            T C = (*oxidant)[i];
            if(C > 1e-12) {
                std::array<T,3> E = (EFields) ? (*EFields)[i] : std::array<T,3>{0.,0.,0.};
                T k = getReactionRate(E, (*oxideFractions)[i], params);
                
                T dFrac = k * C * dt * conversionFactor;
                
                (*oxideFractions)[i] += dFrac;
                if((*oxideFractions)[i] > 1.0) (*oxideFractions)[i] = 1.0;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// MAIN
// -----------------------------------------------------------------------------

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
  cellSet.setCellSetPosition(true); 
  cellSet.setCoverMaterial(coverMaterial);
  cellSet.fromLevelSets(levelSets, matMap, depth);

  cellSet.buildNeighborhood();

  // 1. Load E-Field (Now passing params for spatial mapping)
  addElectricField(cellSet, params);

  // 2. Initialize Fields
  initOxidationFields(cellSet, params);
  cellSet.writeVTU("oxidation_initial.vtu");

  initCellMapping(solution, cellSet);

  // Time Stepping
  T duration = params.get("duration");
  T dx = params.get("gridDelta");
  T D_oxide = params.get("diffusivityOxide");
  
  T dt = std::min(dx * dx / (D_oxide * 2 * D) * params.get("timeStabilityFactor"), duration);
  dt *= 5.0; 
  
  T time = 0.;
  int step = 0;
  
  std::cout << "Starting oxidation simulation..." << std::endl;
  std::cout << "Initial dt: " << dt << std::endl;

  while (time < duration) {
    if (time + dt > duration) dt = duration - time;

    // A. Assemble Matrix 
    assembleMatrix(solution, cellSet, params, dt);
    solution.solver.compute(solution.systemMatrix);

    if (solution.solver.info() != Eigen::Success) {
       cs::Logger::getInstance().addError("Decomposition failed.").print();
       break;
    }

    // B. Solve Diffusion (Pass dt/params for boundary flux)
    assembleRHS(solution, cellSet, params, dt);
    solveDiffusionStep(solution, cellSet);

    // C. Update Oxide State
    updateOxideFraction(cellSet, params, dt);

    time += dt;
    step++;
    
    if(step % 10 == 0) {
        std::cout << "Step " << step << " Time: " << time << std::endl;
        cellSet.writeVTU("oxidation_step_" + std::to_string(step) + ".vtu");
    }
  }

  cellSet.writeVTU("oxidation_final.vtu");
  std::cout << "Done." << std::endl;
}