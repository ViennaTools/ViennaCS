/**
  Electric Field Driven Oxidation Model

  Adapts the diffusion example to model Si -> SiO2 oxidation.
  - Oxidant diffuses from Gas through Oxide to Silicon.
  - Gas/Oxide interface has a prescribed concentration (Dirichlet).
  - Reaction Si + O -> SiO2 depends on E-field magnitude.
  - Tracks oxide growth via "oxideFraction" scalar.
*/

#include <algorithm>
#include <fstream>
#include <iostream> // Added for std::cout
#include <omp.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <csDenseCellSet.hpp>

#include "geometry.hpp"

namespace cs = viennacs;
namespace ls = viennals;

using T = double;
constexpr int D = 3;

// Material IDs used in the simulation
enum MaterialIDs {
  substrate = 0,
  oxide = 1,
  mask = 2,
  ambient = 3
};

template <typename Material> bool isMaterial(Material x, MaterialIDs material) {
  return static_cast<int>(x) == static_cast<int>(material);
}

// -----------------------------------------------------------------------------
// DATA LOADING WITH SPATIAL MAPPING
// -----------------------------------------------------------------------------

class OxidationSimulation {
public:
  OxidationSimulation(const std::string &configFile) {
    params.readConfigFile(configFile);
    omp_set_num_threads(params.get<int>("numThreads"));

    // Generate Geometry and Cell Set
    geometry::makeCellSet<T, D>(cellSet, params, substrate, mask, ambient);
    cellSet.buildNeighborhood();

    // Initialize Fields
    initOxidationFields();
    cellSet.writeVTU("oxidation_initial.vtu");
  }

  void apply() {
    T duration = params.get("duration");
    T dx = params.get("gridDelta");
    T D_oxide = params.get("oxidantDiffusivity");

    // Estimate stable time step (Von Neumann stability for diffusion)
    T dt = std::min(dx * dx / (D_oxide * 2 * D) *
                        params.get("timeStabilityFactor"),
                    duration);

    T time = 0.;
    int step = 0;

    std::cout << "Starting oxidation simulation..." << std::endl;
    std::cout << "Initial dt: " << dt << std::endl;

    while (time < duration) {
      if (time + dt > duration)
        dt = duration - time;

      // 0. Update Active Cells (Dynamic Domain)
      updateActiveCells();

      // Update E-field (Database access simulation)
      updateElectricField(time);

      // Pre-calculate reaction rates for the current time step
      std::vector<T> reactionRates(cellSet.getNumberOfCells(), 0.0);
      if (reactionRates.size() != cellSet.getNumberOfCells()) {
        reactionRates.resize(cellSet.getNumberOfCells());
      }
      calculateReactionRates(reactionRates);

      // Calculate adaptive time step based on Diffusion and Reaction stability
      T max_k = 0.0;
#pragma omp parallel for reduction(max : max_k)
      for (size_t i = 0; i < reactionRates.size(); ++i) {
        if (reactionRates[i] > max_k)
          max_k = reactionRates[i];
      }

      T dt_diff = dx * dx / (D_oxide * 2 * D);
      T dt_react = (max_k > 1e-12) ? (1.0 / max_k) : duration;

      T dt = std::min(dt_diff, dt_react) * params.get("timeStabilityFactor");

      if (time + dt > duration)
        dt = duration - time;

      if (step == 0) {
        std::cout << "Initial dt: " << dt << std::endl;
      }

      // A. Solve Diffusion (Explicit)
      solveDiffusionExplicit(dt, reactionRates);

      // C. Update Oxide State
      updateOxideFraction(dt, reactionRates);

      time += dt;
      step++;

      // Output intermediate results
      if (step % 10 == 0) {
        std::cout << "Step " << step << " Time: " << time << std::endl;
        cellSet.writeVTU("oxidation_step_" + std::to_string(step) + ".vtu");
      }
    }

    // Final Output
    cellSet.writeVTU("oxidation_final.vtu");
    std::cout << "Done." << std::endl;
  }

private:
  cs::util::Parameters params;
  cs::DenseCellSet<T, D> cellSet;
  std::vector<int> activeCells;
  std::vector<bool> isActiveCell;
  std::vector<T> nextOxidant;
  std::vector<T> reactionRates;
  std::vector<T> electricField1D;

  // Updates E-field vectors. In this toy model, it reads from a CSV.
  void updateElectricField(T time) {
    electricField1D.clear();
    std::string csvFileName = params.get<std::string>("EfieldFile");

    std::ifstream EfieldFile(csvFileName);
    if (!EfieldFile.good()) {
      cs::Logger::getInstance()
          .addError("updateElectricField: cannot open file " + csvFileName)
          .print();
      return;
    }

    // 1. Load CSV into a flat buffer
    std::string line;
    while (std::getline(EfieldFile, line)) {
      if (line.empty() || line[0] == '#' || line[0] == '%')
        continue;
      std::stringstream ss(line);
      std::string segment;
      std::vector<T> row;
      while (std::getline(ss, segment, ',')) {
        try {
          row.push_back(std::stod(segment));
        } catch (...) {}
      }

      // For D=2, expect x, E (2 cols). For D=3, expect x, y, E (3 cols).
      if (D == 2 && row.size() >= 2) {
        electricField1D.push_back(row[1]);
      } else if (D == 3 && row.size() >= 3) {
        electricField1D.push_back(row[2]);
      } else if (!row.empty()) {
        // Fallback or error handling if needed, currently just taking last
        electricField1D.push_back(row.back());
      }
    }

    cs::Logger::getInstance()
        .addInfo("Loaded " + std::to_string(electricField1D.size()) +
                 " 1D E-field data points at time " + std::to_string(time))
        .print();
  }

  // Initializes scalar fields for the simulation
  void initOxidationFields() {
    auto oxidant = cellSet.addScalarData("oxidant", 0.);
    auto oxideFraction = cellSet.addScalarData("oxideFraction", 0.);
    auto materials = cellSet.getScalarData("Material");
    const T ambientOxidant = params.get("ambientOxidant");

    // Set initial boundary condition: Cover material (Gas) has fixed oxidant
    // concentration
#pragma omp parallel for
    for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
      if (isMaterial((*materials)[i], ambient)) {
        (*oxidant)[i] = ambientOxidant;
      }
    }
  }

  // Updates the mapping from global cell IDs to linear system row indices.
  void updateActiveCells() {
    auto materials = cellSet.getScalarData("Material");
    auto oxideFractions = cellSet.getScalarData("oxideFraction");

    activeCells.clear();
    isActiveCell.assign(cellSet.getNumberOfCells(), false);

    // Track which cells are "Core Active" (Oxide or Source)
    std::vector<bool> isCoreActive(cellSet.getNumberOfCells(), false);

#pragma omp parallel for
    for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
      if (isMaterial((*materials)[i], substrate) ||
          isMaterial((*materials)[i], oxide)) {
        bool active = (*oxideFractions)[i] > 1e-6;
        if (!active) {
          // Check if touching ambient (Source term)
          for (auto n : cellSet.getNeighbors(i)) {
            if (n >= 0 && isMaterial((*materials)[n], ambient)) {
              active = true;
              break;
            }
          }
        }
        isCoreActive[i] = active;
      }
    }

    // Expand to neighbors to allow front propagation
    for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
      if (isMaterial((*materials)[i], substrate) ||
          isMaterial((*materials)[i], oxide)) {
        // Solve if Core Active OR neighbor of Core Active (Narrow Band)
        bool shouldSolve = isCoreActive[i];
        if (!shouldSolve) {
          for (auto n : cellSet.getNeighbors(i)) {
            if (n >= 0 && isCoreActive[n]) {
              shouldSolve = true;
              break;
            }
          }
        }

        if (shouldSolve) {
          activeCells.push_back(i);
          isActiveCell[i] = true;
        }
      }
    }
  }

  // Calculates the reaction rate k based on the local E-field and oxide
  // fraction.
  T getReactionRate(T Ey, T oxideFrac) {
    const T k_base = params.get("reactionRateConstant");
    const T alpha = params.get("eFieldInfluence");

    // Reaction slows down as Silicon is consumed (oxideFrac approaches 1.0)
    T availableSi = std::max(0.0, 1.0 - oxideFrac);

    // Use absolute value of Ey to ensure positive rate enhancement
    return k_base * (1.0 + alpha * std::abs(Ey)) * availableSi;
  }

  // Calculates diffusion coefficient.
  T getDiffusivity(int material, T oxideFrac, T D_ox) {
    if (isMaterial(material, substrate) || isMaterial(material, oxide)) {
      return D_ox * oxideFrac;
    }
    return 0.0;
  }

  // Pre-calculates reaction rates for all cells
  void calculateReactionRates(std::vector<T> &reactionRates) {
    auto oxideFractions = cellSet.getScalarData("oxideFraction");
    auto materials = cellSet.getScalarData("Material");

    // Define CSV grid parameters for E-field lookup
    const T csvDx = 1.0; // Should match writeCSV.py
    const T csvExtent = 150.0;
    const int nx = static_cast<int>(csvExtent / csvDx);

#pragma omp parallel for
    for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
      if (isMaterial((*materials)[i], substrate) ||
          isMaterial((*materials)[i], oxide)) {

        // E-field lookup based on cell's coordinates
        auto center = cellSet.getCellCenter(i);
        T gen_x = center[0] + (csvExtent / 2.0);
        int ix = static_cast<int>(gen_x / csvDx);
        
        int idx = -1;
        if constexpr (D == 2) {
            idx = ix;
        } else {
            T gen_y = center[1] + (csvExtent / 2.0);
            int iy = static_cast<int>(gen_y / csvDx);
            if (iy >= 0 && iy < nx) idx = iy * nx + ix;
        }

        T Ey = 0.0;
        if (idx >= 0 && static_cast<size_t>(idx) < electricField1D.size()) {
          Ey = electricField1D[idx];
        }

        reactionRates[i] = getReactionRate(Ey, (*oxideFractions)[i]);
      } else {
        reactionRates[i] = 0.0;
      }
    }
  }

  // Solves diffusion using Explicit Finite Difference (Forward Euler)
  void solveDiffusionExplicit(T dt, const std::vector<T> &reactionRates) {
    auto oxidant = cellSet.getScalarData("oxidant");
    auto materials = cellSet.getScalarData("Material");
    auto oxideFractions = cellSet.getScalarData("oxideFraction");

    if (nextOxidant.size() != cellSet.getNumberOfCells()) {
      nextOxidant.resize(cellSet.getNumberOfCells());
    }

    const T ambientOxidant = params.get("ambientOxidant");
    const T D_ox = params.get("oxidantDiffusivity");
    const T dx = cellSet.getGridDelta();
    const T dtdx2 = dt / (dx * dx);

#pragma omp parallel for
    for (int idx = 0; idx < static_cast<int>(activeCells.size()); ++idx) {
      int i = activeCells[idx];
      T C_old = (*oxidant)[i];
      T k = reactionRates[i];
      T diffusion_term = 0.0;
      T D_center = getDiffusivity((*materials)[i], (*oxideFractions)[i], D_ox);

      const auto &neighbors = cellSet.getNeighbors(i);
      for (auto n : neighbors) {
        if (n < 0) continue;
        int mat_n = static_cast<int>((*materials)[n]);
        if (isMaterial(mat_n, mask)) continue;

        if (isMaterial(mat_n, ambient)) {
          // Dirichlet BC
          // Use D_ox for interface flux to allow oxidant entry into Silicon
          diffusion_term += D_ox * (ambientOxidant - C_old);
        } else if (isActiveCell[n]) {
          // Substrate/Oxide neighbor
          T D_neighbor = getDiffusivity(mat_n, (*oxideFractions)[n], D_ox);
          T D_eff = 0.5 * (D_center + D_neighbor);
          diffusion_term += D_eff * ((*oxidant)[n] - C_old);
        }
      }
      
      nextOxidant[i] = C_old + dtdx2 * diffusion_term - dt * k * C_old;
    }

    // Write back
#pragma omp parallel for
    for (int idx = 0; idx < static_cast<int>(activeCells.size()); ++idx) {
      int i = activeCells[idx];
      (*oxidant)[i] = nextOxidant[i];
    }
  }

  // Updates the oxide fraction based on the calculated oxidant concentration.
  void updateOxideFraction(T dt, const std::vector<T> &reactionRates) {
    auto oxidant = cellSet.getScalarData("oxidant");
    auto oxideFractions = cellSet.getScalarData("oxideFraction");
    auto materials = cellSet.getScalarData("Material");

#pragma omp parallel for
    for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
      if (isMaterial((*materials)[i], substrate) ||
          isMaterial((*materials)[i], oxide)) {
        T C = (*oxidant)[i];
        if (C > 1e-12) {
          // Recalculate rate k
          T k = reactionRates[i];

          // Explicit Euler update
          T dFrac = k * C * dt;

          (*oxideFractions)[i] += dFrac;

          // Clamp to 1.0 (100% Oxide)
          if ((*oxideFractions)[i] > 1.0)
            (*oxideFractions)[i] = 1.0;

          if ((*oxideFractions)[i] > 0.5) {
            (*materials)[i] = oxide;
          }
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
// MAIN
// -----------------------------------------------------------------------------

int main(int argc, char **argv) {
  cs::Logger::setLogLevel(cs::LogLevel::DEBUG);

  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <config file>" << std::endl;
    return 1;
  }

  OxidationSimulation simulation(argv[1]);
  simulation.apply();
}