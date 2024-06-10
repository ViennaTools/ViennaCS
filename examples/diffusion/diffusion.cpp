#include <csDenseCellSet.hpp>
#include <lsBooleanOperation.hpp>
#include <lsWriteVisualizationMesh.hpp>

namespace cs = viennacs;
namespace ls = viennals;

using T = double;
constexpr int D = 2;

using cs::util::Parameters;

const int substrateMaterial = 0;
const int maskMaterial = 1;
const int coverMaterial = 2;

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
  const T gridDelta = params.get("gridDelta");
  const T substrateHeight = params.get("substrateHeight");
  const T coverHeight = params.get("coverHeight");
  const T maskHeight = params.get("maskHeight");
  const T holeRadius = params.get("holeRadius");
  ls::BoundaryConditionEnum<D> boundaryConds[D] = {
      ls::BoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY};
  boundaryConds[D - 1] = ls::BoundaryConditionEnum<D>::INFINITE_BOUNDARY;
  T bounds[2 * D] = {-params.get("xExtent") / 2., params.get("xExtent") / 2.,
                     -params.get("yExtent") / 2., params.get("yExtent") / 2.};
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
  addLevelSet(levelSets, bottom, matMap, substrateMaterial);

  origin[D - 1] = substrateHeight;
  auto substrate = levelSetType::New(bounds, boundaryConds, gridDelta);
  makePlane(substrate, origin, normal);
  addLevelSet(levelSets, substrate, matMap, substrateMaterial);

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

    addLevelSet(levelSets, mask, matMap, maskMaterial);
  }

  return levelSets;
}

template <typename Material> bool isMaterial(Material x, int material) {
  return static_cast<int>(x) == material;
}

void addConcentration(cs::DenseCellSet<T, D> &cellSet,
                      const Parameters &params) {
  // Add quantity to be diffused (on top of the cell material)
  auto concentration = cellSet.addScalarData("dopant", 0.);
  auto materials = cellSet.getScalarData("Material");

  // Boundary condition: constant concentration at the top (outside)
#pragma omp parallel for
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isMaterial((*materials)[i], coverMaterial) &&
        cellSet.getCellCenter(i)[D - 1] <
            params.get("substrateHeight") + params.get("gridDelta")) {
      (*concentration)[i] = params.get("boundaryValue");
    }
  }
}

// Explicit diffusion time-step (forward Euler)
void solveDiffusionStep(cs::DenseCellSet<T, D> &cellSet,
                        const Parameters &params, T dt) {
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
  cs::Logger::setLogLevel(cs::LogLevel::INTERMEDIATE);

  Parameters params;
  if (argc > 1) {
    params.readConfigFile(argv[1]);
  } else {
    std::cout << "Usage: " << argv[0] << " <config file>" << std::endl;
    return 1;
  }
  omp_set_num_threads(params.get<int>("numThreads"));

  auto matMap = materialMapType::New();
  auto levelSets = makeStructure(params, matMap);

  cs::DenseCellSet<T, D> cellSet;
  T depth = params.get("substrateHeight") + params.get("coverHeight") + 10.;
  cellSet.setCellSetPosition(true); // isAboveSurface
  cellSet.setCoverMaterial(coverMaterial);
  cellSet.fromLevelSets(levelSets, matMap, depth);

  addConcentration(cellSet, params);
  // We need neighborhood information for solving the diffusion equation
  cellSet.buildNeighborhood();
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
