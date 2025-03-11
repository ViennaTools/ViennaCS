#include <csDenseCellSet.hpp>
#include <lsBooleanOperation.hpp>
#include <lsWriteVisualizationMesh.hpp>

template <typename T> struct Parameters {
  // Domain
  T gridDelta = 2.; // nm
  T xExtent = 80.0; // nm
  T yExtent = 80.0; // nm

  // Geometry
  T substrateHeight = 50.; // nm
  T coverHeight = 30.;
  T maskHeight = 10.;
  T holeRadius = xExtent / 6.;

  // Process
  T duration = 20.;
  T diffusionCoefficient = 1.; // nm²/s
  T velocity = 0.;             // Advection
  T timeStabilityFactor = 0.95;

  int substrateMaterial = 0;
  int maskMaterial = 1;
  int coverMaterial = 2;
};

using T = double;
constexpr int D = 2;

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

auto makeStructure(const Parameters<T> &params, materialMapType matMap) {
  const T gridDelta = params.gridDelta;
  const T substrateHeight = params.substrateHeight;
  const T coverHeight = params.coverHeight;
  const T maskHeight = params.maskHeight;
  const T holeRadius = params.holeRadius;
  ls::BoundaryConditionEnum boundaryConds[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY};
  boundaryConds[D - 1] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;
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
    /*origin[D - 1] = substrateHeight;
    normal[D - 1] = -1;
    makePlane(mask, origin, normal);
    normal[D - 1] = 1.;

    ls::BooleanOperation<T, D>(mask, maskAdd,
                             ls::BooleanOperationEnum::INTERSECT)
        .apply();*/

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

template <typename T> bool isMaterial(T x, int material) {
  return static_cast<int>(x) == material;
}

void addConcentration(cs::DenseCellSet<T, D> &cellSet,
                      const Parameters<T> &params) {
  // Add quantity to be diffused (on top of the cell material)
  auto concentration = cellSet.addScalarData("dopant", 0.);
  auto materials = cellSet.getScalarData("Material");

  // Boundary condition: constant concentration at the top (outside)
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isMaterial((*materials)[i], params.coverMaterial) &&
        cellSet.getCellCenter(i)[D - 1] <
            params.substrateHeight + params.gridDelta) {
      (*concentration)[i] = 1.;
    }
  }
}

// Explicit diffusion time-step (forward Euler)
void solveDiffusionStep(cs::DenseCellSet<T, D> &cellSet,
                        const Parameters<T> &params, T dt) {
  auto data = cellSet.getScalarData("dopant");
  auto materials = cellSet.getScalarData("Material");
  std::vector<T> solution(data->size(), 0.);
  const T C =
      dt * params.diffusionCoefficient / (params.gridDelta * params.gridDelta);

  // #pragma omp parallel for
  for (int e = 0; e < data->size(); e++) {
    if (!isMaterial((*materials)[e], params.substrateMaterial)) {
      solution[e] = (*data)[e];
      continue;
    }

    auto coord = cellSet.getCellCenter(e);
    int numNeighbors = 0;

    auto cellNeighbors = cellSet.getNeighbors(e);
    for (auto n : cellNeighbors) {
      if (n == -1 || isMaterial((*materials)[n], params.maskMaterial))
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
  for (auto ls : levelSets)
    writer.insertNextLevelSet(ls);
  writer.setMaterialMap(matMap);
  writer.apply();
}

int main(int argc, char **argv) {
  omp_set_num_threads(4);
  cs::Logger::setLogLevel(cs::LogLevel::INTERMEDIATE);

  Parameters<T> params;

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

  if (params.velocity != 0.) {
    const T stability = 2 * params.diffusionCoefficient / params.velocity;
    std::cout << "Stability: " << stability << std::endl;
    if (0.5 * stability <= params.gridDelta)
      std::cout << "Unstable parameters. Reduce grid spacing!" << std::endl;
  }

  auto materials = cellSet.getScalarData("Material");

  T duration = params.duration;
  T dt = std::min(params.gridDelta * params.gridDelta /
                      (params.diffusionCoefficient * 2 * D) *
                      params.timeStabilityFactor,
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
