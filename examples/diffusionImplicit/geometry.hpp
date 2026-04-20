#pragma once

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include <vcUtil.hpp>

namespace geometry {

using namespace viennals;

template <class T, int D> using levelSetType = SmartPointer<Domain<T, D>>;
template <class T, int D> using levelSetsType = std::vector<levelSetType<T, D>>;
using materialMapType = SmartPointer<MaterialMap>;

template <class T, int D>
void addLevelSet(levelSetsType<T, D> &levelSets, levelSetType<T, D> levelSet,
                 materialMapType matMap, int material,
                 bool wrapLowerLevelSet = true) {
  if (!levelSets.empty() && wrapLowerLevelSet) {
    BooleanOperation<T, D>(levelSet, levelSets.back(),
                           BooleanOperationEnum::UNION)
        .apply();
  }

  levelSets.push_back(levelSet);
  matMap->insertNextMaterial(material);
}

template <class T, int D>
void makePlane(levelSetType<T, D> &domain, const T *origin, const T *normal) {
  MakeGeometry<T, D>(domain, SmartPointer<Plane<T, D>>::New(origin, normal))
      .apply();
}

template <class T, int D>
void makeBox(levelSetType<T, D> &domain, const T *minPoint, const T *maxPoint) {
  MakeGeometry<T, D>(domain, SmartPointer<Box<T, D>>::New(minPoint, maxPoint))
      .apply();
}

template <class T, int D>
auto makeStructure(util::Parameters &params, materialMapType matMap,
                   int substrateMaterial, int maskMaterial) {
  const T gridDelta = params.get("gridDelta");
  const T substrateHeight = params.get("substrateHeight");
  const T coverHeight = params.get("coverHeight");
  const T maskHeight = params.get("maskHeight");
  const T holeRadius = params.get("holeRadius");
  BoundaryConditionEnum boundaryConds[D] = {
      BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      BoundaryConditionEnum::REFLECTIVE_BOUNDARY};
  boundaryConds[D - 1] = BoundaryConditionEnum::INFINITE_BOUNDARY;
  T bounds[2 * D] = {-params.get("xExtent") / 2., params.get("xExtent") / 2.,
                     -params.get("yExtent") / 2., params.get("yExtent") / 2.};
  bounds[2 * D - 2] = 0.;
  bounds[2 * D - 1] = substrateHeight + maskHeight + gridDelta;

  T origin[D] = {};
  T normal[D] = {};
  normal[D - 1] = 1.;

  levelSetsType<T, D> levelSets;

  // Substrate
  origin[D - 1] = 0.;
  auto bottom = levelSetType<T, D>::New(bounds, boundaryConds, gridDelta);
  makePlane(bottom, origin, normal);
  addLevelSet(levelSets, bottom, matMap, substrateMaterial);

  origin[D - 1] = substrateHeight;
  auto substrate = levelSetType<T, D>::New(bounds, boundaryConds, gridDelta);
  makePlane(substrate, origin, normal);
  addLevelSet(levelSets, substrate, matMap, substrateMaterial);

  // Mask
  if (maskHeight > 0.) {
    auto mask = levelSetType<T, D>::New(bounds, boundaryConds, gridDelta);
    origin[D - 1] = substrateHeight + maskHeight;
    makePlane(mask, origin, normal);

    auto maskAdd = levelSetType<T, D>::New(bounds, boundaryConds, gridDelta);

    T minPoint[D] = {-holeRadius, -holeRadius};
    T maxPoint[D] = {holeRadius, holeRadius};
    minPoint[D - 1] = substrateHeight - gridDelta;
    maxPoint[D - 1] = substrateHeight + maskHeight + gridDelta;

    makeBox(maskAdd, minPoint, maxPoint);

    BooleanOperation<T, D>(mask, maskAdd,
                           BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    addLevelSet(levelSets, mask, matMap, maskMaterial);
  }

  return levelSets;
}

template <class T, int D>
void saveVolumeMesh(std::string name, levelSetsType<T, D> &levelSets,
                    materialMapType matMap) {
  WriteVisualizationMesh<T, D> writer;
  writer.setFileName(name);
  for (const auto &ls : levelSets)
    writer.insertNextLevelSet(ls);
  writer.setMaterialMap(matMap);
  writer.apply();
}

} // namespace geometry
