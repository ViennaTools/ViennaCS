#pragma once

#include "csBVH.hpp"
#include "csTracePath.hpp"

#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

#include <vcUtil.hpp>
#include <vcVectorType.hpp>

#include <bitset>

namespace viennacs {

using namespace viennacore;

/**
  This class represents a cell-based voxel implementation of a volume. The
  depth of the cell set in z-direction can be specified.
*/
template <class T, int D> class DenseCellSet {
private:
  using gridType = SmartPointer<viennals::Mesh<T>>;
  using levelSetsType = std::vector<SmartPointer<viennals::Domain<T, D>>>;
  using materialMapType = SmartPointer<viennals::MaterialMap>;

  levelSetsType levelSets;
  gridType cellGrid = nullptr;
  SmartPointer<viennals::Domain<T, D>> surface = nullptr;
  SmartPointer<BVH<T, D>> bvh = nullptr;
  materialMapType materialMap = nullptr;

  T gridDelta;
  T depth = 0.;
  size_t numberOfCells = 0;
  int BVHlayers = 0;

  std::vector<std::array<int, 2 * D>> cellNeighbors; // -x, x, -y, y, -z, z
  viennahrle::Index<D> minIndex, maxIndex;

  bool cellSetAboveSurface = false;
  int coverMaterial = -1;
  std::bitset<D> periodicBoundary;

  std::vector<T> *fillingFractions_;
  const T eps = 1e-4;

public:
  DenseCellSet() = default;

  DenseCellSet(levelSetsType passedLevelSets,
               materialMapType passedMaterialMap = nullptr, T passedDepth = 0.,
               bool passedCellSetPosition = false)
      : levelSets(passedLevelSets), cellSetAboveSurface(passedCellSetPosition) {
    fromLevelSets(passedLevelSets, passedMaterialMap, passedDepth);
  }

  void fromLevelSets(levelSetsType passedLevelSets,
                     materialMapType passedMaterialMap = nullptr,
                     T passedDepth = 0.) {
    levelSets = passedLevelSets;
    materialMap = passedMaterialMap;

    if (cellGrid == nullptr)
      cellGrid = SmartPointer<viennals::Mesh<T>>::New();

    if (surface == nullptr)
      surface = SmartPointer<viennals::Domain<T, D>>::New(levelSets.back());
    else
      surface->deepCopy(levelSets.back());

    gridDelta = surface->getGrid().getGridDelta();

    depth = passedDepth;
    auto levelSetsInOrder = getLevelSetsInOrder();

    calculateMinMaxIndex(levelSetsInOrder);
    // viennals::ToVoxelMesh<T, D>(levelSetsInOrder, cellGrid).apply();
    // lsToVoxelMesh also saves the extent in the cell grid
    // save the extent of the resulting mesh
    for (unsigned i = 0; i < D; ++i) {
      cellGrid->minimumExtent[i] = std::numeric_limits<T>::max();
      cellGrid->maximumExtent[i] = std::numeric_limits<T>::lowest();
    }

    std::unordered_map<viennahrle::Index<D>, size_t,
                       typename viennahrle::Index<D>::hash>
        pointIdMapping;
    size_t currentPointId = 0;

    // prepare mesh for material ids and filling fractions
    cellGrid->cellData.insertNextScalarData(
        typename viennals::PointData<T>::ScalarDataType(), "Material");
    cellGrid->cellData.insertNextScalarData(
        typename viennals::PointData<T>::ScalarDataType(), "FillingFraction");
    auto &materialIds = *(cellGrid->cellData.getScalarData(0));
    auto &fillingFractions = *(cellGrid->cellData.getScalarData(1));
    const bool useMaterialMap = materialMap != nullptr;

    // set up iterators for all materials
    std::vector<viennahrle::ConstDenseCellIterator<
        typename viennals::Domain<T, D>::DomainType>>
        iterators;
    for (auto it = levelSetsInOrder.begin(); it != levelSetsInOrder.end();
         ++it) {
      iterators.emplace_back((*it)->getDomain(), minIndex);
    }

    // move iterator for lowest material id and then adjust others if they are
    // needed
    for (; iterators.front().getIndices() < maxIndex;
         iterators.front().next()) {
      // go over all materials
      for (unsigned materialId = 0; materialId < levelSetsInOrder.size();
           ++materialId) {

        auto &cellIt = iterators[materialId];

        cellIt.goToIndicesSequential(iterators.front().getIndices());

        // find out whether the centre of the box is inside
        double centerValue = 0.;
        for (int i = 0; i < (1 << D); ++i) {
          centerValue += cellIt.getCorner(i).getValue();
        }

        if (centerValue <= 0.) {
          std::array<unsigned, 1 << D> voxel;
          bool addVoxel = false;
          // now insert all points of voxel into pointList
          for (unsigned i = 0; i < (1 << D); ++i) {
            viennahrle::Index<D> index;
            addVoxel = true;
            for (unsigned j = 0; j < D; ++j) {
              index[j] =
                  cellIt.getIndices(j) + cellIt.getCorner(i).getOffset()[j];
              if (index[j] > maxIndex[j]) {
                addVoxel = false;
                break;
              }
            }
            if (addVoxel) {
              auto pointIdValue = std::make_pair(index, currentPointId);
              auto pointIdPair = pointIdMapping.insert(pointIdValue);
              voxel[i] = pointIdPair.first->second;
              if (pointIdPair.second) {
                ++currentPointId;
              }
            } else {
              break;
            }
          }

          // create element if inside domain bounds
          if (addVoxel) {
            int material = materialId;
            if (useMaterialMap)
              material = indexToMaterial(materialId, materialMap);

            if constexpr (D == 3) {
              // reorder elements for hexas to be ordered correctly
              std::array<unsigned, 8> hexa{voxel[0], voxel[1], voxel[3],
                                           voxel[2], voxel[4], voxel[5],
                                           voxel[7], voxel[6]};
              cellGrid->hexas.push_back(hexa);
            } else {
              std::array<unsigned, 4> tetra{voxel[0], voxel[2], voxel[3],
                                            voxel[1]};
              cellGrid->tetras.push_back(tetra);
            }
            materialIds.push_back(material);
            fillingFractions.push_back(std::max(1., centerValue));
          }
          // jump out of material for loop
          break;
        }
      }
    }

    // now insert points
    cellGrid->nodes.resize(pointIdMapping.size());
    for (auto it = pointIdMapping.begin(); it != pointIdMapping.end(); ++it) {
      std::array<T, 3> coords{};
      for (unsigned i = 0; i < D; ++i) {
        coords[i] = gridDelta * it->first[i];

        // save extent
        if (coords[i] < cellGrid->minimumExtent[i]) {
          cellGrid->minimumExtent[i] = coords[i];
        } else if (coords[i] > cellGrid->maximumExtent[i]) {
          cellGrid->maximumExtent[i] = coords[i];
        }
      }
      cellGrid->nodes[it->second] = coords;
    }

#ifndef NDEBUG
    int db_ls = 0;
    for (auto &ls : levelSetsInOrder) {
      auto mesh = SmartPointer<viennals::Mesh<T>>::New();
      viennals::ToSurfaceMesh<T, D>(ls, mesh).apply();
      viennals::VTKWriter<T>(mesh, "cellSet_debug_" + std::to_string(db_ls++) +
                                       ".vtp")
          .apply();
    }
    viennals::VTKWriter<T>(cellGrid, "cellSet_debug_init.vtu").apply();
#endif

    // adjustMaterialIds();
    fillingFractions_ =
        cellGrid->getCellData().getScalarData("FillingFraction");

    // create filling fractions as default scalar cell data
    numberOfCells = cellGrid->template getElements<(1 << D)>().size();

    // calculate number of BVH layers
    for (unsigned i = 0; i < D; ++i) {
      cellGrid->minimumExtent[i] -= eps;
      cellGrid->maximumExtent[i] += eps;
    }
    auto minExtent = cellGrid->maximumExtent[0] - cellGrid->minimumExtent[0];
    minExtent = std::min(minExtent, cellGrid->maximumExtent[1] -
                                        cellGrid->minimumExtent[1]);
    if constexpr (D == 3)
      minExtent = std::min(minExtent, cellGrid->maximumExtent[2] -
                                          cellGrid->minimumExtent[2]);

    BVHlayers = 0;
    while (minExtent / 2 > gridDelta) {
      BVHlayers++;
      minExtent /= 2;
    }

    bvh = SmartPointer<BVH<T, D>>::New(getBoundingBox(), BVHlayers);
    buildBVH();
  }

  std::array<VectorType<T, D>, 2> getBoundingBox() const {
    if constexpr (D == 3)
      return {cellGrid->minimumExtent, cellGrid->maximumExtent};
    else
      return {Vec2D<T>{cellGrid->minimumExtent[0], cellGrid->minimumExtent[1]},
              Vec2D<T>{cellGrid->maximumExtent[0], cellGrid->maximumExtent[1]}};
  }

  void setPeriodicBoundary(std::array<bool, D> isPeriodic) {
    for (int i = 0; i < D; i++) {
      periodicBoundary[i] = isPeriodic[i];
    }
  }

  std::vector<T> *addScalarData(std::string name, T initValue = 0.) {
    if (cellGrid->getCellData().getScalarData(name, false) != nullptr) {
      auto data = cellGrid->getCellData().getScalarData(name);
      data->resize(numberOfCells, initValue);
      std::fill(data->begin(), data->end(), initValue);
      return data;
    }
    std::vector<T> newData(numberOfCells, initValue);
    cellGrid->getCellData().insertNextScalarData(std::move(newData), name);
    fillingFractions_ =
        cellGrid->getCellData().getScalarData("FillingFraction");
    return cellGrid->getCellData().getScalarData(name);
  }

  T getDepth() const { return depth; }

  T getGridDelta() const { return gridDelta; }

  std::vector<Vec3D<T>> &getNodes() const { return cellGrid->getNodes(); }

  const Vec3D<T> &getNode(unsigned int idx) const {
    return cellGrid->getNodes()[idx];
  }

  std::vector<std::array<unsigned, (1 << D)>> &getElements() const {
    return cellGrid->template getElements<(1 << D)>();
  }

  const std::array<unsigned, (1 << D)> &getElement(unsigned int idx) const {
    return cellGrid->template getElements<(1 << D)>()[idx];
  }

  SmartPointer<viennals::Domain<T, D>> getSurface() { return surface; }

  SmartPointer<viennals::Mesh<T>> getCellGrid() { return cellGrid; }

  levelSetsType getLevelSets() const { return levelSets; }

  size_t getNumberOfCells() const { return numberOfCells; }

  std::vector<T> *getFillingFractions() const { return fillingFractions_; }

  T getFillingFraction(const std::array<T, D> &point) {
    Vec3D<T> point3 = {0., 0., 0.};
    for (int i = 0; i < D; i++)
      point3[i] = point[i];
    auto idx = findIndex(point3);
    if (idx < 0)
      return -1.;

    return getFillingFractions()->at(idx);
  }

  T getAverageFillingFraction(const std::array<T, 3> &point,
                              const T radius) const {
    T sum = 0.;
    int count = 0;
    for (int i = 0; i < numberOfCells; i++) {
      auto &cell = cellGrid->template getElements<(1 << D)>()[i];
      auto node = cellGrid->getNodes()[cell[0]];
      for (int j = 0; j < D; j++)
        node[j] += gridDelta / 2.;
      if (Distance(node, point) < radius) {
        sum += fillingFractions_->at(i);
        count++;
      }
    }
    return sum / count;
  }

  Vec3D<T> getCellCenter(unsigned long idx) const {
    auto center =
        cellGrid
            ->getNodes()[cellGrid->template getElements<(1 << D)>()[idx][0]];
    for (int i = 0; i < D; i++)
      center[i] += gridDelta / 2.;
    return center;
  }

  int getIndex(const Vec3D<T> &point) const { return findIndex(point); }

  std::vector<T> *getScalarData(std::string name) {
    return cellGrid->getCellData().getScalarData(name);
  }

  std::vector<std::string> getScalarDataLabels() const {
    std::vector<std::string> labels;
    auto numScalarData = cellGrid->getCellData().getScalarDataSize();
    for (int i = 0; i < numScalarData; i++) {
      labels.push_back(cellGrid->getCellData().getScalarDataLabel(i));
    }
    return labels;
  }

  // Set whether the cell set should be created below (false) or above (true)
  // the surface.
  void setCellSetPosition(const bool passedCellSetPosition) {
    cellSetAboveSurface = passedCellSetPosition;
  }

  void setCoverMaterial(const int passedCoverMaterial) {
    coverMaterial = passedCoverMaterial;
  }

  bool getCellSetPosition() const { return cellSetAboveSurface; }

  // Sets the filling fraction at given cell index.
  bool setFillingFraction(const int idx, const T fill) {
    if (idx < 0)
      return false;

    getFillingFractions()->at(idx) = fill;
    return true;
  }

  // Sets the filling fraction for cell which contains given point.
  bool setFillingFraction(const Vec3D<T> &point, const T fill) {
    auto idx = findIndex(point);
    return setFillingFraction(idx, fill);
  }

  // Add to the filling fraction at given cell index.
  bool addFillingFraction(const int idx, const T fill) {
    if (idx < 0)
      return false;

    fillingFractions_->at(idx) += fill;
    return true;
  }

  // Add to the filling fraction for cell which contains given point.
  bool addFillingFraction(const Vec3D<T> &point, T fill) {
    auto idx = findIndex(point);
    return addFillingFraction(idx, fill);
  }

  // Add to the filling fraction for cell which contains given point only if the
  // cell has the specified material ID.
  bool addFillingFractionInMaterial(const Vec3D<T> &point, T fill,
                                    int materialId) {
    auto idx = findIndex(point);
    if (getScalarData("Material")->at(idx) == materialId)
      return addFillingFraction(idx, fill);
    else
      return false;
  }

  // Write the cell set as .vtu file
  void writeVTU(const std::string &fileName) {
    viennals::VTKWriter<T>(cellGrid, fileName).apply();
  }

  // Save cell set data in simple text format
  void writeCellSetData(const std::string &fileName) const {
    auto numScalarData = cellGrid->getCellData().getScalarDataSize();

    std::ofstream file(fileName);
    file << numberOfCells << "\n";
    for (int i = 0; i < numScalarData; i++) {
      auto label = cellGrid->getCellData().getScalarDataLabel(i);
      file << label << ",";
    }
    file << "\n";

    for (size_t j = 0; j < numberOfCells; j++) {
      for (int i = 0; i < numScalarData; i++) {
        file << cellGrid->getCellData().getScalarData(i)->at(j) << ",";
      }
      file << "\n";
    }

    file.close();
  }

  // Read cell set data from text
  void readCellSetData(const std::string &fileName) {
    std::ifstream file(fileName);
    std::string line;

    if (!file.is_open()) {
      Logger::getInstance()
          .addWarning("Could not open file " + fileName)
          .print();
      return;
    }

    std::getline(file, line);
    if (std::stoi(line) != numberOfCells) {
      Logger::getInstance().addWarning("Incompatible cell set data.").print();
      return;
    }

    std::vector<std::string> labels;
    std::getline(file, line);
    {
      std::stringstream ss(line);
      std::string label;
      while (std::getline(ss, label, ',')) {
        labels.push_back(label);
      }
    }

    std::vector<std::vector<T> *> cellDataP;
    for (auto &label : labels) {
      auto dataP = getScalarData(label);
      if (dataP == nullptr) {
        addScalarData(label, 0.);
      }
    }

    for (auto &label : labels) {
      cellDataP.push_back(getScalarData(label));
    }

    std::size_t j = 0;
    while (std::getline(file, line)) {
      std::stringstream ss(line);
      std::size_t i = 0;
      std::string value;
      while (std::getline(ss, value, ','))
        cellDataP[i++]->at(j) = std::stod(value);

      j++;
    }
    assert(j == numberOfCells && "Data incompatible");

    file.close();
  }

  // Clear the filling fractions
  void clear() {
    auto ff = getFillingFractions();
    std::fill(ff->begin(), ff->end(), 0.);
  }

  // Update the material IDs of the cell set. This function should be called if
  // the level sets, the cell set is made out of, have changed. This does not
  // work if the surface of the volume has changed. In this case, call the
  // function "updateSurface" first.
  void updateMaterials() {
    auto materialIds = getScalarData("Material");

    // create overlay material
    auto levelSetsInOrder = getLevelSetsInOrder();

    // set up iterators for all materials
    std::vector<viennahrle::ConstDenseCellIterator<
        typename viennals::Domain<T, D>::DomainType>>
        iterators;
    for (const auto &ls : levelSetsInOrder) {
      iterators.emplace_back(ls->getDomain(), minIndex);
    }

    // move iterator for lowest material id and then adjust others if they are
    // needed
    const materialMapType matMapPtr = materialMap;
    unsigned cellIdx = 0;
    for (; iterators.front().getIndices() < maxIndex;
         iterators.front().next()) {
      // go over all materials
      for (unsigned materialId = 0; materialId < levelSetsInOrder.size();
           ++materialId) {

        auto &cellIt = iterators[materialId];
        cellIt.goToIndicesSequential(iterators.front().getIndices());

        // find out whether the centre of the box is inside
        T centerValue = 0.;
        for (int i = 0; i < (1 << D); ++i) {
          centerValue += cellIt.getCorner(i).getValue();
        }

        if (centerValue <= 0.) {
          bool isVoxel = true;
          // check if voxel is in bounds
          for (unsigned i = 0; i < (1 << D) && isVoxel; ++i) {
            viennahrle::Index<D> index;
            for (unsigned j = 0; j < D; ++j) {
              index[j] =
                  cellIt.getIndices(j) + cellIt.getCorner(i).getOffset()[j];
              if (index[j] > maxIndex[j]) {
                isVoxel = false;
                break;
              }
            }
          }

          if (isVoxel) {
            if (matMapPtr) {
              materialIds->at(cellIdx++) =
                  indexToMaterial(materialId, matMapPtr);
            } else {
              materialIds->at(cellIdx++) = materialId;
            }
          }

          // jump out of material for loop
          break;
        }
      }
    }
    assert(cellIdx == numberOfCells &&
           "Cell set changed in `updateMaterials()'");
  }

  // Updates the surface of the cell set. The new surface should be below the
  // old surface as this function can only remove cells from the cell set.
  void updateSurface() {
    auto updateCellGrid = SmartPointer<viennals::Mesh<T>>::New();

    viennals::ToVoxelMesh<T, D> voxelConverter(updateCellGrid);
    {
      auto plane =
          SmartPointer<viennals::Domain<T, D>>::New(surface->getGrid());
      T origin[D] = {0.};
      T normal[D] = {0.};
      origin[D - 1] = depth;
      normal[D - 1] = 1.;

      viennals::MakeGeometry<T, D>(
          plane, SmartPointer<viennals::Plane<T, D>>::New(origin, normal))
          .apply();
      voxelConverter.insertNextLevelSet(plane);
    }
    voxelConverter.insertNextLevelSet(levelSets.back());
    voxelConverter.insertNextLevelSet(surface);
    voxelConverter.apply();

    auto cutMatIds = updateCellGrid->getCellData().getScalarData("Material");
    auto &elements = cellGrid->template getElements<(1 << D)>();

    const auto nCutCells =
        updateCellGrid->template getElements<(1 << D)>().size();

    auto numScalarData = cellGrid->getCellData().getScalarDataSize();

    for (int elIdx = nCutCells - 1; elIdx >= 0; --elIdx) {
      if (cutMatIds->at(elIdx) == 2) {
        for (int i = 0; i < numScalarData; i++) {
          auto data = cellGrid->getCellData().getScalarData(i);
          data->erase(data->begin() + elIdx);
        }
        elements.erase(elements.begin() + elIdx);
      }
    }
    numberOfCells = elements.size();
    surface->deepCopy(levelSets.back());

    buildBVH();
  }

  // Merge a trace path to the cell set.
  void mergePath(TracePath<T> &path, T factor = 1.) {
    auto ff = getFillingFractions();
    if (!path.getData().empty()) {
      for (const auto it : path.getData()) {
        ff->at(it.first) += it.second / factor;
      }
    }

    if (!path.getGridData().empty()) {
      const auto &data = path.getGridData();
      for (size_t idx = 0; idx < numberOfCells; idx++) {
        ff->at(idx) += data[idx] / factor;
      }
    }
  }

  void buildNeighborhood(bool forceRebuild = false) {
    if (!cellNeighbors.empty() && !forceRebuild)
      return;

    Timer timer;
    timer.start();
    const auto &cells = cellGrid->template getElements<(1 << D)>();
    const auto &nodes = cellGrid->getNodes();
    unsigned const numNodes = nodes.size();
    unsigned const numCells = cells.size();
    cellNeighbors.resize(numCells);
    const bool usePeriodicBoundary = periodicBoundary.any();

    std::vector<std::vector<unsigned>> nodeCellConnections(numNodes);

    // for each node, store which cells are connected with the node
    for (unsigned cellIdx = 0; cellIdx < numCells; cellIdx++) {
      for (unsigned cellNodeIdx = 0; cellNodeIdx < (1 << D); cellNodeIdx++) {
        nodeCellConnections[cells[cellIdx][cellNodeIdx]].push_back(cellIdx);
      }
    }

#pragma omp parallel for
    for (int cellIdx = 0; cellIdx < numCells; cellIdx++) {
      auto coord = nodes[cells[cellIdx][0]];
      for (int i = 0; i < D; i++) {
        coord[i] += gridDelta / 2.;
        cellNeighbors[cellIdx][i] = -1;
        cellNeighbors[cellIdx][i + D] = -1;
      }

      if (usePeriodicBoundary) {
        auto onBoundary = isBoundaryCell(coord);
        if (onBoundary.any()) {
          /*look for neighbor cells using BVH*/
          for (std::size_t i = 0; i < 2 * D; i++) {
            auto neighborCoord = coord;
            bool minBoundary = i % 2 == 0;

            if (onBoundary.test(i)) {
              // wrap around boundary
              if (!minBoundary) {
                neighborCoord[i / 2] =
                    cellGrid->minimumExtent[i / 2] + gridDelta / 2.;
              } else {
                neighborCoord[i / 2] =
                    cellGrid->maximumExtent[i / 2] - gridDelta / 2.;
              }
            } else {
              neighborCoord[i / 2] += minBoundary ? -gridDelta : gridDelta;
            }
            cellNeighbors[cellIdx][i] = findIndex(neighborCoord);
          }
          continue;
        }
      }

      for (unsigned cellNodeIdx = 0; cellNodeIdx < (1 << D); cellNodeIdx++) {
        auto &cellsAtNode = nodeCellConnections[cells[cellIdx][cellNodeIdx]];

        for (const auto &neighborCell : cellsAtNode) {
          if (neighborCell != cellIdx) {

            auto neighborCoord = getCellCenter(neighborCell);

            if (Distance(coord, neighborCoord) < gridDelta + eps) {

              for (int i = 0; i < D; i++) {
                if (coord[i] - neighborCoord[i] > gridDelta / 2.) {
                  cellNeighbors[cellIdx][i * 2] = neighborCell;
                } else if (coord[i] - neighborCoord[i] < -gridDelta / 2.) {
                  cellNeighbors[cellIdx][i * 2 + 1] = neighborCell;
                }
              }
            }
          }
        }
      }
    }
    timer.finish();
    Logger::getInstance()
        .addTiming("Building cell set neighborhood structure took",
                   timer.currentDuration * 1e-9)
        .print();
  }

  const std::array<int, 2 * D> &getNeighbors(unsigned long cellIdx) const {
    assert(cellIdx < numberOfCells && "Cell idx out of bounds");
    return cellNeighbors[cellIdx];
  }

  bool isPointInCell(const Vec3D<T> &point, unsigned int cellIdx) const {
    const auto &elem = getElement(cellIdx);
    const auto &cellMin = getNode(elem[0]);
    return isPointInCell(point, cellMin);
  }

private:
  int findIndex(const Vec3D<T> &point) const {
    const auto &elems = cellGrid->template getElements<(1 << D)>();
    const auto &nodes = cellGrid->getNodes();
    int idx = -1;

    auto cellIds = bvh->getCellIds(point);
    if (!cellIds)
      return idx;
    for (const auto cellId : *cellIds) {
      if (isPointInCell(point, nodes[elems[cellId][0]])) {
        idx = cellId;
        break;
      }
    }
    return idx;
  }

  void adjustMaterialIds() {
    // This assumes the current material IDs correspond to the level-set index
    auto matIds = getScalarData("Material");
    if (!materialMap)
      return;

#pragma omp parallel for
    for (int i = 0; i < matIds->size(); i++) {
      matIds->at(i) =
          indexToMaterial(static_cast<int>(matIds->at(i)), materialMap);
    }
  }

  int indexToMaterial(int index, const materialMapType &materialMap) {
    // This takes into account the added coverMaterial layer
    if (!cellSetAboveSurface)
      --index;
    if (index >= 0 && index < materialMap->getNumberOfLayers())
      return materialMap->getMaterialId(index);
    return coverMaterial;
  }

  auto getLevelSetsInOrder() {
    std::vector<SmartPointer<viennals::Domain<T, D>>> levelSetsInOrder;
    auto plane = SmartPointer<viennals::Domain<T, D>>::New(surface->getGrid());
    {
      T origin[D] = {0.};
      T normal[D] = {0.};
      origin[D - 1] = depth;
      normal[D - 1] = 1.;
      viennals::MakeGeometry<T, D>(
          plane, SmartPointer<viennals::Plane<T, D>>::New(origin, normal))
          .apply();
    }
    if (!cellSetAboveSurface)
      levelSetsInOrder.push_back(plane);
    for (const auto &ls : levelSets)
      levelSetsInOrder.push_back(ls);
    if (cellSetAboveSurface)
      levelSetsInOrder.push_back(plane);

    return levelSetsInOrder;
  }

  int findSurfaceHitPoint(Vec3D<T> &hitPoint, const Vec3D<T> &direction) {
    // find surface hit point
    auto idx = findIndex(hitPoint);

    if (idx > 0)
      return idx;

    auto moveDirection = multNew(direction, gridDelta / 2.);
    size_t sanityCounter = 0;
    while (idx < 0) {
      add(hitPoint, moveDirection);
      if (++sanityCounter > 100 || !checkBoundsPeriodic(hitPoint)) {
        return -1;
      }
      idx = findIndex(hitPoint);
    }

    return idx;
  }

  bool isPointInCell(const Vec3D<T> &point, const Vec3D<T> &cellMin) const {
    if constexpr (D == 3)
      return point[0] >= cellMin[0] && point[0] <= (cellMin[0] + gridDelta) &&
             point[1] >= cellMin[1] && point[1] <= (cellMin[1] + gridDelta) &&
             point[2] >= cellMin[2] && point[2] <= (cellMin[2] + gridDelta);
    else
      return point[0] >= cellMin[0] && point[0] <= (cellMin[0] + gridDelta) &&
             point[1] >= cellMin[1] && point[1] <= (cellMin[1] + gridDelta);
  }

  void buildBVH() {
    Timer timer;
    timer.start();
    auto &elems = cellGrid->template getElements<(1 << D)>();
    auto &nodes = cellGrid->getNodes();
    bvh->clearCellIds();

    for (size_t elemIdx = 0; elemIdx < elems.size(); elemIdx++) {
      for (size_t n = 0; n < (1 << D); n++) {
        auto &node = nodes[elems[elemIdx][n]];
        auto cell = bvh->getCellIds(node);
        if (cell == nullptr) {
          Logger::getInstance().addError("BVH building error.").print();
        }
        cell->insert(elemIdx);
      }
    }
    timer.finish();
    Logger::getInstance()
        .addTiming("Building cell set BVH took", timer.currentDuration * 1e-9)
        .print();
  }

  void
  calculateMinMaxIndex(const std::vector<SmartPointer<viennals::Domain<T, D>>>
                           &levelSetsInOrder) {
    // set to zero
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::numeric_limits<viennahrle::IndexType>::max();
      maxIndex[i] = std::numeric_limits<viennahrle::IndexType>::lowest();
    }
    for (unsigned l = 0; l < levelSetsInOrder.size(); ++l) {
      auto &grid = levelSetsInOrder[l]->getGrid();
      auto &domain = levelSetsInOrder[l]->getDomain();
      for (unsigned i = 0; i < D; ++i) {
        minIndex[i] = std::min(minIndex[i], (grid.isNegBoundaryInfinite(i))
                                                ? domain.getMinRunBreak(i)
                                                : grid.getMinBounds(i));

        maxIndex[i] = std::max(maxIndex[i], (grid.isPosBoundaryInfinite(i))
                                                ? domain.getMaxRunBreak(i)
                                                : grid.getMaxBounds(i));
      }
    }
  }

  std::bitset<2 * D> isBoundaryCell(const Vec3D<T> &cellCoord) {
    std::bitset<2 * D> onBoundary;
    for (int i = 0; i < 2 * D; i += 2) {
      if (!periodicBoundary[i / 2])
        continue;
      if (cellCoord[i / 2] - cellGrid->minimumExtent[i / 2] < gridDelta) {
        /* cell is at min boundary */
        onBoundary.set(i);
      }
      if (cellGrid->maximumExtent[i / 2] - cellCoord[i / 2] < gridDelta) {
        /* cell is at max boundary */
        onBoundary.set(i + 1);
      }
    }
    return onBoundary;
  }
};

} // namespace viennacs
