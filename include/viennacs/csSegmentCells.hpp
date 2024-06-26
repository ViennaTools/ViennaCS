#pragma once

#include "csDenseCellSet.hpp"

namespace viennacs {

using namespace viennacore;

template <class NumericType, int D> class SegmentCells {
  SmartPointer<DenseCellSet<NumericType, D>> cellSet = nullptr;
  std::string cellTypeString = "CellType";
  int bulkMaterial = -1;

public:
  SegmentCells(const SmartPointer<DenseCellSet<NumericType, D>> &passedCellSet)
      : cellSet(passedCellSet) {}

  SegmentCells(const SmartPointer<DenseCellSet<NumericType, D>> &passedCellSet,
               const std::string cellTypeString, const int passedBulkMaterial)
      : cellSet(passedCellSet), cellTypeString(cellTypeString),
        bulkMaterial(passedBulkMaterial) {}

  void
  setCellSet(const SmartPointer<DenseCellSet<NumericType, D>> &passedCellSet) {
    cellSet = passedCellSet;
  }

  void setCellTypeString(const std::string passedCellTypeString) {
    cellTypeString = passedCellTypeString;
  }

  void setBulkMaterial(const int passedBulkMaterial) {
    bulkMaterial = passedBulkMaterial;
  }

  void apply() {
    auto cellType = cellSet->addScalarData(cellTypeString, -1.);
    cellSet->buildNeighborhood();
    auto materials = cellSet->getScalarData("Material");

#pragma omp parallel for
    for (int i = 0; i < materials->size(); ++i) {
      if (static_cast<int>(materials->at(i)) != bulkMaterial) {
        auto neighbors = cellSet->getNeighbors(i);
        for (auto n : neighbors) {
          if (n >= 0 && static_cast<int>(materials->at(n)) == bulkMaterial) {
            cellType->at(i) = 0.;
            break;
          }
        }
      } else {
        cellType->at(i) = 1.;
      }
    }
  }
};

} // namespace viennacs
