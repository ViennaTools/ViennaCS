// apply the model to the geometry
#pragma once

#include <csDenseCellSet.hpp>
#include <csImplantModel.hpp>

#include <vcLogger.hpp>
#include <random>

namespace viennacs {

using namespace viennacore;

template <class NumericType, int D> class Implant {
  SmartPointer<DenseCellSet<NumericType, D>> cellSet_;
  SmartPointer<ImplantModel<NumericType, D>> model_;
  std::vector<int> maskMaterials;

public:
  Implant() = default;

  void setCellSet(SmartPointer<DenseCellSet<NumericType, D>> passedCellSet) {
    cellSet_ = passedCellSet;
  }

  void setImplantModel(
      SmartPointer<ImplantModel<NumericType, D>> passedImplantModel) {
    model_ = passedImplantModel;
  }

  template <class... Mats> void setMaskMaterials(Mats... mats) {
    maskMaterials = {mats...};
  }

  void apply() {
    if (!model_) {
      Logger::getInstance()
          .addWarning("No implant model passed to Implant.")
          .print();
      return;
    }
    if (!cellSet_) {
      Logger::getInstance().addWarning("No cellSet passed to Implant.").print();
      return;
    }

    // apply the implant model to the cellSet_
    std::cout << "applying the model bli bla blub" << std::endl;
    model_->getDepthProfile(1);
    int totalNumberOfIons = 50000; // idk if this is the best way to determine the dose
    auto halfXlength = cellSet_->getBoundingBox()[1][1] * 1 / 2; // the beam hits exactly in the middle, orthogonal to the y plane
    auto yLength = cellSet_->getBoundingBox()[1][0];
    auto gridDelta = cellSet_->getGridDelta();
    int numberOfHorizontalCells = yLength / gridDelta;
    int halfNumberOfVerticalCells = halfXlength / gridDelta;
    auto concentration = cellSet_->getScalarData("aaa");
    // iterate over all horizontal cells
    for (int i = 0; i < numberOfHorizontalCells; i++){
        // iterate over vertical cells
        for (int j = -halfNumberOfVerticalCells; j < halfNumberOfVerticalCells; j++){
            NumericType depth = i * gridDelta;
            NumericType lateralDisplacement = j * gridDelta;
            std::array<double, 3> coords{depth, halfXlength + lateralDisplacement, 0};
            auto index = cellSet_->getIndex(coords);
            (*concentration)[index] = model_->getDepthProfile(depth) * model_->getLateralProfile(lateralDisplacement, depth);
        }
    }
  }
};

} // namespace viennacs