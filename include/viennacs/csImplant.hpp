// apply the model to the geometry
#pragma once

#include <csDenseCellSet.hpp>
#include <csImplantModel.hpp>

#include <vcLogger.hpp>
#include <random>
#include <cmath>

namespace viennacs {

using namespace viennacore;

template <class NumericType, int D> class Implant {
  SmartPointer<DenseCellSet<NumericType, D>> cellSet_;
  SmartPointer<ImplantModel<NumericType, D>> model_;
  util::Parameters params_;
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

  void setParameters(util::Parameters passedParameters){
      params_ = passedParameters;
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
    //model_->getDepthProfile(1, params_);
    int totalNumberOfIons = 50000; // idk if this is the best way to determine the dose
    NumericType xLength = cellSet_->getBoundingBox()[1][1]; // the beam hits exactly in the middle, orthogonal to the y plane
    NumericType yLength = cellSet_->getBoundingBox()[1][0];
    auto gridDelta = cellSet_->getGridDelta();
    int numberOfVerticalCells = xLength / gridDelta;
    int numberOfHorizontalCells = yLength / gridDelta;
    //int halfNumberOfVerticalCells = halfXlength / gridDelta;
    auto concentration = cellSet_->getScalarData("concentration");


    double angle = -45;  // angle in degrees
    // Convert angle to radians (C++ trig functions use radians)
    double radians = angle * M_PI / 180.0;

    // iterate over all horizontal cells
    for (int i = 0; i < numberOfHorizontalCells * 2; i++){
        // iterate over all vertical cells
        for (int j = 0; j < numberOfVerticalCells * 2; j++){
            NumericType x = i * gridDelta * 0.5; // going in half steps to avoid weird stripes
            NumericType y = j * gridDelta * 0.5;
            NumericType shifted_y = y - 0.25 * yLength;
            NumericType depth = pow(std::cos(radians) * x, 2) + pow(std::sin(radians) * shifted_y, 2);
            NumericType lateralDisplacement = std::sin(radians) * x + std::cos(radians) * shifted_y;
            std::array<double, 3> coords{x, y, 0};
            auto index = cellSet_->getIndex(coords);
            (*concentration)[index] = model_->getDepthProfile(depth, params_) * model_->getLateralProfile(lateralDisplacement, depth, params_);
        }
    }
  }
};

} // namespace viennacs