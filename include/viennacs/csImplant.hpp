#pragma once

#include <csDenseCellSet.hpp>
#include <csImplantModel.hpp>

#include <vcLogger.hpp>

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
  }
};

} // namespace viennacs