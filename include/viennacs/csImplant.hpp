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
    int totalNumberOfIons = 50000; // idk if this is the best way to determine the dose
    auto xCoord = cellSet_->getBoundingBox()[1][1] * 1/2; // the beam hits exactly in the middle, orthogonal to the y plane
    auto length = cellSet_->getBoundingBox()[1][0];
    auto gridDelta = cellSet_->getGridDelta();
    int numberOfHorizontalCells = length/gridDelta;
    std::cout << numberOfHorizontalCells << std::endl;

    // doing some normal distribution
    std::random_device rd{};
    std::mt19937 gen{rd()};
    NumericType mu = 10. / gridDelta;
    NumericType sigma = 2. / gridDelta;
    NumericType mu_lateral = 0. / gridDelta;
    NumericType sigma_lateral = 0.05/ gridDelta;
    std::normal_distribution d{mu, sigma};
    std::normal_distribution d_lateral{mu_lateral, sigma_lateral};
    auto xPosIon = [&d, &gen]{ return std::round(d(gen));};
    auto lateralDisplacementIon = [&d_lateral, &gen]{return std::round(d_lateral(gen));};
    auto concentration = cellSet_->getScalarData("aaa");
    for (int i = 0; i<totalNumberOfIons; i++){
        std::array<double, 3> coords{0, xCoord, 0};
        coords[0] = xPosIon() * gridDelta;
        coords[1] += lateralDisplacementIon() * gridDelta;
        auto index = cellSet_->getIndex(coords);
        (*concentration)[index]++;
    }
  }
};

} // namespace viennacs