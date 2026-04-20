// apply the model to the geometry
#pragma once

#include "csImplantModel.hpp"
#include <csDenseCellSet.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <unordered_set>
#include <vcLogger.hpp>

namespace viennacs {

using namespace viennacore;

enum class ImplantDoseControl {
  Off,
  WaferDose,
  BeamDose,
};

template <class NumericType, int D> class Implant {
  SmartPointer<DenseCellSet<NumericType, D>> cellSet_;
  SmartPointer<ImplantModel<NumericType, D>> model_;
  std::vector<int> maskMaterials;
  NumericType angle_ = NumericType(0);
  NumericType dosePerCm2_ = NumericType(0);
  NumericType lengthUnitInCm_ = NumericType(1e-7);
  ImplantDoseControl doseControl_ = ImplantDoseControl::Off;
  bool writeBeamHits_ = false;
  bool outputConcentrationInCm3_ = false;
  std::string concentrationLabel_ = "concentration";
  std::string beamHitsLabel_ = "beamHits";

public:
  Implant() = default;

  void setCellSet(SmartPointer<DenseCellSet<NumericType, D>> passedCellSet) {
    cellSet_ = passedCellSet;
  }

  void setImplantAngle(const NumericType angle) { angle_ = angle; }

  void setDose(const NumericType dosePerCm2) { dosePerCm2_ = dosePerCm2; }

  void setLengthUnitInCm(const NumericType lengthUnitInCm) {
    lengthUnitInCm_ = lengthUnitInCm;
  }

  void setDoseControl(const ImplantDoseControl doseControl) {
    doseControl_ = doseControl;
  }

  void enableBeamHits(const bool enable = true) { writeBeamHits_ = enable; }

  void setConcentrationLabel(const std::string &label) {
    concentrationLabel_ = label;
  }

  void setBeamHitsLabel(const std::string &label) { beamHitsLabel_ = label; }

  void setOutputConcentrationInCm3(const bool enable = true) {
    outputConcentrationInCm3_ = enable;
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

    // now we apply the implant model to the cell set
    auto boundingBox = cellSet_->getBoundingBox();
    auto gridDelta = cellSet_->getGridDelta();
    auto concentration =
        cellSet_->getCellGrid()->getCellData().getScalarData(concentrationLabel_, true);
    if (!concentration)
      concentration = cellSet_->addScalarData(concentrationLabel_, 0.);
    auto beamHits = writeBeamHits_ ? cellSet_->addScalarData(beamHitsLabel_, 0.)
                                   : nullptr;
    auto material = cellSet_->getScalarData("Material");

    const NumericType minX = boundingBox[0][0];
    const NumericType maxX = boundingBox[1][0];
    const NumericType minY = boundingBox[0][1];
    const NumericType maxY = boundingBox[1][1];
    const NumericType minZ = boundingBox[0][2];
    const NumericType maxZ = boundingBox[1][2];

    NumericType xLength = std::abs(maxX - minX);
    NumericType yLength = std::abs(maxY - minY);
    NumericType zLength = std::abs(maxZ - minZ);

    int numberOfcellsXdirection = xLength / gridDelta;
    int numberOfcellsYdirection = yLength / gridDelta;
    int numberOfcellsZdirection = zLength / gridDelta;

    const auto relevant_Cells =
        std::max(1, static_cast<int>(std::ceil(model_->getMaxDepth() /
                                               std::max(gridDelta, NumericType(1e-9)))));
    const auto relevant_lateral_Cells =
        std::max(1, static_cast<int>(std::ceil(model_->getMaxLateralRange() /
                                               std::max(gridDelta, NumericType(1e-9)))));
    const std::unordered_set<int> maskMaterialSet(maskMaterials.begin(),
                                                  maskMaterials.end());

    double radians = angle_ * M_PI / 180;
    NumericType doseWeight = NumericType(1);
    if (doseControl_ != ImplantDoseControl::Off && dosePerCm2_ > NumericType(0)) {
      auto effectiveDosePerCm2 = dosePerCm2_;
      if (doseControl_ == ImplantDoseControl::BeamDose) {
        effectiveDosePerCm2 *= std::max(NumericType(0), std::cos(radians));
      }
      const auto dosePerLengthUnitsSquared =
          effectiveDosePerCm2 * lengthUnitInCm_ * lengthUnitInCm_;
      NumericType beamMeasure = gridDelta;
      if constexpr (D == 3)
        beamMeasure *= gridDelta;
      doseWeight = dosePerLengthUnitsSquared * beamMeasure;
    }

    if constexpr (D == 3) {
      // iterate over all the beams that hit the xy-plane from the z direction:
      for (int i = 0; i < numberOfcellsXdirection; i++) {
        for (int h = 0; h < numberOfcellsYdirection; h++) {
          NumericType initialX = minX + i * gridDelta + gridDelta;
          NumericType initialY = minY + h * gridDelta + gridDelta;
          NumericType initialZ = maxZ - gridDelta;
          std::array<NumericType, 3> initialCoords{initialX, initialY,
                                                   initialZ};
          int initialIndex;
          numberOfcellsZdirection = zLength / gridDelta;
          while (true) {
            initialCoords[0] = initialX;
            initialCoords[1] = initialY;
            initialCoords[2] = initialZ;
            initialIndex = cellSet_->getIndex(initialCoords);
            if (initialIndex == -1) {
              initialZ -= gridDelta * std::cos(radians);
              initialX -= gridDelta * std::sin(radians);
              numberOfcellsZdirection -= 1;
              if (initialZ < minZ - gridDelta || initialX < minX - xLength ||
                  initialX > maxX + xLength) {
                break;
              }
              continue;
            } else {
              const auto materialId = static_cast<int>((*material)[initialIndex]);
              if (materialId == 0) {
                initialZ -= gridDelta * std::cos(radians);
                initialX -= gridDelta * std::sin(radians);
                numberOfcellsZdirection -= 1;
                if (initialZ < minZ - gridDelta || initialX < minX - xLength ||
                    initialX > maxX + xLength) {
                  break;
                }
                continue;
              }
              if (maskMaterialSet.count(materialId) > 0) {
                break;
              }
              // iterate over all the cells in y [z] direction (depth), and
              // then iterate over all the cells in x [&y] direction to get
              // the lateral displacement
              for (int j = 0; j < relevant_Cells + 1; j++) {
                for (int k = 0; k < relevant_lateral_Cells + 1; k++) {
                  for (int l = 0; l < relevant_lateral_Cells + 1; l++) {
                    NumericType zCord = j * gridDelta;
                    NumericType xCord = initialX + k * gridDelta;
                    NumericType yCord = initialY + l * gridDelta;
                    NumericType shifted_xCord =
                        xCord - (relevant_lateral_Cells * gridDelta) / 2;
                    NumericType shifted_yCord =
                        yCord - (relevant_lateral_Cells * gridDelta) / 2;
                    NumericType shifted_zCord = initialZ - zCord;
                    // std::cout << "coord: [" << shifted_xCord << ", " <<
                    // shifted_yCord << "]" <<std::endl;
                    NumericType depth =
                        std::cos(radians) * zCord +
                        std::sin(radians) * (initialX - shifted_xCord);
                    NumericType lateralDisplacement =
                        std::sqrt(std::pow((std::cos(radians) *
                                                (initialX - shifted_xCord) -
                                            std::sin(radians) * zCord),
                                           2) +
                                  std::pow(initialY - shifted_yCord, 2));
                    std::array<NumericType, 3> coords{
                        shifted_xCord, shifted_yCord, shifted_zCord};
                    auto index = cellSet_->getIndex(coords);
                    if (index != -1) {
                      const auto contribution =
                          doseWeight * model_->getProfile(depth, lateralDisplacement);
                      (*concentration)[index] +=
                          contribution;
                      if (beamHits != nullptr && contribution > NumericType(0)) {
                        (*beamHits)[index] = NumericType(1);
                      }
                    }
                  }
                }
              }
              break;
            }
          }
        }
      }
    } else {
      // iterate over all the beams that hit the x-plane from the y direction:
      for (int i = 0; i < numberOfcellsXdirection; i++) {
        NumericType initialX = minX + i * gridDelta + gridDelta;
        NumericType initialY = maxY - gridDelta;
        std::array<NumericType, 3> initialCoords{initialX, initialY, 0};
        int initialIndex;
        numberOfcellsYdirection = yLength / gridDelta;
        while (true) {
          initialCoords[0] = initialX;
          initialCoords[1] = initialY;
          initialIndex = cellSet_->getIndex(initialCoords);
          if (initialIndex == -1) {
            initialY -= gridDelta * std::cos(radians);
            initialX -= gridDelta * std::sin(radians);
            numberOfcellsYdirection -= 1;
            if (initialY < minY - gridDelta || initialX < minX - xLength ||
                initialX > maxX + xLength) {
              break;
            }
            continue;
          } else {
            const auto materialId = static_cast<int>((*material)[initialIndex]);
            if (materialId == 0) {
              initialY -= gridDelta * std::cos(radians);
              initialX -= gridDelta * std::sin(radians);
              numberOfcellsYdirection -= 1;
              if (initialY < minY - gridDelta || initialX < minX - xLength ||
                  initialX > maxX + xLength) {
                break;
              }
              continue;
            }
            if (maskMaterialSet.count(materialId) > 0) {
              break;
            }

            // iterate over all the cells in y direction (depth), and then
            // iterate over all the cells in x direction to get the lateral
            // displacement
            for (int j = 0; j < relevant_Cells + 1; j++) {
              for (int k = 0; k < relevant_lateral_Cells + 1; k++) {
                NumericType yCord = j * gridDelta;
                NumericType xCord = initialX + k * gridDelta;
                NumericType shifted_xCord =
                    xCord - (relevant_lateral_Cells * gridDelta) / 2;
                NumericType shifted_yCord = initialY - yCord;
                // std::cout << "coord: [" << shifted_xCord << ", " <<
                // shifted_yCord << "]" <<std::endl;
                NumericType depth =
                    std::cos(radians) * yCord +
                    std::sin(radians) * (initialX - shifted_xCord);
                NumericType lateralDisplacement =
                    std::abs(std::cos(radians) * (initialX - shifted_xCord) -
                             std::sin(radians) * yCord);
                std::array<NumericType, 3> coords{shifted_xCord, shifted_yCord,
                                                  0};
                auto index = cellSet_->getIndex(coords);
                if (index != -1) {
                  const auto contribution =
                      doseWeight * model_->getProfile(depth, lateralDisplacement);
                  (*concentration)[index] +=
                      contribution;
                  if (beamHits != nullptr && contribution > NumericType(0)) {
                    (*beamHits)[index] = NumericType(1);
                  }
                }
              }
            }
            break;
          }
        }
      }
    }

    if (outputConcentrationInCm3_) {
      const auto scale = NumericType(1) /
                         std::pow(std::max(lengthUnitInCm_, NumericType(1e-30)),
                                  3);
      for (auto &value : *concentration)
        value *= scale;
    }
  }
};
} // namespace viennacs
