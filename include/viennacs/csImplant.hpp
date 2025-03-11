// apply the model to the geometry
#pragma once

#include "csImplantModel.hpp"
#include <csDenseCellSet.hpp>

#include <cmath>
#include <random>
#include <vcLogger.hpp>

namespace viennacs {

using namespace viennacore;

template <class NumericType, int D> class Implant {
  SmartPointer<DenseCellSet<NumericType, D>> cellSet_;
  SmartPointer<ImplantModel<NumericType, D>> model_;
  std::vector<int> maskMaterials;
  NumericType angle_;

public:
  Implant() = default;

  void setCellSet(SmartPointer<DenseCellSet<NumericType, D>> passedCellSet) {
    cellSet_ = passedCellSet;
  }

  void setImplantAngle(const NumericType angle) { angle_ = angle; }

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
    auto concentration = cellSet_->getScalarData("concentration");
    auto material = cellSet_->getScalarData("Material");

    NumericType xLength = std::abs(boundingBox[1][0] - boundingBox[0][0]);
    NumericType yLength = std::abs(boundingBox[1][1] - boundingBox[0][1]);
    NumericType zLength = std::abs(boundingBox[1][2] - boundingBox[0][2]);

    int numberOfcellsXdirection = xLength / gridDelta;
    int numberOfcellsYdirection = yLength / gridDelta;
    int numberOfcellsZdirection = zLength / gridDelta;

    double radians = angle_ * M_PI / 180;

    if constexpr (D == 3) {
      // iterate over all the beams that hit the xy-plane from the z direction:
      for (int i = 0; i < numberOfcellsXdirection; i++) {
        for (int h = 0; h < numberOfcellsYdirection; h++) {
          NumericType initialX = i * gridDelta - xLength / 2 + gridDelta;
          NumericType initialY = h * gridDelta - yLength / 2 + gridDelta;
          NumericType initialZ = zLength;
          std::array<NumericType, 3> initialCoords{initialX, initialY,
                                                   initialZ};
          int initialIndex;
          numberOfcellsZdirection = zLength / gridDelta;
          do {
            initialCoords[0] = initialX;
            initialCoords[1] = initialY;
            initialCoords[2] = initialZ;
            initialIndex = cellSet_->getIndex(initialCoords);
            if (initialIndex == -1) {
              initialZ -= gridDelta * std::cos(radians);
              initialX -= gridDelta * std::sin(radians);
              numberOfcellsZdirection -= 1;
              if (std::abs(initialCoords[2]) > zLength) {
                break;
              }
            } else {
              if ((*material)[initialIndex] == 0.) {
                break;
              }
              NumericType sigma = 1;         // params.get("sigma");
              NumericType lateral_sigma = 1; // params.get("lateral_sigma");
              int relevant_lateral_Cells = 3 * lateral_sigma / gridDelta;
              int relevant_Cells = 3 * sigma / gridDelta;
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
                      (*concentration)[index] +=
                          model_->getDepthProfile(depth) *
                          model_->getLateralProfile(lateralDisplacement, depth);
                    }
                  }
                }
              }
            }
          } while (initialIndex == -1);
        }
      }
    } else {
      // iterate over all the beams that hit the x-plane from the y direction:
      for (int i = 0; i < numberOfcellsXdirection; i++) {
        NumericType initialX = i * gridDelta - xLength / 2 + gridDelta;
        NumericType initialY = yLength - gridDelta;
        std::array<NumericType, 3> initialCoords{initialX, initialY, 0};
        int initialIndex;
        numberOfcellsYdirection = yLength / gridDelta;
        do {
          initialCoords[0] = initialX;
          initialCoords[1] = initialY;
          initialIndex = cellSet_->getIndex(initialCoords);
          if (initialIndex == -1) {
            initialY -= gridDelta * std::cos(radians);
            initialX -= gridDelta * std::sin(radians);
            numberOfcellsYdirection -= 1;
            if (std::abs(initialCoords[1]) > yLength) {
              break;
            }
          } else {
            if ((*material)[initialIndex] == 0.) {
              break;
            } else {
              // iterate over all the cells in y direction (depth), and then
              // iterate over all the cells in x direction to get the lateral
              // displacement
              for (int j = 0; j < numberOfcellsYdirection + 2; j++) {
                for (int k = 0; k < numberOfcellsXdirection + 1; k++) {
                  NumericType yCord = j * gridDelta;
                  NumericType xCord = k * gridDelta;
                  NumericType shifted_xCord = xCord - xLength / 2;
                  NumericType shifted_yCord = initialY - yCord;
                  // std::cout << "coord: [" << shifted_xCord << ", " <<
                  // shifted_yCord << "]" <<std::endl;
                  NumericType depth =
                      std::cos(radians) * yCord +
                      std::sin(radians) * (initialX - shifted_xCord);
                  NumericType lateralDisplacement =
                      std::abs(std::cos(radians) * (initialX - shifted_xCord) -
                               std::sin(radians) * yCord);
                  std::array<NumericType, 3> coords{shifted_xCord,
                                                    shifted_yCord, 0};
                  auto index = cellSet_->getIndex(coords);
                  if (index != -1) {
                    (*concentration)[index] +=
                        model_->getDepthProfile(depth) *
                        model_->getLateralProfile(lateralDisplacement, depth);
                  }
                }
              }
            }
          }
        } while (initialIndex == -1);
      }
    }
  }
};
} // namespace viennacs