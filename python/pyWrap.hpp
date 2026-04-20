/*
  This file is used to generate the python module of ViennaCS.
  It uses pybind11 to create the modules.
*/

#define PYBIND11_DETAILED_ERROR_MESSAGES
#define VIENNACS_PYTHON_BUILD

// correct module name macro
#define TOKENPASTE_INTERNAL(x, y, z) x##y##z
#define TOKENPASTE(x, y, z) TOKENPASTE_INTERNAL(x, y, z)
#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)
#define VIENNACS_MODULE_VERSION STRINGIZE(VIENNACS_VERSION)

#include <optional>
#include <vector>

#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// all header files which define API functions
#include <csAtomicLayerProcess.hpp>
#include <csConstants.hpp>
#include <csDenseCellSet.hpp>
#include <csImplant.hpp>
#include <csImplantModel.hpp>
#include <csMeanFreePath.hpp>
#include <csSegmentCells.hpp>
#include <csVersion.hpp>
#include <models/csImplantGaussian.hpp>
#include <models/csImplantPearson.hpp>

using namespace viennacs;

// always use double for python export
typedef double T;
namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(Types, SmartPointer<Types>)

template <int D> void bindAPI(py::module &module) {
  // DenseCellSet
  py::class_<DenseCellSet<T, D>, SmartPointer<DenseCellSet<T, D>>>(
      module, "DenseCellSet")
      .def(py::init())
      .def("fromLevelSets", &DenseCellSet<T, D>::fromLevelSets,
           py::arg("levelSets"), py::arg("materialMap") = nullptr,
           py::arg("depth") = 0.)
      .def("getBoundingBox", &DenseCellSet<T, D>::getBoundingBox)
      .def(
          "addScalarData",
          [](DenseCellSet<T, D> &cellSet, std::string name, T initValue) {
            cellSet.addScalarData(name, initValue);
            // discard return value
          },
          "Add a scalar value to be stored and modified in each cell.")
      .def(
          "addVectorData",
          [](DenseCellSet<T, D> &cellSet, std::string name,
             std::array<T, 3> initValue) {
            cellSet.addVectorData(name, initValue);
            // discard return value
          },
          py::arg("name"), py::arg("initValue") = std::array<T, 3>{0., 0., 0.},
          "Add a vector value (3 components) to be stored and modified in each "
          "cell.")
      .def("getDepth", &DenseCellSet<T, D>::getDepth,
           "Get the depth of the cell set.")
      .def("getGridDelta", &DenseCellSet<T, D>::getGridDelta,
           "Get the cell size.")
      .def("getNodes", &DenseCellSet<T, D>::getNodes,
           "Get the nodes of the cell set which correspond to the corner "
           "points of the cells.")
      .def("getNode", &DenseCellSet<T, D>::getNode,
           "Get the node at the given index.")
      .def("getElements", &DenseCellSet<T, D>::getElements,
           "Get elements (cells). The indicies in the elements correspond to "
           "the corner nodes.")
      .def("getElement", &DenseCellSet<T, D>::getElement,
           "Get the element at the given index.")
      .def("getSurface", &DenseCellSet<T, D>::getSurface,
           "Get the surface level-set.")
      .def("getCellGrid", &DenseCellSet<T, D>::getCellGrid,
           "Get the underlying mesh of the cell set.")
      .def("getLevelSets", &DenseCellSet<T, D>::getLevelSets,
           "Get the level sets used to construct the cell set.")
      .def("getNumberOfCells", &DenseCellSet<T, D>::getNumberOfCells,
           "Get the number of cells.")
      .def("getFillingFraction", &DenseCellSet<T, D>::getFillingFraction,
           "Get the filling fraction of the cell containing the point.")
      .def("getFillingFractions", &DenseCellSet<T, D>::getFillingFractions,
           "Get the filling fractions of all cells.")
      .def("getAverageFillingFraction",
           &DenseCellSet<T, D>::getAverageFillingFraction,
           "Get the average filling at a point in some radius.")
      .def("getCellCenter", &DenseCellSet<T, D>::getCellCenter,
           "Get the center of a cell with given index")
      .def("getScalarData", &DenseCellSet<T, D>::getScalarData,
           "Get the data stored at each cell. WARNING: This function only "
           "returns a copy of the data")
      .def("getVectorData", &DenseCellSet<T, D>::getVectorData,
           "Get the vector data stored at each cell. WARNING: This function "
           "only "
           "returns a copy of the data")
      .def("setScalarData", &DenseCellSet<T, D>::setScalarData, py::arg("name"),
           py::arg("newData"),
           "Overwrite the scalar data associated with 'name' with a new array.")
      .def("setVectorData", &DenseCellSet<T, D>::setVectorData, py::arg("name"),
           py::arg("newData"),
           "Overwrite the vector data associated with 'name' with a new array.")
      .def("getScalarDataLabels", &DenseCellSet<T, D>::getScalarDataLabels,
           "Get the labels of the scalar data stored in the cell set.")
      .def("getIndex", &DenseCellSet<T, D>::getIndex,
           "Get the index of the cell containing the given point.")
      .def("setCellSetPosition", &DenseCellSet<T, D>::setCellSetPosition,
           "Set whether the cell set should be created below (false) or above "
           "(true) the surface.")
      .def("getCellSetPosition", &DenseCellSet<T, D>::getCellSetPosition,
           "Get whether the cell set is created below or above the surface.")
      .def(
          "setCoverMaterial", &DenseCellSet<T, D>::setCoverMaterial,
          "Set the material of the cells which are above or below the surface.")
      .def("setPeriodicBoundary", &DenseCellSet<T, D>::setPeriodicBoundary,
           "Enable periodic boundary conditions in specified dimensions.")
      .def("setFillingFraction",
           py::overload_cast<const int, const T>(
               &DenseCellSet<T, D>::setFillingFraction),
           "Sets the filling fraction at given cell index.")
      .def("setFillingFraction",
           py::overload_cast<const std::array<T, 3> &, const T>(
               &DenseCellSet<T, D>::setFillingFraction),
           "Sets the filling fraction for cell which contains given point.")
      .def("addFillingFraction",
           py::overload_cast<const int, const T>(
               &DenseCellSet<T, D>::addFillingFraction),
           "Add to the filling fraction at given cell index.")
      .def("addFillingFraction",
           py::overload_cast<const std::array<T, 3> &, const T>(
               &DenseCellSet<T, D>::addFillingFraction),
           "Add to the filling fraction for cell which contains given point.")
      .def("addFillingFractionInMaterial",
           &DenseCellSet<T, D>::addFillingFractionInMaterial,
           "Add to the filling fraction for cell which contains given point "
           "only if the cell has the specified material ID.")
      .def("writeVTU", &DenseCellSet<T, D>::writeVTU,
           "Write the cell set as .vtu file")
      .def("writeCellSetData", &DenseCellSet<T, D>::writeCellSetData,
           "Save cell set data in simple text format.")
      .def("readCellSetData", &DenseCellSet<T, D>::readCellSetData,
           "Read cell set data from text.")
      .def("clear", &DenseCellSet<T, D>::clear, "Clear the filling fractions.")
      .def("updateMaterials", &DenseCellSet<T, D>::updateMaterials,
           "Update the material IDs of the cell set. This function should be "
           "called if the level sets, the cell set is made out of, have "
           "changed. This does not work if the surface of the volume has "
           "changed. In this case, call the function 'updateSurface' first.")
      .def("updateSurface", &DenseCellSet<T, D>::updateSurface,
           "Updates the surface of the cell set. The new surface should be "
           "below the old surface as this function can only remove cells from "
           "the cell set.")
      .def("buildNeighborhood", &DenseCellSet<T, D>::buildNeighborhood,
           "Generate fast neighbor access for each cell.",
           py::arg("forceRebuild") = false)
      .def("getNeighbors", &DenseCellSet<T, D>::getNeighbors,
           "Get the neighbor indices for a cell.");

  // SegmentCells
  py::class_<SegmentCells<T, D>, SmartPointer<SegmentCells<T, D>>>(
      module, "SegmentCells")
      .def(py::init<SmartPointer<DenseCellSet<T, D>>>())
      .def(py::init<SmartPointer<DenseCellSet<T, D>>, std::string, int>(),
           py::arg("cellSet"), py::arg("cellTypeString") = "CellType",
           py::arg("bulkMaterial") = 1)
      .def("setCellSet", &SegmentCells<T, D>::setCellSet,
           "Set the cell set in the segmenter.")
      .def("setCellTypeString", &SegmentCells<T, D>::setCellTypeString,
           "Set the cell type string in the segmenter.")
      .def("setBulkMaterial", &SegmentCells<T, D>::setBulkMaterial,
           "Set the bulk material in the segmenter.")
      .def("apply", &SegmentCells<T, D>::apply,
           "Segment the cells into surface, material, and gas cells.");

  // Mean Free Path
  py::class_<MeanFreePath<T, D>, SmartPointer<MeanFreePath<T, D>>>(
      module, "MeanFreePath")
      .def(py::init<SmartPointer<DenseCellSet<T, D>>>())
      .def("setNumRaysPerCell", &MeanFreePath<T, D>::setNumRaysPerCell)
      .def("setReflectionLimit", &MeanFreePath<T, D>::setReflectionLimit)
      .def("setRngSeed", &MeanFreePath<T, D>::setRngSeed)
      .def("setMaterial", &MeanFreePath<T, D>::setMaterial)
      .def("setBulkLambda", &MeanFreePath<T, D>::setBulkLambda)
      .def("enableSmoothing", &MeanFreePath<T, D>::enableSmoothing)
      .def("disableSmoothing", &MeanFreePath<T, D>::disableSmoothing)
      .def("apply", &MeanFreePath<T, D>::apply);

  py::class_<typename AtomicLayerProcess<T, D>::Precursor>(module, "Precursor")
      .def(py::init<>())
      .def_readwrite("name", &AtomicLayerProcess<T, D>::Precursor::name)
      .def_readwrite("meanThermalVelocity",
                     &AtomicLayerProcess<T, D>::Precursor::meanThermalVelocity)
      .def_readwrite("adsorptionRate",
                     &AtomicLayerProcess<T, D>::Precursor::adsorptionRate)
      .def_readwrite("desorptionRate",
                     &AtomicLayerProcess<T, D>::Precursor::desorptionRate)
      .def_readwrite("duration", &AtomicLayerProcess<T, D>::Precursor::duration)
      .def_readwrite("inFlux", &AtomicLayerProcess<T, D>::Precursor::inFlux);

  // Atomic Layer Process
  py::class_<AtomicLayerProcess<T, D>, SmartPointer<AtomicLayerProcess<T, D>>>(
      module, "AtomicLayerProcess")
      .def(py::init<SmartPointer<DenseCellSet<T, D>>, const bool>(),
           py::arg("cellSet"), py::arg("etch") = false)
      .def("setFirstPrecursor",
           py::overload_cast<std::string, T, T, T, T, T>(
               &AtomicLayerProcess<T, D>::setFirstPrecursor))
      .def("setFirstPrecursor",
           py::overload_cast<
               const typename AtomicLayerProcess<T, D>::Precursor &>(
               &AtomicLayerProcess<T, D>::setFirstPrecursor))
      .def("setSecondPrecursor",
           py::overload_cast<std::string, T, T, T, T, T>(
               &AtomicLayerProcess<T, D>::setSecondPrecursor))
      .def("setSecondPrecursor",
           py::overload_cast<
               const typename AtomicLayerProcess<T, D>::Precursor &>(
               &AtomicLayerProcess<T, D>::setSecondPrecursor))
      .def("setPurgeParameters", &AtomicLayerProcess<T, D>::setPurgeParameters)
      .def("setReactionOrder", &AtomicLayerProcess<T, D>::setReactionOrder)
      .def("setMaxLambda", &AtomicLayerProcess<T, D>::setMaxLambda)
      .def("setStabilityFactor", &AtomicLayerProcess<T, D>::setStabilityFactor)
      .def("setMaxTimeStep", &AtomicLayerProcess<T, D>::setMaxTimeStep)
      .def("setPrintInterval", &AtomicLayerProcess<T, D>::setPrintInterval)
      .def("apply", &AtomicLayerProcess<T, D>::apply);

  py::class_<ImplantModel<T, D>, SmartPointer<ImplantModel<T, D>>>(
      module, "ImplantModel")
      .def("getDepthProfile", &ImplantModel<T, D>::getDepthProfile)
      .def("getLateralProfile", &ImplantModel<T, D>::getLateralProfile)
      .def("getMaxDepth", &ImplantModel<T, D>::getMaxDepth)
      .def("getMaxLateralRange", &ImplantModel<T, D>::getMaxLateralRange);

  py::class_<ImplantGaussian<T, D>, ImplantModel<T, D>,
             SmartPointer<ImplantGaussian<T, D>>>(module, "ImplantGaussian")
      .def(py::init<T, T, T, T>(), py::arg("mu"), py::arg("sigma"),
           py::arg("lateralSigma"), py::arg("lateralMu"));

  py::class_<ImplantPearsonIV<T, D>, ImplantModel<T, D>,
             SmartPointer<ImplantPearsonIV<T, D>>>(module,
                                                   "ImplantPearsonIV")
      .def(py::init<const constants::PearsonIVParameters<T> &, T, T>(),
           py::arg("params"), py::arg("lateralMu"),
           py::arg("lateralSigma"));

  py::class_<ImplantPearsonIVChanneling<T, D>, ImplantModel<T, D>,
             SmartPointer<ImplantPearsonIVChanneling<T, D>>>(
      module, "ImplantPearsonIVChanneling")
      .def(py::init<const constants::PearsonIVParameters<T> &, T, T, T, T, T,
                    T>(),
           py::arg("params"), py::arg("lateralMu"),
           py::arg("lateralSigma"), py::arg("tailFraction"),
           py::arg("tailStartDepth"), py::arg("tailDecayLength"),
           py::arg("tailBlendWidth") = T(0));

  py::class_<ImplantDualPearsonIV<T, D>, ImplantModel<T, D>,
             SmartPointer<ImplantDualPearsonIV<T, D>>>(
      module, "ImplantDualPearsonIV")
      .def(py::init<const constants::PearsonIVParameters<T> &,
                    const constants::PearsonIVParameters<T> &, T, T, T>(),
           py::arg("headParams"), py::arg("tailParams"),
           py::arg("headFraction"), py::arg("lateralMu"),
           py::arg("lateralSigma"));

  py::class_<Implant<T, D>>(module, "Implant")
      .def(py::init<>())
      .def("setCellSet", &Implant<T, D>::setCellSet)
      .def("setImplantAngle", &Implant<T, D>::setImplantAngle)
      .def("setImplantModel", &Implant<T, D>::setImplantModel)
      .def(
          "setMaskMaterials",
          [](Implant<T, D> &implant, const std::vector<int> &materials) {
            switch (materials.size()) {
            case 0:
              implant.setMaskMaterials();
              break;
            case 1:
              implant.setMaskMaterials(materials[0]);
              break;
            case 2:
              implant.setMaskMaterials(materials[0], materials[1]);
              break;
            case 3:
              implant.setMaskMaterials(materials[0], materials[1],
                                       materials[2]);
              break;
            case 4:
              implant.setMaskMaterials(materials[0], materials[1],
                                       materials[2], materials[3]);
              break;
            case 5:
              implant.setMaskMaterials(materials[0], materials[1],
                                       materials[2], materials[3],
                                       materials[4]);
              break;
            case 6:
              implant.setMaskMaterials(materials[0], materials[1],
                                       materials[2], materials[3],
                                       materials[4], materials[5]);
              break;
            case 7:
              implant.setMaskMaterials(materials[0], materials[1],
                                       materials[2], materials[3],
                                       materials[4], materials[5],
                                       materials[6]);
              break;
            case 8:
              implant.setMaskMaterials(materials[0], materials[1],
                                       materials[2], materials[3],
                                       materials[4], materials[5],
                                       materials[6], materials[7]);
              break;
            default:
              throw std::runtime_error(
                  "setMaskMaterials currently supports up to 8 materials.");
            }
          },
          py::arg("materials"),
          "Set the material IDs to be treated as mask materials.")
      .def("apply", &Implant<T, D>::apply);
}
