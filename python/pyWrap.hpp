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
#include <csAnneal.hpp>
#include <csAtomicLayerProcess.hpp>
#include <csDenseCellSet.hpp>
#include <csImplant.hpp>
#include <csImplantModel.hpp>
#include <csMeanFreePath.hpp>
#include <csNetDoping.hpp>
#include <csSegmentCells.hpp>
#include <csSheetResistance.hpp>
#include <csVersion.hpp>

using namespace viennacs;

// always use double for python export
typedef double T;
namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(Types, SmartPointer<Types>)

template <int D> void bindAPI(py::module &module) {

  // EmbeddedBoundaryPoint struct (read-only; produced by DenseCellSet)
  using EBP = typename DenseCellSet<T, D>::EmbeddedBoundaryPoint;
  py::class_<EBP>(module, "EmbeddedBoundaryPoint")
      .def(py::init<>())
      .def_readwrite("coordinate", &EBP::coordinate)
      .def_readwrite("normal", &EBP::normal)
      .def_readwrite("signedDistance", &EBP::signedDistance)
      .def_readwrite("edgeFraction", &EBP::edgeFraction)
      .def_readwrite("levelSetIndex", &EBP::levelSetIndex)
      .def_readwrite("negativeMaterial", &EBP::negativeMaterial)
      .def_readwrite("positiveMaterial", &EBP::positiveMaterial)
      .def_readwrite("adjacentCell", &EBP::adjacentCell)
      .def_readwrite("axis", &EBP::axis);

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
           "Get the neighbor indices for a cell.")
      .def("enableEmbeddedBoundaries",
           &DenseCellSet<T, D>::enableEmbeddedBoundaries,
           py::arg("enable") = true,
           "Enable sub-grid embedded boundary point generation. Must be called "
           "before fromLevelSets().")
      .def("embeddedBoundariesEnabled",
           &DenseCellSet<T, D>::embeddedBoundariesEnabled,
           "Return true if embedded boundary generation was enabled.")
      .def("hasEmbeddedBoundaries", &DenseCellSet<T, D>::hasEmbeddedBoundaries,
           "Return true if the cell set has any embedded boundary points.")
      .def("numEmbeddedBoundaryPoints",
           &DenseCellSet<T, D>::numEmbeddedBoundaryPoints,
           "Return the total number of embedded boundary points.")
      .def("getEmbeddedBoundaryPoints",
           &DenseCellSet<T, D>::getEmbeddedBoundaryPoints,
           "Return the list of all EmbeddedBoundaryPoint objects.")
      .def("getEmbeddedBoundaryPointIds",
           &DenseCellSet<T, D>::getEmbeddedBoundaryPointIds, py::arg("cellIdx"),
           "Return the indices into getEmbeddedBoundaryPoints() for the given "
           "cell.")
      .def("getFaceBoundaryPointId",
           &DenseCellSet<T, D>::getFaceBoundaryPointId, py::arg("cellIdx"),
           py::arg("faceIdx"),
           "Index into getEmbeddedBoundaryPoints() for the boundary point on "
           "face faceIdx (axis*2 + (offset>0?1:0)) of cellIdx, or -1 if none.")
      .def("getFaceBoundaryDistance",
           &DenseCellSet<T, D>::getFaceBoundaryDistance, py::arg("cellIdx"),
           py::arg("faceIdx"),
           "Distance from the cell center to the boundary point on face "
           "faceIdx. Returns gridDelta/2 when no point exists.")
      .def("getMinFaceBoundaryDistance",
           &DenseCellSet<T, D>::getMinFaceBoundaryDistance,
           "Minimum face-boundary distance across all cells. Returns "
           "gridDelta/2 when no embedded boundaries are present.");

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
      .def("getProfile", &ImplantModel<T, D>::getProfile)
      .def("getMaxDepth", &ImplantModel<T, D>::getMaxDepth)
      .def("getMaxLateralRange", &ImplantModel<T, D>::getMaxLateralRange);

  py::class_<Implant<T, D>>(module, "Implant")
      .def(py::init<>())
      .def("setCellSet", &Implant<T, D>::setCellSet)
      .def("setImplantAngle", &Implant<T, D>::setImplantAngle)
      .def("setDose", &Implant<T, D>::setDose)
      .def("setLengthUnitInCm", &Implant<T, D>::setLengthUnitInCm)
      .def("setDoseControl", &Implant<T, D>::setDoseControl)
      .def("enableBeamHits", &Implant<T, D>::enableBeamHits,
           py::arg("enable") = true)
      .def("setConcentrationLabel", &Implant<T, D>::setConcentrationLabel)
      .def("setBeamHitsLabel", &Implant<T, D>::setBeamHitsLabel)
      .def("setDamageLabel", &Implant<T, D>::setDamageLabel)
      .def("setLastDamageLabel", &Implant<T, D>::setLastDamageLabel)
      .def("setDamageFactor", &Implant<T, D>::setDamageFactor)
      .def("setOutputConcentrationInCm3",
           &Implant<T, D>::setOutputConcentrationInCm3,
           py::arg("enable") = true)
      .def("setImplantModel", &Implant<T, D>::setImplantModel)
      .def("setDamageModel", &Implant<T, D>::setDamageModel)
      .def(
          "setMaskMaterials",
          [](Implant<T, D> &implant, const std::vector<int> &materials) {
            implant.setMaskMaterials(materials);
          },
          py::arg("materials"),
          "Set the material IDs to be treated as mask materials.")
      .def(
          "setScreenMaterials",
          [](Implant<T, D> &implant, const std::vector<int> &materials) {
            implant.setScreenMaterials(materials);
          },
          py::arg("materials"),
          "Set the material IDs to be treated as screen/cap materials.")
      .def(
          "setVoidMaterials",
          [](Implant<T, D> &implant, const std::vector<int> &materials) {
            implant.setVoidMaterials(materials);
          },
          py::arg("materials"),
          "Set material IDs treated as vacuum/void (beam passes through). "
          "Default: {0}.")
      .def("setVoidMaterial", &Implant<T, D>::setVoidMaterial,
           py::arg("material"),
           "Set a single material ID treated as vacuum/void (beam passes "
           "through). Default: 0.")
      .def("apply", &Implant<T, D>::apply);

  py::class_<typename Anneal<T, D>::DefectDiagnosticsRow>(
      module, "DefectDiagnosticsRow")
      .def(py::init<>())
      .def_readwrite("step", &Anneal<T, D>::DefectDiagnosticsRow::step)
      .def_readwrite("time_s", &Anneal<T, D>::DefectDiagnosticsRow::time_s)
      .def_readwrite("temperature_K",
                     &Anneal<T, D>::DefectDiagnosticsRow::temperature_K)
      .def_readwrite("I_mean", &Anneal<T, D>::DefectDiagnosticsRow::I_mean)
      .def_readwrite("V_mean", &Anneal<T, D>::DefectDiagnosticsRow::V_mean)
      .def_readwrite("I_min", &Anneal<T, D>::DefectDiagnosticsRow::I_min)
      .def_readwrite("I_max", &Anneal<T, D>::DefectDiagnosticsRow::I_max)
      .def_readwrite("V_min", &Anneal<T, D>::DefectDiagnosticsRow::V_min)
      .def_readwrite("V_max", &Anneal<T, D>::DefectDiagnosticsRow::V_max)
      .def_readwrite("I_over_Ieq_mean",
                     &Anneal<T, D>::DefectDiagnosticsRow::I_over_Ieq_mean)
      .def_readwrite("V_over_Veq_mean",
                     &Anneal<T, D>::DefectDiagnosticsRow::V_over_Veq_mean)
      .def_readwrite("IV_over_IeqVeq_mean",
                     &Anneal<T, D>::DefectDiagnosticsRow::IV_over_IeqVeq_mean)
      .def_readwrite("Ieq", &Anneal<T, D>::DefectDiagnosticsRow::Ieq)
      .def_readwrite("Veq", &Anneal<T, D>::DefectDiagnosticsRow::Veq);

  py::class_<Anneal<T, D>>(module, "Anneal")
      .def(py::init<>())
      .def("setCellSet", &Anneal<T, D>::setCellSet)
      .def("setSpeciesLabel", &Anneal<T, D>::setSpeciesLabel)
      .def("setMaterialLabel", &Anneal<T, D>::setMaterialLabel)
      .def("setDuration", &Anneal<T, D>::setDuration)
      .def("setTimeStep", &Anneal<T, D>::setTimeStep)
      .def("setStabilityFactor", &Anneal<T, D>::setStabilityFactor)
      .def("setClampNonNegative", &Anneal<T, D>::setClampNonNegative,
           py::arg("enable") = true)
      .def("setMode", &Anneal<T, D>::setMode)
      .def("setImplicitSolverOptions", &Anneal<T, D>::setImplicitSolverOptions)
      .def("setDiffusionCoefficient", &Anneal<T, D>::setDiffusionCoefficient)
      .def("setArrheniusParameters", &Anneal<T, D>::setArrheniusParameters)
      .def("setTemperature", &Anneal<T, D>::setTemperature)
      .def("clearTemperatureSchedule", &Anneal<T, D>::clearTemperatureSchedule)
      .def("addIsothermalStep", &Anneal<T, D>::addIsothermalStep)
      .def("addRampStep", &Anneal<T, D>::addRampStep)
      .def("setTemperatureSchedule", &Anneal<T, D>::setTemperatureSchedule,
           py::arg("durations"), py::arg("temperatures"),
           "Set a temperature schedule from lists of durations (s) and "
           "temperatures (K). "
           "N temperatures → N isothermal steps; N+1 temperatures → N ramp "
           "steps.")
      .def("setDiffusionMaterials", &Anneal<T, D>::setDiffusionMaterials)
      .def("setBlockingMaterials", &Anneal<T, D>::setBlockingMaterials)
      .def("enableDefectCoupling", &Anneal<T, D>::enableDefectCoupling,
           py::arg("enable") = true)
      .def("setDamageLabels", &Anneal<T, D>::setDamageLabels)
      .def("setDefectLabels", &Anneal<T, D>::setDefectLabels)
      .def("resetDefectInitialization",
           &Anneal<T, D>::resetDefectInitialization)
      .def("setDefectSourceWeights", &Anneal<T, D>::setDefectSourceWeights)
      .def("setDefectPartition", &Anneal<T, D>::setDefectPartition)
      .def("setDefectPartitionFactors",
           &Anneal<T, D>::setDefectPartitionFactors)
      .def("setDefectDiffusivities", &Anneal<T, D>::setDefectDiffusivities)
      .def("setDefectReactionRates", &Anneal<T, D>::setDefectReactionRates)
      .def("enableDefectEquilibrium", &Anneal<T, D>::enableDefectEquilibrium,
           py::arg("enable") = true)
      .def("setDefectEquilibrium", &Anneal<T, D>::setDefectEquilibrium)
      .def("setDefectEquilibriumArrhenius",
           &Anneal<T, D>::setDefectEquilibriumArrhenius)
      .def("clearEquilibriumArrhenius",
           &Anneal<T, D>::clearEquilibriumArrhenius)
      .def("setDefectEnhancedDiffusion",
           &Anneal<T, D>::setDefectEnhancedDiffusion)
      .def("setTEDFromDamageFactor", &Anneal<T, D>::setTEDFromDamageFactor,
           py::arg("damageFactor"), py::arg("coefficientScale") = T(0.5),
           py::arg("normalization") = T(1e20))
      .def("enableDefectClustering", &Anneal<T, D>::enableDefectClustering,
           py::arg("enable") = true)
      .def("setDefectClusterLabel", &Anneal<T, D>::setDefectClusterLabel)
      .def("setDefectClusterKinetics", &Anneal<T, D>::setDefectClusterKinetics)
      .def("setDefectClusterInitFraction",
           &Anneal<T, D>::setDefectClusterInitFraction)
      .def("enableDiagnostics", &Anneal<T, D>::enableDiagnostics,
           py::arg("enable") = true)
      .def("setDiagnosticsMaterialFilter",
           &Anneal<T, D>::setDiagnosticsMaterialFilter,
           py::arg("materialId") = -1)
      .def("clearDefectDiagnostics", &Anneal<T, D>::clearDefectDiagnostics)
      .def("getDefectDiagnostics", &Anneal<T, D>::getDefectDiagnostics)
      .def("diffusivity", &Anneal<T, D>::diffusivity)
      .def("enableSolidActivation", &Anneal<T, D>::enableSolidActivation,
           py::arg("enable") = true,
           "Enable the solid solubility activation model. When active, writes "
           "'active_concentration' field: C_A+ = C_SS*C/(C_SS+C).")
      .def("setSolidSolubilityArrhenius",
           &Anneal<T, D>::setSolidSolubilityArrhenius, py::arg("C0"),
           py::arg("Ea_eV"),
           "Set the solid solubility Arrhenius parameters manually. "
           "C0 must be in nm⁻³ (same units as the concentration field).")
      .def("setActiveLabel", &Anneal<T, D>::setActiveLabel, py::arg("label"),
           "Set the cell-set field name for the active concentration output "
           "(default: 'active_concentration').")
      .def("setSurfaceBoundaryCondition",
           &Anneal<T, D>::setSurfaceBoundaryCondition, py::arg("condition"),
           "Set the embedded boundary condition applied at level-set surfaces "
           "(default: zero-flux Neumann). Only used when the cell set has "
           "embedded boundaries.")
      .def("setSourceField",
           py::overload_cast<std::vector<T>>(&Anneal<T, D>::setSourceField),
           py::arg("source"),
           "Per-cell volumetric source term S added to dc/dt = D∇²c + S. "
           "Pass an empty list to clear.")
      .def("clearSourceField", &Anneal<T, D>::clearSourceField,
           "Remove any previously set concentration source field.")
      .def("apply", &Anneal<T, D>::apply)
      .def("applyActivation", &Anneal<T, D>::applyActivation,
           "Apply only the solid-activation model (C_active = C_SS·C/(C_SS+C)) "
           "without running the diffusion solver.  Equivalent to Sentaurus "
           "'diffuse time=0'.  Requires setCellSet(), enableSolidActivation(), "
           "and setSolidSolubilityArrhenius() to be configured first.");

  // ── SheetResistance
  // ─────────────────────────────────────────────────────────
  py::class_<SheetResistance<T, D>>(
      module, "SheetResistance",
      "Compute sheet resistance (Rsh, Ω/□) from a concentration field stored "
      "in\n"
      "a DenseCellSet.\n\n"
      "Default configuration targets ViennaPS nm-unit domains:\n"
      "  length unit = 1e-7 (nm → cm),  conc unit = 1e21 (nm⁻³ → cm⁻³),\n"
      "  depth axis  = D−1  (y for 2-D, z for 3-D),\n"
      "  surface position = 0  (depth = surface − coordinate).\n\n"
      "Example::\n\n"
      "  sr = SheetResistance()\n"
      "  sr.setCellSet(domain.getCellSet())\n"
      "  sr.setConcentrationLabel(\"P_active\")\n"
      "  rsh = sr.computeElectron()   # Masetti n-type (P in Si)")
      .def(py::init<>())
      .def("setCellSet", &SheetResistance<T, D>::setCellSet, py::arg("cellSet"),
           "Attach the DenseCellSet to analyse.")
      .def("setConcentrationLabel",
           &SheetResistance<T, D>::setConcentrationLabel, py::arg("label"),
           "Name of the scalar field holding the active concentration "
           "(default: 'active_concentration').")
      .def("setDepthAxis", &SheetResistance<T, D>::setDepthAxis,
           py::arg("axis"), "Cell-centre axis index for depth  (default: D−1).")
      .def("setSurfacePosition", &SheetResistance<T, D>::setSurfacePosition,
           py::arg("surfacePosition"),
           "Wafer-surface coordinate along the depth axis. Depth is computed "
           "as surfacePosition minus the cell-centre coordinate.")
      .def("setLengthUnit", &SheetResistance<T, D>::setLengthUnit,
           py::arg("lu_cm"),
           "Conversion factor from domain length unit to cm "
           "(default: 1e-7 for nm domains). "
           "Also updates the concentration unit to stay consistent.")
      .def("setConcentrationUnit", &SheetResistance<T, D>::setConcentrationUnit,
           py::arg("unit"),
           "Multiplicative factor to convert the cell-set concentration to "
           "cm⁻³ (default: 1e21 for nm⁻³ fields).")
      .def("computeElectron", &SheetResistance<T, D>::computeElectron,
           "Rsh [Ω/□] using the Masetti-Severi electron mobility model "
           "(n-type, e.g. P-doped Si).")
      .def("computeHole", &SheetResistance<T, D>::computeHole,
           "Rsh [Ω/□] using the Masetti-Severi hole mobility model "
           "(p-type, e.g. B-doped Si).");

  // ── NetDoping
  // ────────────────────────────────────────────────────────────────
  py::class_<NetDoping<T, D>>(
      module, "NetDoping",
      "Compute net doping (Σ donors − Σ acceptors) and extract the\n"
      "metallurgical junction depth from a DenseCellSet.\n\n"
      "Donor and acceptor labels typically refer to active-concentration\n"
      "fields written by the Anneal solid-activation model (e.g. 'P_active',\n"
      "'B_active'), but total-concentration fields work too.\n\n"
      "Example::\n\n"
      "  nd = NetDoping()\n"
      "  nd.setCellSet(domain.getCellSet())\n"
      "  nd.addDonorLabel('P_active')\n"
      "  nd.addAcceptorLabel('B_active')\n"
      "  nd.apply()                      # writes 'net_doping' field\n"
      "  xj = nd.junctionDepth()         # nm from surface")
      .def(py::init<>())
      .def("setCellSet", &NetDoping<T, D>::setCellSet, py::arg("cellSet"),
           "Attach the DenseCellSet to analyse.")
      .def("addDonorLabel", &NetDoping<T, D>::addDonorLabel, py::arg("label"),
           "Append one donor (n-type) concentration field name.")
      .def("addAcceptorLabel", &NetDoping<T, D>::addAcceptorLabel,
           py::arg("label"),
           "Append one acceptor (p-type) concentration field name.")
      .def("setDonorLabels", &NetDoping<T, D>::setDonorLabels,
           py::arg("labels"), "Replace the donor label list.")
      .def("setAcceptorLabels", &NetDoping<T, D>::setAcceptorLabels,
           py::arg("labels"), "Replace the acceptor label list.")
      .def("setOutputLabel", &NetDoping<T, D>::setOutputLabel, py::arg("label"),
           "Name of the output scalar field (default: 'net_doping').")
      .def("setDepthAxis", &NetDoping<T, D>::setDepthAxis, py::arg("axis"),
           "Cell-centre axis index for depth (default: D−1).")
      .def("setSurfacePosition", &NetDoping<T, D>::setSurfacePosition,
           py::arg("surfacePosition"),
           "Wafer-surface coordinate along the depth axis. Depth is computed "
           "as surfacePosition minus the cell-centre coordinate.")
      .def("apply", &NetDoping<T, D>::apply,
           "Compute net_doping = Σ donors − Σ acceptors and write to the "
           "output field.")
      .def("junctionDepth", &NetDoping<T, D>::junctionDepth,
           "Shallowest depth [domain length units] where net_doping changes "
           "sign.  Returns inf if no junction exists.")
      .def("junctionDepths", &NetDoping<T, D>::junctionDepths,
           "All junction depths [domain length units], sorted ascending.  "
           "Useful for retrograde profiles with multiple sign changes.")
      .def("junctionCount", &NetDoping<T, D>::junctionCount,
           "Number of metallurgical junctions (net_doping sign changes).")
      .def("lateralJunctionPosition", &NetDoping<T, D>::lateralJunctionPosition,
           py::arg("atDepth"),
           "Shallowest lateral position [domain length units] where net_doping "
           "changes sign at the given depth.  Use for vertical (lateral) PN "
           "junctions where P and B are implanted side by side.  Returns inf "
           "if no lateral junction exists at that depth.")
      .def("lateralJunctionPositions",
           &NetDoping<T, D>::lateralJunctionPositions, py::arg("atDepth"),
           "All lateral junction positions at the given depth, sorted "
           "ascending.");
}
