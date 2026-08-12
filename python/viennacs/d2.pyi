"""
2D bindings
"""
from __future__ import annotations
import collections.abc
import typing
import viennacs._core
import viennacs.d2
__all__: list[str] = ['Anneal', 'AtomicLayerProcess', 'DefectDiagnosticsRow', 'DenseCellSet', 'EmbeddedBoundaryPoint', 'Implant', 'ImplantModel', 'MeanFreePath', 'NetDoping', 'Precursor', 'SegmentCells', 'SheetResistance']
class Anneal:
    def __init__(self: viennacs.d2.Anneal) -> None:
        ...
    def addIsothermalStep(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def addRampStep(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat) -> None:
        ...
    def apply(self: viennacs.d2.Anneal) -> None:
        ...
    def applyActivation(self: viennacs.d2.Anneal) -> None:
        """
        Apply only the solid-activation model (C_active = C_SS·C/(C_SS+C)) without running the diffusion solver.  Equivalent to Sentaurus 'diffuse time=0'.  Requires setCellSet(), enableSolidActivation(), and setSolidSolubilityArrhenius() to be configured first.
        """
    def clearDefectDiagnostics(self: viennacs.d2.Anneal) -> None:
        ...
    def clearEquilibriumArrhenius(self: viennacs.d2.Anneal) -> None:
        ...
    def clearSourceField(self: viennacs.d2.Anneal) -> None:
        """
        Remove any previously set concentration source field.
        """
    def clearTemperatureSchedule(self: viennacs.d2.Anneal) -> None:
        ...
    def diffusivity(self: viennacs.d2.Anneal) -> float:
        ...
    def enableDefectClustering(self: viennacs.d2.Anneal, enable: bool = True) -> None:
        ...
    def enableDefectCoupling(self: viennacs.d2.Anneal, enable: bool = True) -> None:
        ...
    def enableDefectEquilibrium(self: viennacs.d2.Anneal, enable: bool = True) -> None:
        ...
    def enableDiagnostics(self: viennacs.d2.Anneal, enable: bool = True) -> None:
        ...
    def enableSolidActivation(self: viennacs.d2.Anneal, enable: bool = True) -> None:
        """
        Enable the solid solubility activation model. When active, writes 'active_concentration' field: C_A+ = C_SS*C/(C_SS+C).
        """
    def getDefectDiagnostics(self: viennacs.d2.Anneal) -> list[viennacs.d2.DefectDiagnosticsRow]:
        ...
    def resetDefectInitialization(self: viennacs.d2.Anneal) -> None:
        ...
    def setActiveLabel(self: viennacs.d2.Anneal, label: str) -> None:
        """
        Set the cell-set field name for the active concentration output (default: 'active_concentration').
        """
    def setArrheniusParameters(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def setBlockingMaterials(self: viennacs.d2.Anneal, arg0: collections.abc.Sequence[typing.SupportsInt]) -> None:
        ...
    def setCellSet(self: viennacs.d2.Anneal, arg0: viennacs.d2.DenseCellSet) -> None:
        ...
    def setClampNonNegative(self: viennacs.d2.Anneal, enable: bool = True) -> None:
        ...
    def setDamageLabels(self: viennacs.d2.Anneal, arg0: str, arg1: str) -> None:
        ...
    def setDefectClusterInitFraction(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat) -> None:
        ...
    def setDefectClusterKinetics(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat) -> None:
        ...
    def setDefectClusterLabel(self: viennacs.d2.Anneal, arg0: str) -> None:
        ...
    def setDefectDiffusivities(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def setDefectEnhancedDiffusion(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def setDefectEquilibrium(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def setDefectEquilibriumArrhenius(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat, arg3: typing.SupportsFloat) -> None:
        ...
    def setDefectLabels(self: viennacs.d2.Anneal, arg0: str, arg1: str) -> None:
        ...
    def setDefectPartition(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def setDefectPartitionFactors(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def setDefectReactionRates(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat) -> None:
        ...
    def setDefectSourceWeights(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def setDiagnosticsMaterialFilter(self: viennacs.d2.Anneal, materialId: typing.SupportsInt = -1) -> None:
        ...
    def setDiffusionCoefficient(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat) -> None:
        ...
    def setDiffusionMaterials(self: viennacs.d2.Anneal, arg0: collections.abc.Sequence[typing.SupportsInt]) -> None:
        ...
    def setDuration(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat) -> None:
        ...
    def setImplicitSolverOptions(self: viennacs.d2.Anneal, arg0: typing.SupportsInt, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat) -> None:
        ...
    def setMaterialLabel(self: viennacs.d2.Anneal, arg0: str) -> None:
        ...
    def setMode(self: viennacs.d2.Anneal, arg0: viennacs._core.DiffusionSolverMode) -> None:
        ...
    def setSolidSolubilityArrhenius(self: viennacs.d2.Anneal, C0: typing.SupportsFloat, Ea_eV: typing.SupportsFloat) -> None:
        """
        Set the solid solubility Arrhenius parameters manually. C0 must be in nm⁻³ (same units as the concentration field).
        """
    def setSourceField(self: viennacs.d2.Anneal, source: collections.abc.Sequence[typing.SupportsFloat]) -> None:
        """
        Per-cell volumetric source term S added to dc/dt = D∇²c + S. Pass an empty list to clear.
        """
    def setSpeciesLabel(self: viennacs.d2.Anneal, arg0: str) -> None:
        ...
    def setStabilityFactor(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat) -> None:
        ...
    def setSurfaceBoundaryCondition(self: viennacs.d2.Anneal, condition: viennacs._core.BoundaryCondition) -> None:
        """
        Set the embedded boundary condition applied at level-set surfaces (default: zero-flux Neumann). Only used when the cell set has embedded boundaries.
        """
    def setTEDFromDamageFactor(self: viennacs.d2.Anneal, damageFactor: typing.SupportsFloat, coefficientScale: typing.SupportsFloat = 0.5, normalization: typing.SupportsFloat = 1e+20) -> None:
        ...
    def setTemperature(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat) -> None:
        ...
    def setTemperatureSchedule(self: viennacs.d2.Anneal, durations: collections.abc.Sequence[typing.SupportsFloat], temperatures: collections.abc.Sequence[typing.SupportsFloat]) -> None:
        """
        Set a temperature schedule from lists of durations (s) and temperatures (K). N temperatures → N isothermal steps; N+1 temperatures → N ramp steps.
        """
    def setTimeStep(self: viennacs.d2.Anneal, arg0: typing.SupportsFloat) -> None:
        ...
class AtomicLayerProcess:
    def __init__(self: viennacs.d2.AtomicLayerProcess, cellSet: viennacs.d2.DenseCellSet, etch: bool = False) -> None:
        ...
    def apply(self: viennacs.d2.AtomicLayerProcess) -> None:
        ...
    @typing.overload
    def setFirstPrecursor(self: viennacs.d2.AtomicLayerProcess, arg0: str, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat, arg3: typing.SupportsFloat, arg4: typing.SupportsFloat, arg5: typing.SupportsFloat) -> None:
        ...
    @typing.overload
    def setFirstPrecursor(self: viennacs.d2.AtomicLayerProcess, arg0: viennacs.d2.Precursor) -> None:
        ...
    def setMaxLambda(self: viennacs.d2.AtomicLayerProcess, arg0: typing.SupportsFloat) -> None:
        ...
    def setMaxTimeStep(self: viennacs.d2.AtomicLayerProcess, arg0: typing.SupportsFloat) -> None:
        ...
    def setPrintInterval(self: viennacs.d2.AtomicLayerProcess, arg0: typing.SupportsFloat) -> None:
        ...
    def setPurgeParameters(self: viennacs.d2.AtomicLayerProcess, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> None:
        ...
    def setReactionOrder(self: viennacs.d2.AtomicLayerProcess, arg0: typing.SupportsFloat) -> None:
        ...
    @typing.overload
    def setSecondPrecursor(self: viennacs.d2.AtomicLayerProcess, arg0: str, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat, arg3: typing.SupportsFloat, arg4: typing.SupportsFloat, arg5: typing.SupportsFloat) -> None:
        ...
    @typing.overload
    def setSecondPrecursor(self: viennacs.d2.AtomicLayerProcess, arg0: viennacs.d2.Precursor) -> None:
        ...
    def setStabilityFactor(self: viennacs.d2.AtomicLayerProcess, arg0: typing.SupportsFloat) -> None:
        ...
class DefectDiagnosticsRow:
    def __init__(self: viennacs.d2.DefectDiagnosticsRow) -> None:
        ...
    @property
    def IV_over_IeqVeq_mean(self) -> float:
        ...
    @IV_over_IeqVeq_mean.setter
    def IV_over_IeqVeq_mean(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def I_max(self) -> float:
        ...
    @I_max.setter
    def I_max(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def I_mean(self) -> float:
        ...
    @I_mean.setter
    def I_mean(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def I_min(self) -> float:
        ...
    @I_min.setter
    def I_min(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def I_over_Ieq_mean(self) -> float:
        ...
    @I_over_Ieq_mean.setter
    def I_over_Ieq_mean(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def Ieq(self) -> float:
        ...
    @Ieq.setter
    def Ieq(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def V_max(self) -> float:
        ...
    @V_max.setter
    def V_max(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def V_mean(self) -> float:
        ...
    @V_mean.setter
    def V_mean(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def V_min(self) -> float:
        ...
    @V_min.setter
    def V_min(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def V_over_Veq_mean(self) -> float:
        ...
    @V_over_Veq_mean.setter
    def V_over_Veq_mean(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def Veq(self) -> float:
        ...
    @Veq.setter
    def Veq(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def step(self) -> int:
        ...
    @step.setter
    def step(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def temperature_K(self) -> float:
        ...
    @temperature_K.setter
    def temperature_K(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def time_s(self) -> float:
        ...
    @time_s.setter
    def time_s(self, arg0: typing.SupportsFloat) -> None:
        ...
class DenseCellSet:
    def __init__(self: viennacs.d2.DenseCellSet) -> None:
        ...
    @typing.overload
    def addFillingFraction(self: viennacs.d2.DenseCellSet, arg0: typing.SupportsInt, arg1: typing.SupportsFloat) -> bool:
        """
        Add to the filling fraction at given cell index.
        """
    @typing.overload
    def addFillingFraction(self: viennacs.d2.DenseCellSet, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"], arg1: typing.SupportsFloat) -> bool:
        """
        Add to the filling fraction for cell which contains given point.
        """
    def addFillingFractionInMaterial(self: viennacs.d2.DenseCellSet, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"], arg1: typing.SupportsFloat, arg2: typing.SupportsInt) -> bool:
        """
        Add to the filling fraction for cell which contains given point only if the cell has the specified material ID.
        """
    def addScalarData(self: viennacs.d2.DenseCellSet, arg0: str, arg1: typing.SupportsFloat) -> None:
        """
        Add a scalar value to be stored and modified in each cell.
        """
    def addVectorData(self: viennacs.d2.DenseCellSet, name: str, initValue: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"] = [0.0, 0.0, 0.0]) -> None:
        """
        Add a vector value (3 components) to be stored and modified in each cell.
        """
    def buildNeighborhood(self: viennacs.d2.DenseCellSet, forceRebuild: bool = False) -> None:
        """
        Generate fast neighbor access for each cell.
        """
    def clear(self: viennacs.d2.DenseCellSet) -> None:
        """
        Clear the filling fractions.
        """
    def embeddedBoundariesEnabled(self: viennacs.d2.DenseCellSet) -> bool:
        """
        Return true if embedded boundary generation was enabled.
        """
    def enableEmbeddedBoundaries(self: viennacs.d2.DenseCellSet, enable: bool = True) -> None:
        """
        Enable sub-grid embedded boundary point generation. Must be called before fromLevelSets().
        """
    def fromLevelSets(self: viennacs.d2.DenseCellSet, levelSets: collections.abc.Sequence[..., ...], materialMap: ... = None, depth: typing.SupportsFloat = 0.0) -> None:
        ...
    def getAverageFillingFraction(self: viennacs.d2.DenseCellSet, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"], arg1: typing.SupportsFloat) -> float:
        """
        Get the average filling at a point in some radius.
        """
    def getBoundingBox(self: viennacs.d2.DenseCellSet) -> typing.Annotated[list[typing.Annotated[list[float], "FixedSize(2)"]], "FixedSize(2)"]:
        ...
    def getCellCenter(self: viennacs.d2.DenseCellSet, arg0: typing.SupportsInt) -> typing.Annotated[list[float], "FixedSize(3)"]:
        """
        Get the center of a cell with given index
        """
    def getCellGrid(self: viennacs.d2.DenseCellSet) -> ...:
        """
        Get the underlying mesh of the cell set.
        """
    def getCellSetPosition(self: viennacs.d2.DenseCellSet) -> bool:
        """
        Get whether the cell set is created below or above the surface.
        """
    def getDepth(self: viennacs.d2.DenseCellSet) -> float:
        """
        Get the depth of the cell set.
        """
    def getElement(self: viennacs.d2.DenseCellSet, arg0: typing.SupportsInt) -> typing.Annotated[list[int], "FixedSize(4)"]:
        """
        Get the element at the given index.
        """
    def getElements(self: viennacs.d2.DenseCellSet) -> list[typing.Annotated[list[int], "FixedSize(4)"]]:
        """
        Get elements (cells). The indicies in the elements correspond to the corner nodes.
        """
    def getEmbeddedBoundaryPointIds(self: viennacs.d2.DenseCellSet, cellIdx: typing.SupportsInt) -> list[int]:
        """
        Return the indices into getEmbeddedBoundaryPoints() for the given cell.
        """
    def getEmbeddedBoundaryPoints(self: viennacs.d2.DenseCellSet) -> list[viennacs.d2.EmbeddedBoundaryPoint]:
        """
        Return the list of all EmbeddedBoundaryPoint objects.
        """
    def getFaceBoundaryDistance(self: viennacs.d2.DenseCellSet, cellIdx: typing.SupportsInt, faceIdx: typing.SupportsInt) -> float:
        """
        Distance from the cell center to the boundary point on face faceIdx. Returns gridDelta/2 when no point exists.
        """
    def getFaceBoundaryPointId(self: viennacs.d2.DenseCellSet, cellIdx: typing.SupportsInt, faceIdx: typing.SupportsInt) -> int:
        """
        Index into getEmbeddedBoundaryPoints() for the boundary point on face faceIdx (axis*2 + (offset>0?1:0)) of cellIdx, or -1 if none.
        """
    def getFillingFraction(self: viennacs.d2.DenseCellSet, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(2)"]) -> float:
        """
        Get the filling fraction of the cell containing the point.
        """
    def getFillingFractions(self: viennacs.d2.DenseCellSet) -> list[float]:
        """
        Get the filling fractions of all cells.
        """
    def getGridDelta(self: viennacs.d2.DenseCellSet) -> float:
        """
        Get the cell size.
        """
    def getIndex(self: viennacs.d2.DenseCellSet, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"]) -> int:
        """
        Get the index of the cell containing the given point.
        """
    def getLevelSets(self: viennacs.d2.DenseCellSet) -> list[..., ...]:
        """
        Get the level sets used to construct the cell set.
        """
    def getMinFaceBoundaryDistance(self: viennacs.d2.DenseCellSet) -> float:
        """
        Minimum face-boundary distance across all cells. Returns gridDelta/2 when no embedded boundaries are present.
        """
    def getNeighbors(self: viennacs.d2.DenseCellSet, arg0: typing.SupportsInt) -> typing.Annotated[list[int], "FixedSize(4)"]:
        """
        Get the neighbor indices for a cell.
        """
    def getNode(self: viennacs.d2.DenseCellSet, arg0: typing.SupportsInt) -> typing.Annotated[list[float], "FixedSize(3)"]:
        """
        Get the node at the given index.
        """
    def getNodes(self: viennacs.d2.DenseCellSet) -> list[typing.Annotated[list[float], "FixedSize(3)"]]:
        """
        Get the nodes of the cell set which correspond to the corner points of the cells.
        """
    def getNumberOfCells(self: viennacs.d2.DenseCellSet) -> int:
        """
        Get the number of cells.
        """
    def getScalarData(self: viennacs.d2.DenseCellSet, arg0: str) -> list[float]:
        """
        Get the data stored at each cell. WARNING: This function only returns a copy of the data
        """
    def getScalarDataLabels(self: viennacs.d2.DenseCellSet) -> list[str]:
        """
        Get the labels of the scalar data stored in the cell set.
        """
    def getSurface(self: viennacs.d2.DenseCellSet) -> ...:
        """
        Get the surface level-set.
        """
    def getVectorData(self: viennacs.d2.DenseCellSet, arg0: str) -> list[typing.Annotated[list[float], "FixedSize(3)"]]:
        """
        Get the vector data stored at each cell. WARNING: This function only returns a copy of the data
        """
    def hasEmbeddedBoundaries(self: viennacs.d2.DenseCellSet) -> bool:
        """
        Return true if the cell set has any embedded boundary points.
        """
    def numEmbeddedBoundaryPoints(self: viennacs.d2.DenseCellSet) -> int:
        """
        Return the total number of embedded boundary points.
        """
    def readCellSetData(self: viennacs.d2.DenseCellSet, arg0: str) -> None:
        """
        Read cell set data from text.
        """
    def setCellSetPosition(self: viennacs.d2.DenseCellSet, arg0: bool) -> None:
        """
        Set whether the cell set should be created below (false) or above (true) the surface.
        """
    def setCoverMaterial(self: viennacs.d2.DenseCellSet, arg0: typing.SupportsInt) -> None:
        """
        Set the material of the cells which are above or below the surface.
        """
    @typing.overload
    def setFillingFraction(self: viennacs.d2.DenseCellSet, arg0: typing.SupportsInt, arg1: typing.SupportsFloat) -> bool:
        """
        Sets the filling fraction at given cell index.
        """
    @typing.overload
    def setFillingFraction(self: viennacs.d2.DenseCellSet, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"], arg1: typing.SupportsFloat) -> bool:
        """
        Sets the filling fraction for cell which contains given point.
        """
    def setPeriodicBoundary(self: viennacs.d2.DenseCellSet, arg0: typing.Annotated[collections.abc.Sequence[bool], "FixedSize(2)"]) -> None:
        """
        Enable periodic boundary conditions in specified dimensions.
        """
    def setScalarData(self: viennacs.d2.DenseCellSet, name: str, newData: collections.abc.Sequence[typing.SupportsFloat]) -> None:
        """
        Overwrite the scalar data associated with 'name' with a new array.
        """
    def setVectorData(self: viennacs.d2.DenseCellSet, name: str, newData: collections.abc.Sequence[typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"]]) -> None:
        """
        Overwrite the vector data associated with 'name' with a new array.
        """
    def updateMaterials(self: viennacs.d2.DenseCellSet) -> None:
        """
        Update the material IDs of the cell set. This function should be called if the level sets, the cell set is made out of, have changed. This does not work if the surface of the volume has changed. In this case, call the function 'updateSurface' first.
        """
    def updateSurface(self: viennacs.d2.DenseCellSet) -> None:
        """
        Updates the surface of the cell set. The new surface should be below the old surface as this function can only remove cells from the cell set.
        """
    def writeCellSetData(self: viennacs.d2.DenseCellSet, arg0: str) -> None:
        """
        Save cell set data in simple text format.
        """
    def writeVTU(self: viennacs.d2.DenseCellSet, arg0: str) -> None:
        """
        Write the cell set as .vtu file
        """
class EmbeddedBoundaryPoint:
    def __init__(self: viennacs.d2.EmbeddedBoundaryPoint) -> None:
        ...
    @property
    def adjacentCell(self) -> int:
        ...
    @adjacentCell.setter
    def adjacentCell(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def axis(self) -> int:
        ...
    @axis.setter
    def axis(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def coordinate(self) -> typing.Annotated[list[float], "FixedSize(3)"]:
        ...
    @coordinate.setter
    def coordinate(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"]) -> None:
        ...
    @property
    def edgeFraction(self) -> float:
        ...
    @edgeFraction.setter
    def edgeFraction(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def levelSetIndex(self) -> int:
        ...
    @levelSetIndex.setter
    def levelSetIndex(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def negativeMaterial(self) -> int:
        ...
    @negativeMaterial.setter
    def negativeMaterial(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def normal(self) -> typing.Annotated[list[float], "FixedSize(3)"]:
        ...
    @normal.setter
    def normal(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"]) -> None:
        ...
    @property
    def positiveMaterial(self) -> int:
        ...
    @positiveMaterial.setter
    def positiveMaterial(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def signedDistance(self) -> float:
        ...
    @signedDistance.setter
    def signedDistance(self, arg0: typing.SupportsFloat) -> None:
        ...
class Implant:
    def __init__(self: viennacs.d2.Implant) -> None:
        ...
    def apply(self: viennacs.d2.Implant) -> None:
        ...
    def enableBeamHits(self: viennacs.d2.Implant, enable: bool = True) -> None:
        ...
    def setBeamHitsLabel(self: viennacs.d2.Implant, arg0: str) -> None:
        ...
    def setCellSet(self: viennacs.d2.Implant, arg0: viennacs.d2.DenseCellSet) -> None:
        ...
    def setConcentrationLabel(self: viennacs.d2.Implant, arg0: str) -> None:
        ...
    def setDamageFactor(self: viennacs.d2.Implant, arg0: typing.SupportsFloat) -> None:
        ...
    def setDamageLabel(self: viennacs.d2.Implant, arg0: str) -> None:
        ...
    def setDamageModel(self: viennacs.d2.Implant, arg0: viennacs.d2.ImplantModel) -> None:
        ...
    def setDose(self: viennacs.d2.Implant, arg0: typing.SupportsFloat) -> None:
        ...
    def setDoseControl(self: viennacs.d2.Implant, arg0: viennacs._core.ImplantDoseControl) -> None:
        ...
    def setImplantAngle(self: viennacs.d2.Implant, arg0: typing.SupportsFloat) -> None:
        ...
    def setImplantModel(self: viennacs.d2.Implant, arg0: viennacs.d2.ImplantModel) -> None:
        ...
    def setLastDamageLabel(self: viennacs.d2.Implant, arg0: str) -> None:
        ...
    def setLengthUnitInCm(self: viennacs.d2.Implant, arg0: typing.SupportsFloat) -> None:
        ...
    def setMaskMaterials(self: viennacs.d2.Implant, materials: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
        Set the material IDs to be treated as mask materials.
        """
    def setOutputConcentrationInCm3(self: viennacs.d2.Implant, enable: bool = True) -> None:
        ...
    def setScreenMaterials(self: viennacs.d2.Implant, materials: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
        Set the material IDs to be treated as screen/cap materials.
        """
    def setVoidMaterial(self: viennacs.d2.Implant, material: typing.SupportsInt) -> None:
        """
        Set a single material ID treated as vacuum/void (beam passes through). Default: 0.
        """
    def setVoidMaterials(self: viennacs.d2.Implant, materials: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
        Set material IDs treated as vacuum/void (beam passes through). Default: {0}.
        """
class ImplantModel:
    def getDepthProfile(self: viennacs.d2.ImplantModel, arg0: typing.SupportsFloat) -> float:
        ...
    def getLateralProfile(self: viennacs.d2.ImplantModel, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> float:
        ...
    def getMaxDepth(self: viennacs.d2.ImplantModel) -> float:
        ...
    def getMaxLateralRange(self: viennacs.d2.ImplantModel) -> float:
        ...
    def getProfile(self: viennacs.d2.ImplantModel, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat) -> float:
        ...
class MeanFreePath:
    def __init__(self: viennacs.d2.MeanFreePath, arg0: viennacs.d2.DenseCellSet) -> None:
        ...
    def apply(self: viennacs.d2.MeanFreePath) -> None:
        ...
    def disableSmoothing(self: viennacs.d2.MeanFreePath) -> None:
        ...
    def enableSmoothing(self: viennacs.d2.MeanFreePath) -> None:
        ...
    def setBulkLambda(self: viennacs.d2.MeanFreePath, arg0: typing.SupportsFloat) -> None:
        ...
    def setMaterial(self: viennacs.d2.MeanFreePath, arg0: typing.SupportsInt) -> None:
        ...
    def setNumRaysPerCell(self: viennacs.d2.MeanFreePath, arg0: typing.SupportsFloat) -> None:
        ...
    def setReflectionLimit(self: viennacs.d2.MeanFreePath, arg0: typing.SupportsInt) -> None:
        ...
    def setRngSeed(self: viennacs.d2.MeanFreePath, arg0: typing.SupportsInt) -> None:
        ...
class NetDoping:
    """
    Compute net doping (Σ donors − Σ acceptors) and extract the
    metallurgical junction depth from a DenseCellSet.
    
    Donor and acceptor labels typically refer to active-concentration
    fields written by the Anneal solid-activation model (e.g. 'P_active',
    'B_active'), but total-concentration fields work too.
    
    Example::
    
      nd = NetDoping()
      nd.setCellSet(domain.getCellSet())
      nd.addDonorLabel('P_active')
      nd.addAcceptorLabel('B_active')
      nd.apply()                      # writes 'net_doping' field
      xj = nd.junctionDepth()         # nm from surface
    """
    def __init__(self: viennacs.d2.NetDoping) -> None:
        ...
    def addAcceptorLabel(self: viennacs.d2.NetDoping, label: str) -> None:
        """
        Append one acceptor (p-type) concentration field name.
        """
    def addDonorLabel(self: viennacs.d2.NetDoping, label: str) -> None:
        """
        Append one donor (n-type) concentration field name.
        """
    def apply(self: viennacs.d2.NetDoping) -> None:
        """
        Compute net_doping = Σ donors − Σ acceptors and write to the output field.
        """
    def junctionCount(self: viennacs.d2.NetDoping) -> int:
        """
        Number of metallurgical junctions (net_doping sign changes).
        """
    def junctionDepth(self: viennacs.d2.NetDoping) -> float:
        """
        Shallowest depth [domain length units] where net_doping changes sign.  Returns inf if no junction exists.
        """
    def junctionDepths(self: viennacs.d2.NetDoping) -> list[float]:
        """
        All junction depths [domain length units], sorted ascending.  Useful for retrograde profiles with multiple sign changes.
        """
    def lateralJunctionPosition(self: viennacs.d2.NetDoping, atDepth: typing.SupportsFloat) -> float:
        """
        Shallowest lateral position [domain length units] where net_doping changes sign at the given depth.  Use for vertical (lateral) PN junctions where P and B are implanted side by side.  Returns inf if no lateral junction exists at that depth.
        """
    def lateralJunctionPositions(self: viennacs.d2.NetDoping, atDepth: typing.SupportsFloat) -> list[float]:
        """
        All lateral junction positions at the given depth, sorted ascending.
        """
    def setAcceptorLabels(self: viennacs.d2.NetDoping, labels: collections.abc.Sequence[str]) -> None:
        """
        Replace the acceptor label list.
        """
    def setCellSet(self: viennacs.d2.NetDoping, cellSet: viennacs.d2.DenseCellSet) -> None:
        """
        Attach the DenseCellSet to analyse.
        """
    def setDepthAxis(self: viennacs.d2.NetDoping, axis: typing.SupportsInt) -> None:
        """
        Cell-centre axis index for depth (default: D−1).
        """
    def setDonorLabels(self: viennacs.d2.NetDoping, labels: collections.abc.Sequence[str]) -> None:
        """
        Replace the donor label list.
        """
    def setOutputLabel(self: viennacs.d2.NetDoping, label: str) -> None:
        """
        Name of the output scalar field (default: 'net_doping').
        """
    def setSurfacePosition(self: viennacs.d2.NetDoping, surfacePosition: typing.SupportsFloat) -> None:
        """
        Wafer-surface coordinate along the depth axis. Depth is computed as surfacePosition minus the cell-centre coordinate.
        """
class Precursor:
    name: str
    def __init__(self: viennacs.d2.Precursor) -> None:
        ...
    @property
    def adsorptionRate(self) -> float:
        ...
    @adsorptionRate.setter
    def adsorptionRate(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def desorptionRate(self) -> float:
        ...
    @desorptionRate.setter
    def desorptionRate(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def duration(self) -> float:
        ...
    @duration.setter
    def duration(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def inFlux(self) -> float:
        ...
    @inFlux.setter
    def inFlux(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def meanThermalVelocity(self) -> float:
        ...
    @meanThermalVelocity.setter
    def meanThermalVelocity(self, arg0: typing.SupportsFloat) -> None:
        ...
class SegmentCells:
    @typing.overload
    def __init__(self: viennacs.d2.SegmentCells, arg0: viennacs.d2.DenseCellSet) -> None:
        ...
    @typing.overload
    def __init__(self: viennacs.d2.SegmentCells, cellSet: viennacs.d2.DenseCellSet, cellTypeString: str = 'CellType', bulkMaterial: typing.SupportsInt = 1) -> None:
        ...
    def apply(self: viennacs.d2.SegmentCells) -> None:
        """
        Segment the cells into surface, material, and gas cells.
        """
    def setBulkMaterial(self: viennacs.d2.SegmentCells, arg0: typing.SupportsInt) -> None:
        """
        Set the bulk material in the segmenter.
        """
    def setCellSet(self: viennacs.d2.SegmentCells, arg0: viennacs.d2.DenseCellSet) -> None:
        """
        Set the cell set in the segmenter.
        """
    def setCellTypeString(self: viennacs.d2.SegmentCells, arg0: str) -> None:
        """
        Set the cell type string in the segmenter.
        """
class SheetResistance:
    """
    Compute sheet resistance (Rsh, Ω/□) from a concentration field stored in
    a DenseCellSet.
    
    Default configuration targets ViennaPS nm-unit domains:
      length unit = 1e-7 (nm → cm),  conc unit = 1e21 (nm⁻³ → cm⁻³),
      depth axis  = D−1  (y for 2-D, z for 3-D),
      surface position = 0  (depth = surface − coordinate).
    
    Example::
    
      sr = SheetResistance()
      sr.setCellSet(domain.getCellSet())
      sr.setConcentrationLabel("P_active")
      rsh = sr.computeElectron()   # Masetti n-type (P in Si)
    """
    def __init__(self: viennacs.d2.SheetResistance) -> None:
        ...
    def computeElectron(self: viennacs.d2.SheetResistance) -> float:
        """
        Rsh [Ω/□] using the Masetti-Severi electron mobility model (n-type, e.g. P-doped Si).
        """
    def computeHole(self: viennacs.d2.SheetResistance) -> float:
        """
        Rsh [Ω/□] using the Masetti-Severi hole mobility model (p-type, e.g. B-doped Si).
        """
    def setCellSet(self: viennacs.d2.SheetResistance, cellSet: viennacs.d2.DenseCellSet) -> None:
        """
        Attach the DenseCellSet to analyse.
        """
    def setConcentrationLabel(self: viennacs.d2.SheetResistance, label: str) -> None:
        """
        Name of the scalar field holding the active concentration (default: 'active_concentration').
        """
    def setConcentrationUnit(self: viennacs.d2.SheetResistance, unit: typing.SupportsFloat) -> None:
        """
        Multiplicative factor to convert the cell-set concentration to cm⁻³ (default: 1e21 for nm⁻³ fields).
        """
    def setDepthAxis(self: viennacs.d2.SheetResistance, axis: typing.SupportsInt) -> None:
        """
        Cell-centre axis index for depth  (default: D−1).
        """
    def setLengthUnit(self: viennacs.d2.SheetResistance, lu_cm: typing.SupportsFloat) -> None:
        """
        Conversion factor from domain length unit to cm (default: 1e-7 for nm domains). Also updates the concentration unit to stay consistent.
        """
    def setSurfacePosition(self: viennacs.d2.SheetResistance, surfacePosition: typing.SupportsFloat) -> None:
        """
        Wafer-surface coordinate along the depth axis. Depth is computed as surfacePosition minus the cell-centre coordinate.
        """
