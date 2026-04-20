"""
3D bindings
"""
from __future__ import annotations
import collections.abc
import typing
__all__: list[str] = ['AtomicLayerProcess', 'DenseCellSet', 'MeanFreePath', 'Precursor', 'SegmentCells']
class AtomicLayerProcess:
    def __init__(self, cellSet: DenseCellSet, etch: bool = False) -> None:
        ...
    def apply(self) -> None:
        ...
    @typing.overload
    def setFirstPrecursor(self, arg0: str, arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: typing.SupportsFloat | typing.SupportsIndex, arg3: typing.SupportsFloat | typing.SupportsIndex, arg4: typing.SupportsFloat | typing.SupportsIndex, arg5: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def setFirstPrecursor(self, arg0: Precursor) -> None:
        ...
    def setMaxLambda(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def setMaxTimeStep(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def setPrintInterval(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def setPurgeParameters(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def setReactionOrder(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def setSecondPrecursor(self, arg0: str, arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: typing.SupportsFloat | typing.SupportsIndex, arg3: typing.SupportsFloat | typing.SupportsIndex, arg4: typing.SupportsFloat | typing.SupportsIndex, arg5: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def setSecondPrecursor(self, arg0: Precursor) -> None:
        ...
    def setStabilityFactor(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class DenseCellSet:
    def __init__(self) -> None:
        ...
    @typing.overload
    def addFillingFraction(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> bool:
        """
        Add to the filling fraction at given cell index.
        """
    @typing.overload
    def addFillingFraction(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"], arg1: typing.SupportsFloat | typing.SupportsIndex) -> bool:
        """
        Add to the filling fraction for cell which contains given point.
        """
    def addFillingFractionInMaterial(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"], arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Add to the filling fraction for cell which contains given point only if the cell has the specified material ID.
        """
    def addScalarData(self, arg0: str, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Add a scalar value to be stored and modified in each cell.
        """
    def addVectorData(self, name: str, initValue: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"] = [0.0, 0.0, 0.0]) -> None:
        """
        Add a vector value (3 components) to be stored and modified in each cell.
        """
    def buildNeighborhood(self, forceRebuild: bool = False) -> None:
        """
        Generate fast neighbor access for each cell.
        """
    def clear(self) -> None:
        """
        Clear the filling fractions.
        """
    def fromLevelSets(self, levelSets: collections.abc.Sequence[..., ...], materialMap: ... = None, depth: typing.SupportsFloat | typing.SupportsIndex = 0.0) -> None:
        ...
    def getAverageFillingFraction(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"], arg1: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
        Get the average filling at a point in some radius.
        """
    def getBoundingBox(self) -> typing.Annotated[list[typing.Annotated[list[float], "FixedSize(3)"]], "FixedSize(2)"]:
        ...
    def getCellCenter(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[list[float], "FixedSize(3)"]:
        """
        Get the center of a cell with given index
        """
    def getCellGrid(self) -> ...:
        """
        Get the underlying mesh of the cell set.
        """
    def getDepth(self) -> float:
        """
        Get the depth of the cell set.
        """
    def getElement(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[list[int], "FixedSize(8)"]:
        """
        Get the element at the given index.
        """
    def getElements(self) -> list[typing.Annotated[list[int], "FixedSize(8)"]]:
        """
        Get elements (cells). The indicies in the elements correspond to the corner nodes.
        """
    def getFillingFraction(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"]) -> float:
        """
        Get the filling fraction of the cell containing the point.
        """
    def getFillingFractions(self) -> list[float]:
        """
        Get the filling fractions of all cells.
        """
    def getGridDelta(self) -> float:
        """
        Get the cell size.
        """
    def getIndex(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"]) -> int:
        """
        Get the index of the cell containing the given point.
        """
    def getNeighbors(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[list[int], "FixedSize(6)"]:
        """
        Get the neighbor indices for a cell.
        """
    def getNode(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[list[float], "FixedSize(3)"]:
        """
        Get the node at the given index.
        """
    def getNodes(self) -> list[typing.Annotated[list[float], "FixedSize(3)"]]:
        """
        Get the nodes of the cell set which correspond to the corner points of the cells.
        """
    def getNumberOfCells(self) -> int:
        """
        Get the number of cells.
        """
    def getScalarData(self, arg0: str) -> list[float]:
        """
        Get the data stored at each cell. WARNING: This function only returns a copy of the data
        """
    def getScalarDataLabels(self) -> list[str]:
        """
        Get the labels of the scalar data stored in the cell set.
        """
    def getSurface(self) -> ...:
        """
        Get the surface level-set.
        """
    def getVectorData(self, arg0: str) -> list[typing.Annotated[list[float], "FixedSize(3)"]]:
        """
        Get the vector data stored at each cell. WARNING: This function only returns a copy of the data
        """
    def readCellSetData(self, arg0: str) -> None:
        """
        Read cell set data from text.
        """
    def setCellSetPosition(self, arg0: bool) -> None:
        """
        Set whether the cell set should be created below (false) or above (true) the surface.
        """
    def setCoverMaterial(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Set the material of the cells which are above or below the surface.
        """
    @typing.overload
    def setFillingFraction(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> bool:
        """
        Sets the filling fraction at given cell index.
        """
    @typing.overload
    def setFillingFraction(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"], arg1: typing.SupportsFloat | typing.SupportsIndex) -> bool:
        """
        Sets the filling fraction for cell which contains given point.
        """
    def setPeriodicBoundary(self, arg0: typing.Annotated[collections.abc.Sequence[bool], "FixedSize(3)"]) -> None:
        """
        Enable periodic boundary conditions in specified dimensions.
        """
    def setScalarData(self, name: str, newData: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        """
        Overwrite the scalar data associated with 'name' with a new array.
        """
    def setVectorData(self, name: str, newData: collections.abc.Sequence[typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"]]) -> None:
        """
        Overwrite the vector data associated with 'name' with a new array.
        """
    def updateMaterials(self) -> None:
        """
        Update the material IDs of the cell set. This function should be called if the level sets, the cell set is made out of, have changed. This does not work if the surface of the volume has changed. In this case, call the function 'updateSurface' first.
        """
    def updateSurface(self) -> None:
        """
        Updates the surface of the cell set. The new surface should be below the old surface as this function can only remove cells from the cell set.
        """
    def writeCellSetData(self, arg0: str) -> None:
        """
        Save cell set data in simple text format.
        """
    def writeVTU(self, arg0: str) -> None:
        """
        Write the cell set as .vtu file
        """
class MeanFreePath:
    def __init__(self, arg0: DenseCellSet) -> None:
        ...
    def apply(self) -> None:
        ...
    def disableSmoothing(self) -> None:
        ...
    def enableSmoothing(self) -> None:
        ...
    def setBulkLambda(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def setMaterial(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def setNumRaysPerCell(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def setReflectionLimit(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def setRngSeed(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class Precursor:
    name: str
    def __init__(self) -> None:
        ...
    @property
    def adsorptionRate(self) -> float:
        ...
    @adsorptionRate.setter
    def adsorptionRate(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def desorptionRate(self) -> float:
        ...
    @desorptionRate.setter
    def desorptionRate(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def duration(self) -> float:
        ...
    @duration.setter
    def duration(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def inFlux(self) -> float:
        ...
    @inFlux.setter
    def inFlux(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def meanThermalVelocity(self) -> float:
        ...
    @meanThermalVelocity.setter
    def meanThermalVelocity(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class SegmentCells:
    @typing.overload
    def __init__(self, arg0: DenseCellSet) -> None:
        ...
    @typing.overload
    def __init__(self, cellSet: DenseCellSet, cellTypeString: str = 'CellType', bulkMaterial: typing.SupportsInt | typing.SupportsIndex = 1) -> None:
        ...
    def apply(self) -> None:
        """
        Segment the cells into surface, material, and gas cells.
        """
    def setBulkMaterial(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Set the bulk material in the segmenter.
        """
    def setCellSet(self, arg0: DenseCellSet) -> None:
        """
        Set the cell set in the segmenter.
        """
    def setCellTypeString(self, arg0: str) -> None:
        """
        Set the cell type string in the segmenter.
        """
