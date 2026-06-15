"""
ViennaCS cell-set framework and solver kernels.
"""
from __future__ import annotations
import typing
from viennacs import d2
from viennacs import d3
__all__: list[str] = ['BoundaryCondition', 'BoundaryConditionType', 'DiffusionSolverMode', 'Dirichlet', 'ImplantDoseControl', 'Neumann', 'Robin', 'd2', 'd3', 'setNumThreads', 'version']
class BoundaryCondition:
    type: BoundaryConditionType
    @staticmethod
    def dirichlet(boundaryValue: typing.SupportsFloat) -> BoundaryCondition:
        """
        Create a Dirichlet condition with the given surface value.
        """
    @staticmethod
    def neumann(outwardFlux: typing.SupportsFloat = 0.0) -> BoundaryCondition:
        """
        Create a Neumann condition with the given outward flux (default 0 = zero flux).
        """
    @staticmethod
    def robin(transfer: typing.SupportsFloat, exteriorValue: typing.SupportsFloat) -> BoundaryCondition:
        """
        Create a Robin condition: flux = transfer*(c - exteriorValue).
        """
    def __init__(self) -> None:
        ...
    @property
    def transferCoefficient(self) -> float:
        ...
    @transferCoefficient.setter
    def transferCoefficient(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def value(self) -> float:
        ...
    @value.setter
    def value(self, arg0: typing.SupportsFloat) -> None:
        ...
class BoundaryConditionType:
    """
    Members:
    
      Neumann
    
      Dirichlet
    
      Robin
    """
    Dirichlet: typing.ClassVar[BoundaryConditionType]  # value = <BoundaryConditionType.Dirichlet: 1>
    Neumann: typing.ClassVar[BoundaryConditionType]  # value = <BoundaryConditionType.Neumann: 0>
    Robin: typing.ClassVar[BoundaryConditionType]  # value = <BoundaryConditionType.Robin: 2>
    __members__: typing.ClassVar[dict[str, BoundaryConditionType]]  # value = {'Neumann': <BoundaryConditionType.Neumann: 0>, 'Dirichlet': <BoundaryConditionType.Dirichlet: 1>, 'Robin': <BoundaryConditionType.Robin: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class DiffusionSolverMode:
    """
    Members:

      Explicit

      GaussSeidel
    """
    Explicit: typing.ClassVar[DiffusionSolverMode]  # value = <DiffusionSolverMode.Explicit: 0>
    GaussSeidel: typing.ClassVar[DiffusionSolverMode]  # value = <DiffusionSolverMode.GaussSeidel: 1>
    __members__: typing.ClassVar[dict[str, DiffusionSolverMode]]  # value = {'Explicit': <DiffusionSolverMode.Explicit: 0>, 'GaussSeidel': <DiffusionSolverMode.GaussSeidel: 1>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class ImplantDoseControl:
    """
    Members:
    
      Off
    
      WaferDose
    
      BeamDose
    """
    BeamDose: typing.ClassVar[ImplantDoseControl]  # value = <ImplantDoseControl.BeamDose: 2>
    Off: typing.ClassVar[ImplantDoseControl]  # value = <ImplantDoseControl.Off: 0>
    WaferDose: typing.ClassVar[ImplantDoseControl]  # value = <ImplantDoseControl.WaferDose: 1>
    __members__: typing.ClassVar[dict[str, ImplantDoseControl]]  # value = {'Off': <ImplantDoseControl.Off: 0>, 'WaferDose': <ImplantDoseControl.WaferDose: 1>, 'BeamDose': <ImplantDoseControl.BeamDose: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
def setNumThreads(arg0: typing.SupportsInt) -> None:
    ...
Dirichlet: BoundaryConditionType  # value = <BoundaryConditionType.Dirichlet: 1>
Neumann: BoundaryConditionType  # value = <BoundaryConditionType.Neumann: 0>
Robin: BoundaryConditionType  # value = <BoundaryConditionType.Robin: 2>
__version__: str = '2.0.1'
version: str = '2.0.1'
