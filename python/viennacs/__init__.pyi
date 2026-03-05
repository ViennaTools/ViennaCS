"""

ViennaCS
========

ViennaCS is a header-only C++ cell set library, which adds the possibility
of using volumetric representations on top of existing level-set functionalities
for surfaces. Combined with ray tracing techniques, this enables the simulation
of particle scattering and ion implantation."
"""
from __future__ import annotations
import sys as _sys
from viennacs._core import setNumThreads
from viennacs.d2 import AtomicLayerProcess
from viennacs.d2 import DenseCellSet
from viennacs.d2 import MeanFreePath
from viennacs.d2 import Precursor
from viennacs.d2 import SegmentCells
from . import _core
from . import d2
from . import d3
__all__: list[str] = ['AtomicLayerProcess', 'DenseCellSet', 'MeanFreePath', 'PROXY_DIM', 'Precursor', 'SegmentCells', 'd2', 'd3', 'getDimension', 'setDimension', 'setNumThreads']
def __dir__():
    ...
def __getattr__(name):
    ...
def _windows_dll_path():
    ...
def getDimension() -> int:
    """
    Get the current dimension of the simulation.
    
        Returns
        -------
        int
            The currently set dimension (2 or 3).
        
    """
def setDimension(d: int):
    """
    Set the dimension of the simulation (2 or 3).
    
        Parameters
        ----------
        d: int
            Dimension of the simulation (2 or 3).
        
    """
PROXY_DIM: int = 2
__version__: str = 'VIENNACS_VERSION'
_C = _core
