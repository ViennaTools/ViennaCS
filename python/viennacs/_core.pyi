"""
ViennaCS is a header-only C++ cell set library, which adds the possibility of using volumetric representations on top of existing level-set functionalities for surfaces. Combined with ray tracing techniques, this enables the simulation of particle scattering and ion implantation.
"""
from __future__ import annotations
import typing
from viennacs import d2
from viennacs import d3
__all__: list[str] = ['d2', 'd3', 'setNumThreads']
def setNumThreads(arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
    ...
__version__: str = 'VIENNACS_VERSION'
