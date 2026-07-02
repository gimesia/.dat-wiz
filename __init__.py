"""Standalone k-space viewer + NUFFT reconstruction tool.

Inspect raw Siemens ``.dat`` k-space, browse ICE dimensions, and reconstruct a
selected slice with NUFFT using the trajectory from a companion Pulseq ``.seq``.

Run with:  ``python -m src.kspace_viewer.main --dat <.dat> --seq <.seq>``
"""

from .loader import LoadedData, load_dat
from .trajectory import Trajectory, extract_trajectory
from .recon import ReconEngine

__all__ = [
    "LoadedData",
    "load_dat",
    "Trajectory",
    "extract_trajectory",
    "ReconEngine",
]
