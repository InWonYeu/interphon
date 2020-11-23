"""
InterPhon: A python package to calculate 2-D interface phonon within 3D electronic structure framework.

This package includes following five sub-packages:

error -> Collection of error modules defined by developer to guide user.
inout -> Collection of parser modules to read and write files in different DFT formats.
core -> Collection of core modules responsible for control tower of processes.
util -> Collection of util modules to extend the utility.
analysis -> Collection of analysis modules to characterize the phonon.
"""

__version__ = "2020.11.04"
__all__ = ["error", "inout", "core", "util", "analysis"]
