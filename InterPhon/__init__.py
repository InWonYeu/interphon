"""
InterPhon: A Python Package for Ab initio Interface Phonon Calculations within a 3D Electronic Structure Framework.

This package includes following five sub-packages:

error -> Collection of error modules defined by developer to guide users.
inout -> Collection of parser modules to read and write files in different DFT formats.
core -> Collection of core modules responsible for the control tower of processes.
util -> Collection of util modules to extend the functionality.
analysis -> Collection of analysis modules to characterize phonons.
"""

__version__ = "0.2.1"
__all__ = ["error", "inout", "core", "util", "analysis"]
