"""
InterPhon: A Python Package for Ab initio Interface Phonon Calculations within a 3D Electronic Structure Framework.

This package consists of following five sub-packages:

error -> Collection of error modules defined by developer to guide users.
inout -> Collection of inout modules to read and write files in different DFT formats.
core -> Collection of core modules responsible for the central processes.
util -> Collection of util modules to extend the functionality.
analysis -> Collection of analysis modules to characterize phonon properties.
"""

__version__ = "1.4.1"
__all__ = ["error", "inout", "core", "util", "analysis"]
