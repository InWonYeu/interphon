"""
InterPhon error implementations.

This sub-package includes a following module:

error.py -> Collection of error class defined by developer to guide user.
"""

from .error import Insufficient_PBC_Error, Insufficient_DIM_Error, \
    Mismatch_DIM_and_PBC_Error, Mismatch_DIM_post_Error, \
    Invalid_Line_Kpath_Error, Not_Specified_Kpath_Error

__all__ = ["error",
           "Insufficient_PBC_Error",
           "Insufficient_DIM_Error",
           "Mismatch_DIM_and_PBC_Error",
           "Mismatch_DIM_post_Error",
           "Invalid_Line_Kpath_Error",
           "Not_Specified_Kpath_Error"]
