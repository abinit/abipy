"""
Functions providing access to file data for unit tests and tutorias.
Preferred way to import the module is via the import syntax:

import abipy.data as data
"""
from __future__ import print_function, division

import os

from os.path import join as pj, dirname, abspath

__all__ = [

]

dirpath = dirname(__file__)

_CIF_DIRPATH = abspath( pj(dirname(__file__), "cifs") )

_PSEUDOS_DIRPATH = abspath( pj(dirname(__file__), "pseudos") )

try:
    _DATA_NCFILES = [pj(dirpath, f) for f in os.listdir(dirpath) if f.endswith(".nc")]
except:
    _DATA_NCFILES = []

def ncfile(filename):
    """Return the absolute path of file filename located in data. None if not found."""
    for f in _DATA_NCFILES:
        if os.path.basename(f) == filename:
            return f

    return None

def ncfiles_with_ext(ext):
    """Return a list with the absolute path of the files with extension ext."""
    ncfiles = []
    for filename in _DATA_NCFILES:
        f = filename.rstrip(".nc").rstrip("-etsf")
        if f.endswith("_"+ext):
            ncfiles.append(filename)

    return ncfiles

def ref_file(filename):
    """Returns the absolute path of filename in tests/data directory."""
    return os.path.join(dirpath, filename)


def cif_file(filename):
    """Returns the absolute path of the CIF file in tests/data/cifs."""
    return os.path.join(_CIF_DIRPATH, filename)

def _pseudo(filename):
    """Returns the absolute path of a pseudopotential file in tests/data/pseudos."""
    return os.path.join(_PSEUDOS_DIRPATH, filename)

def pseudos(*filenames):
    """Returns the absolute path of a pseudopotential file in tests/data/pseudos."""
    return [_pseudo(f) for f in filenames]

##########################################################################################

WFK_NCFILES = ncfiles_with_ext("WFK")

DEN_NCFILES = ncfiles_with_ext("DEN")

ALL_NCFILES = WFK_NCFILES + DEN_NCFILES

