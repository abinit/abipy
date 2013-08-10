"""
Functions providing access to file data for unit tests and tutorias.
Preferred way to import the module is via the import syntax:

import abipy.data as data
"""
from __future__ import print_function, division

import os

from os.path import join as pj

__all__ = [

]

dirpath = os.path.dirname(__file__)

_CIF_DIRPATH = os.path.abspath(pj(os.path.dirname(__file__), "cifs"))

_PSEUDOS_DIRPATH = os.path.abspath(pj(os.path.dirname(__file__), "pseudos"))

def cif_file(filename):
    """Returns the absolute path of the CIF file in tests/data/cifs."""
    return os.path.join(_CIF_DIRPATH, filename)


def pseudos(*filenames):
    """Returns the absolute path of a pseudopotential file in tests/data/pseudos."""
    return [os.path.join(_PSEUDOS_DIRPATH, f) for f in filenames]


def find_ncfiles(top):
    """
    Find all netcdf files starting from the top-level directory top.
    Filenames must be unique.

    Returns:
        dictionary with mapping: basename --> absolute path.
    """
    ncfiles = {}
    for dirpath, dirnames, filenames in os.walk(top):
        for basename in filenames:
            apath = os.path.join(dirpath, basename)
            if basename.endswith(".nc"):
                if basename in ncfiles:
                    raise ValueError("Found duplicated basename %s" % basename)
                ncfiles[basename] = apath 

    return ncfiles

_DATA_NCFILES = find_ncfiles(top=os.path.dirname(__file__))

def ref_file(basename):
    """Returns the absolute path of basename in tests/data directory."""
    if basename in _DATA_NCFILES:
        return _DATA_NCFILES[basename]
    else:
        return os.path.join(dirpath, basename)


def ncfiles_with_ext(ext):
    """Return a list with the absolute path of the files with extension ext."""
    ncfiles = []
    for basename, path in _DATA_NCFILES.items():
        f = basename.rstrip(".nc").rstrip("-etsf")
        if f.endswith("_"+ext):
            ncfiles.append(path)

    return ncfiles

WFK_NCFILES = ncfiles_with_ext("WFK")

DEN_NCFILES = ncfiles_with_ext("DEN")

GSR_NCFILES = ncfiles_with_ext("GSR")

ALL_NCFILES = _DATA_NCFILES.values()
