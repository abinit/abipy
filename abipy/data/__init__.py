"""
Functions providing access to file data for unit tests and tutorials.
Preferred way to import the module is via the import syntax:

import abipy.data as data
"""
from __future__ import print_function, division

import os

from pymatgen.io.abinitio.pseudos import PseudoParser, PseudoTable
from abipy.data.ucells import structure_from_ucell

__all__ = [
    "cif_file",
    "pseudo",
    "pseudos",
    "ref_file",
    "ref_files",
    "structure_from_ucell",
]

dirpath = os.path.dirname(__file__)

_CIF_DIRPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "cifs"))

_PSEUDOS_DIRPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "pseudos"))


def cif_file(filename):
    """Returns the absolute path of the CIF file in tests/data/cifs."""
    return os.path.join(_CIF_DIRPATH, filename)

pseudo_dir = _PSEUDOS_DIRPATH


def pseudo(filename):
    """Returns a `Pseudo` object."""
    filepath = os.path.join(_PSEUDOS_DIRPATH, filename)
    return PseudoParser().parse(filepath)


def pseudos(*filenames):
    """Returns a PseudoTable constructed from the input filenames  located in tests/data/pseudos."""
    pseudos = map(pseudo, filenames)
    return PseudoTable(pseudos)


def find_ncfiles(top):
    """
    Find all netcdf files starting from the top-level directory top.
    Filenames must be unique. Directories whose start with "tmp_" are
    excluded from the search.

    Returns:
        dictionary with mapping: basename --> absolute path.
    """
    SILENT = 0
    ncfiles = {}
    for dirpath, dirnames, filenames in os.walk(top):

        if "tmp_" in dirpath:
            continue

        for basename in filenames:
            apath = os.path.join(dirpath, basename)
            if basename.endswith(".nc"):

                if basename in ncfiles:
                    err_msg =  "Found duplicated basename %s\n" % basename
                    err_msg += "Stored: %s, new %s\n" % (ncfiles[basename], apath)

                    if not SILENT:
                        SILENT += 1
                        #import warnings
                        #warnings.warn(err_msg)
                        #raise ValueError(err_msg)
                else:
                    ncfiles[basename] = apath 

    return ncfiles

_DATA_NCFILES = find_ncfiles(top=os.path.dirname(__file__))


def ref_file(basename):
    """Returns the absolute path of basename in tests/data directory."""
    if basename in _DATA_NCFILES:
        return _DATA_NCFILES[basename]
    else:
        return os.path.join(dirpath, basename)


def ref_files(*basenames):
    return map(ref_file, basenames)


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

SIGRES_NCFILES = ncfiles_with_ext("SIGRES")

ALL_NCFILES = _DATA_NCFILES.values()
