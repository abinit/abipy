"""
Functions providing access to file data for unit tests and tutorials.
Preferred way to import the module is via the import syntax:

import abipy.data as data
"""
from __future__ import print_function, division, unicode_literals

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
    pseudos = list(map(pseudo, filenames))
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
                        import warnings
                        warnings.warn(err_msg)
                        raise ValueError(err_msg)
                        SILENT += 1

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
    """List with the absolute path of basenames in tests/data directory."""
    return list(map(ref_file, basenames))


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



class AbinitFilesGenerator(object):
    """This class generates the output files used in the unit tests and in the examples."""
    # Subclasses must define the following class attributes:
    # List of pseudos in (basenames in abipy/data/pseudos
    #pseudos = ["14si.pspnc"]

    # Mapping old_name --> new_name for the output files that must be preserved.
    #files_to_save = {
    #    "out_DS1_DEN-etsf.nc": "si_DEN-etsf.nc",
    #    "out_DS2_GSR.nc": "si_nscf_GSR.nc",
    #}

    def __init__(self, workdir=".", finalize=True, verbose=1):
        """
        Args:
            workdir:
                Working directory.
            finalize:
                True if self.finalize is called
            verbose:
                Verbosity level.
        """
        if not hasattr(self, "pseudos") or not hasattr(self, "files_to_save"):
            raise ValueError("pseudos and files_to_save are not defined")

        from monty.os.path import which
        if which("abinit") is None:
            raise RuntimeError("Cannot find abinit in $PATH")

        # Absolute path for the pseudopotentials.
        self.pseudos = [p.filepath for p in pseudos(*self.pseudos)]
        self.filesfile = "\n".join(["run.abi", "run.abo", "in", "out","tmp"] + self.pseudos)
        self.workdir = os.path.abspath(workdir)
        self.finalize, self.verbose = finalize, verbose

        self.files_to_keep = \
            ["run.abi", "run.abo", os.path.basename(__file__)] + list(self.files_to_save.keys())

    def __str__(self):
        lines = []
        app = lines.append
        app("workdir = %s" % self.workdir)
        app("pseudos = %s" % self.pseudos)

        return "\n".join(lines)

    def run(self):
        """Run Abinit and rename output files. Return 0 if success"""
        os.chdir(self.workdir)
        process = self._run()
        process.wait()

        if process.returncode != 0: 
            print("returncode == %s" % process.returncode)
            print(process.stderr.readlines())
            return process.returncode

        if self.finalize:
            self._finalize()

        return 0

    def _run(self):
        from subprocess import Popen, PIPE
        with open(os.path.join(self.workdir, "run.files"), "w") as fh:
            fh.write(self.filesfile)

        cmd = "abinit < run.files > run.log"
        return Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE, cwd=self.workdir)

    def _finalize(self):
        all_files = os.listdir(self.workdir)

        # Remove files
        garbage = [f for f in all_files if f not in self.files_to_keep]
        for f in garbage: 
            if f.endswith(".py"): continue
            if self.verbose: print("Will remove file %s" % f)
            os.remove(f)

        # Rename files.
        for old, new in self.files_to_save.items():
            if self.verbose: print("Will rename %s --> %s" % (old, new))
            os.rename(old, new)
