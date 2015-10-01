"""
Functions providing access to file data for unit tests and tutorials.
Preferred way to import the module is via the import syntax:

import abipy.abidata as abidata
"""
from __future__ import print_function, division, unicode_literals

import os

from pymatgen.io.abinit.pseudos import PseudoParser, PseudoTable
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

_VARIABLES_DIRPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "variables"))

_MPDATA_DIRPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "mpdata"))


_SCRIPTS_DIRPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "examples"))
_SCRIPTS = None


def pyscript(basename):
    """Return the absolute name of one of the scripts in the `examples` directory from its basename."""

    global _SCRIPTS 
    if _SCRIPTS is None:
        # Build mapping basename --> path.
        from monty.os.path import find_exts
        pypaths = find_exts(_SCRIPTS_DIRPATH, ".py", exclude_dirs="_*|.*")
        _SCRIPTS = {}
        for p in pypaths:
            k = os.path.basename(p)
            assert k not in _SCRIPTS
            _SCRIPTS[k] = p

    return _SCRIPTS[basename]


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


def var_file(filename):
    """Returns a yml file located in data/variables."""
    return os.path.join(_VARIABLES_DIRPATH, filename)


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
                        #raise ValueError(err_msg)
                        SILENT += 1

                else:
                    ncfiles[basename] = apath 

    return ncfiles

_DATA_NCFILES = find_ncfiles(top=os.path.join(os.path.dirname(__file__), "refs"))


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


def mp_structures():
    """ Returns a dictionary containg the structures stored in mpdata/mp_structures. """
    import json
    from monty.json import MontyDecoder

    with open(os.path.join(_MPDATA_DIRPATH, 'mp_structures.json'), 'r') as f:
        return json.load(f, cls=MontyDecoder)


WFK_NCFILES = ncfiles_with_ext("WFK")

DEN_NCFILES = ncfiles_with_ext("DEN")

GSR_NCFILES = ncfiles_with_ext("GSR")

SIGRES_NCFILES = ncfiles_with_ext("SIGRES")

ALL_NCFILES = _DATA_NCFILES.values()


class FilesGenerator(object):
    """This class generates the output files used in the unit tests and in the examples."""

    def __init__(self, **kwargs):
        """
        Args:
            workdir: Working directory.
            finalize: True if self.finalize is called
            verbose: Verbosity level.
        """
        if not hasattr(self, "files_to_save"):
            raise ValueError("files_to_save are not defined")

        from monty.os.path import which
        if which(self.executable) is None:
            raise RuntimeError("Cannot find %s in $PATH" % self.executable)

        self.workdir = os.path.abspath(kwargs.pop("workdir", "."))
        self.finalize = kwargs.pop("finalize", True)
        self.verbose = kwargs.pop("verbose", 1)

        self.files_to_keep = set([os.path.basename(__file__), "run.abi", "run.abo"] + 
                list(self.files_to_save.keys()))

    def __str__(self):
        lines = []
        app = lines.append
        app("%s: workdir:%s" % (self.__class__.__name__, self.workdir))

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
            fh.write(self.make_filesfile_str())

        cmd = self.executable + " < run.files > run.log"
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


class AbinitFilesGenerator(FilesGenerator):
    # Subclasses must define the following class attributes:
    # List of pseudos in (basenames in abipy/data/pseudos
    #pseudos = ["14si.pspnc"]

    # Mapping old_name --> new_name for the output files that must be preserved.
    #files_to_save = {
    #    "out_DS1_DEN-etsf.nc": "si_DEN-etsf.nc",
    #    "out_DS2_GSR.nc": "si_nscf_GSR.nc",
    #}
    executable = "abinit"

    def __init__(self, **kwargs):
        super(AbinitFilesGenerator, self).__init__(**kwargs)

        # Add Absolute paths for the pseudopotentials.
        self.pseudos = [p.filepath for p in pseudos(*self.pseudos)]

    def make_filesfile_str(self):
        s = "\n".join(["run.abi", "run.abo", "in", "out","tmp"] + self.pseudos)
        return s


class AnaddbFilesGenerator(FilesGenerator):
    # Subclasses must define the following class attributes:

    # 1) Mapping old_name --> new_name for the output files that must be preserved.
    #files_to_save = {
    #    "trf2_5.out_PHBST.nc": "trf2_5.out_PHBST.nc",
    #    "trf2_5.out_PHDOS.nc": "trf2_5.out_PHDOS.nc",
    #}

    # 2) Input DDB (mandatory)
    in_ddb = None

    # 3) output DDB (optional)
    out_ddb = "dummy_out_ddb"

    # 3) Input GKK (optional)
    in_gkk = "dummy_in_gkk"

    # 4) base name for elphon output files (optional)
    elph_basename = "dummy_elph_basename"

    # 5) file containing ddk filenames for elphon/transport
    in_ddk = "dummy_in_ddk"

    executable = "anaddb"

    def __init__(self, **kwargs):
        super(AnaddbFilesGenerator, self).__init__(**kwargs)

        if self.in_ddb is None:
            raise ValueError("in_ddb must be specified")

        self.files_to_keep.update([
            self.in_ddb,
            self.out_ddb,
            self.in_gkk,
            self.elph_basename,
            self.in_ddk,
        ])

    def make_filesfile_str(self):
        # 1) formatted input file
        # 2) formatted output file e.g. t13.out
        # 3) input derivative database e.g. t13.ddb.in
        # 4) output molecular dynamics e.g. t13.md
        # 5) input elphon matrix elements  (GKK file) :
        # 6) base name for elphon output files e.g. t13
        # 7) file containing ddk filenames for elphon/transport
        return "\n".join([
            "run.abi", 
            "run.abo", 
            self.in_ddb,
            self.out_ddb,
            self.in_gkk,
            self.elph_basename,
            self.in_ddk,
        ])


if __name__ == "__main__":
    print(pyscript("plot_spectral_functions.py"))
