from pymatgen.util.io_utils import which
from pymatgen.io.abinitio.eos import EOS
from pymatgen.io.abinitio.wrappers import Mrgscr, Mrgddb, Mrggkk, Anaddb
from pymatgen.io.abinitio import qadapters
from pymatgen.io.abinitio.tasks import * 
from pymatgen.io.abinitio.workflows import *
from pymatgen.io.abinitio.flows import *

from abipy.core import constants
from abipy.core.structure import Structure, StructureModifier
from abipy.htc.input import AbiInput, LdauParams, LexxParams
from abipy.electrons import ElectronDosPlotter, ElectronBandsPlotter, SIGRES_Plotter
from abipy.phonons import PhononBands, PHDOS_Reader, PHDOS_File

FloatWithUnit = constants.FloatWithUnit
ArrayWithUnit = constants.ArrayWithUnit


def _straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def abifile_subclass_from_filename(filename):
    from abipy.iotools.files import AbinitFile, AbinitLogFile, AbinitOutputFile
    from abipy.electrons import SIGRES_File, GSR_File, MDF_File
    from abipy.waves import WFK_File
    #from abipy.phonons import PHDOS_File, PHBST_File

    ext2ncfile = {
        "SIGRES.nc": SIGRES_File,
        "WFK-etsf.nc": WFK_File,
        "MDF.nc" : MDF_File,
        "GSR.nc": GSR_File
        #"PHDOS.nc": PHDOS_File,
        #"PHBST.nc": PHBST_File,
    }

    #if filename.endswith(".abi"):
    #    return AbinitInputFile
                                                                                        
    if filename.endswith(".abo"):
        return AbinitOutputFile
    
    if filename.endswith(".log"):
        return AbinitLogFile

    # CIF files.
    if filename.endswith(".cif"):
        from abipy.core.structure import Structure
        return Structure.from_file(filename)

    ext = filename.split("_")[-1]
    try:
        return ext2ncfile[ext]
    except KeyError:
        raise KeyError("No class has been registered for extension %s" % ext)


def abiopen(filepath):
    """
    Factory function that opens any file supported by abipy.

    Args:
        filepath:
            string with the filename. 
    """
    cls = abifile_subclass_from_filename(filepath)
    return cls.from_file(filepath)


def software_stack(with_wx=True):
    """
    Import all the hard dependencies.
    Returns a dict with the version.
    """
    import numpy, scipy, netCDF4

    d = dict(
        numpy=numpy.version.version,
        scipy=scipy.version.version,
        netCDF4=netCDF4.getlibversion(),
    )

    if with_wx:
        import wx
        d["wx"] = wx=wx.version()

    return d


def abicheck():
    """
    This function tests if the most important ABINIT executables
    can be found in $PATH and whether the python modules needed
    at run-time can be imported.

    Raises:
        RuntimeError if not all the dependencies are fulfilled.
    """
    # executables must be in $PATH. Unfortunately we cannot 
    # test if the binaries work (e.g. shared libraries) 
    # A possible approach would be to execute "exe -h"
    # but supporting argv in Fortran is not trivial.
    executables = [
        "abinit",
        "mrgddb",
        "mrggkk",
        "anaddb",
    ]

    err_lines = []
    app = err_lines.append
    for exe in executables:
        if which(exe) is None:
            app("Cannot found %s in $PATH" % exe)

    try:    
        software_stack(with_wx=False)
    except:
        app(_straceback())

    if err_lines:
        header = "The environment on the local machine is not properly setup\n"
        raise RuntimeError(header + "\n".join(err_lines))

    return 0
            


