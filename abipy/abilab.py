"""
This module gathers the most important classes and helper functions used for scripting.
"""
from __future__ import print_function, division, unicode_literals

from monty.os.path import which
from pymatgen.core.units import *
from pymatgen.io.abinitio.eos import EOS
from pymatgen.io.abinitio.pseudos import PseudoTable
from pymatgen.io.abinitio.wrappers import Mrgscr, Mrgddb, Mrggkk
#from pymatgen.io.abinitio.tasks import (TaskManager, ScfTask, NscfTask, RelaxTask, DDK_Task,
#    PhononTask, G_Task, HaydockBseTask, OpticTask, AnaddbTask)
#from pymatgen.io.abinitio.workflows import (Workflow, IterativeWorkflow, BandStructureWorkflow,
#    RelaxWorkflow, DeltaFactorWorkflow, G0W0_Workflow, SigmaConvWorkflow, BSEMDF_Workflow,
#    PhononWorkflow)
from pymatgen.io.abinitio.tasks import *
from pymatgen.io.abinitio.workflows import *
from pymatgen.io.abinitio.flows import *
from pymatgen.io.abinitio.launcher import PyFlowScheduler

from abipy.tools.prettytable import PrettyTable
from abipy.core.structure import Structure, StructureModifier
from abipy.htc.input import AbiInput, LdauParams, LexxParams, input_gen, AnaddbInput
from abipy.electrons import ElectronDosPlotter, ElectronBandsPlotter
from abipy.electrons.gsr import GSR_File, GsrRobot, SigresRobot 
from abipy.electrons.gw import SIGRES_File, SIGRES_Plotter 
from abipy.electrons.bse import MDF_File
from abipy.electrons.scissors import ScissorsBuilder
from abipy.phonons import PHBST_File, PhononBands, PHDOS_File, PHDOS_Reader
from abipy.core.mixins import AbinitFile, AbinitLogFile, AbinitOutputFile
from abipy.waves import WFK_File

# Tools for unit conversion
#from abipy.core import constants
#FloatWithUnit = constants.FloatWithUnit
#ArrayWithUnit = constants.ArrayWithUnit

def _straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def abifile_subclass_from_filename(filename):
    """Returns the appropriate class associated to the given filename."""
    ext2ncfile = {
        "SIGRES.nc": SIGRES_File,
        "WFK-etsf.nc": WFK_File,
        "MDF.nc": MDF_File,
        "GSR.nc": GSR_File,
        "PHBST.nc": PHBST_File,
    }

    #if filename.endswith(".abi"):
    #    return AbinitInputFile
                                                                                        
    if filename.endswith(".abo"): return AbinitOutputFile
    if filename.endswith(".log"): return AbinitLogFile

    # CIF files.
    if filename.endswith(".cif"):
        return Structure

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


def software_stack():
    """
    Import all the hard dependencies.
    Returns a dict with the version.
    """
    # Mandatory
    import numpy, scipy

    d = dict(
        numpy=numpy.version.version,
        scipy=scipy.version.version,
    )

    # Optional but strongly suggested.
    try:
        import netCDF4, matplotlib
        d.update(dict(
            netCDF4=netCDF4.getlibversion(),
            matplotlib="Version: %s, backend: %s" % (matplotlib.__version__, matplotlib.get_backend()),
            ))
    except ImportError:
        pass

    # Optional (GUIs).
    try:
        import wx
        d["wx"] = wx.version()
    except ImportError:
        pass

    return d


def abicheck():
    """
    This function tests if the most important ABINIT executables
    can be found in $PATH and whether the python modules needed
    at run-time can be imported.

    Raises:
        RuntimeError if not all the dependencies are fulfilled.
    """
    import os
    # Executables must be in $PATH. Unfortunately we cannot
    # test the version of the binaries.
    # A possible approach would be to execute "exe -v"
    # but supporting argv in Fortran is not trivial.
    # Dynamic linking is tested by calling `ldd exe`
    executables = [
        "abinit",
        "mrgddb",
        "mrggkk",
        "anaddb",
    ]

    has_ldd = which("ldd") is not None

    err_lines = []
    app = err_lines.append
    for exe in executables:
        exe_path = which(exe)

        if exe_path is None:
            app("Cannot find %s in $PATH" % exe)
        else:
            if has_ldd and os.system("ldd %s > /dev/null " % exe_path) != 0:
                app("Missing shared library dependencies for %s" % exe)

    try:    
        software_stack()
    except:
        app(_straceback())

    if err_lines:
        header = "The environment on the local machine is not properly setup\n"
        raise RuntimeError(header + "\n".join(err_lines))

    return 0


def flow_main(main):
    """
    This decorator is used to decorate main functions producing `AbinitFlows`.
    It adds the initialization of the logger and an argument parser that allows one to select 
    the loglevel, the workdir of the flow as well as the YAML file with the parameters of the `TaskManager`.
    The main function shall have the signature:

        main(options)

    where options in the container with the command line options generated by `ArgumentParser`.

    Args:
        main:
            main function.
    """
    from functools import wraps

    @wraps(main)
    def wrapper(*args, **kwargs):
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('--loglevel', default="ERROR", type=str,
                            help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

        parser.add_argument("-w", '--workdir', default="", type=str, help="Working directory of the flow.")

        parser.add_argument("-m", '--manager', default="", type=str,
                            help="YAML file with the parameters of the task manager")

        parser.add_argument("--prof", action="store_true", default=False, help="Profile code wth cProfile ")

        options = parser.parse_args()

        # loglevel is bound to the string value obtained from the command line argument. 
        # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
        import logging
        numeric_level = getattr(logging, options.loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % options.loglevel)
        logging.basicConfig(level=numeric_level)

        if options.prof:
            import pstats, cProfile
            cProfile.runctx("main(options)", globals(), locals(), "Profile.prof")
            s = pstats.Stats("Profile.prof")
            s.strip_dirs().sort_stats("time").print_stats()
            return 0
        else:
            return main(options)

    return wrapper
