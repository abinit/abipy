"""
This module gathers the most important classes and helper functions used for scripting.
"""
from __future__ import print_function, division, unicode_literals

from monty.os.path import which
from pymatgen.core.units import *
from pymatgen.io.abinit.eos import EOS
from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable
from pymatgen.io.abinit.wrappers import Mrgscr, Mrgddb, Mrggkk
from pymatgen.io.abinit.tasks import *
from pymatgen.io.abinit.works import *
from pymatgen.io.abinit.flows import (Flow, G0W0WithQptdmFlow, bandstructure_flow, 
    g0w0_flow, phonon_flow, phonon_conv_flow)
# Need new version of pymatgen.
try:
    from pymatgen.io.abinit.flows import PhononFlow
except ImportError:
    pass

from pymatgen.io.abinit.launcher import PyFlowScheduler, BatchLauncher

from abipy.core.structure import Lattice, Structure, StructureModifier
from abipy.htc.input import AbiInput, LdauParams, LexxParams, input_gen
from abipy.abio.robots import GsrRobot, SigresRobot, MdfRobot, DdbRobot, abirobot
from abipy.abio.inputs import AbinitInput, MultiDataset, AnaddbInput, OpticInput
from abipy.abio.factories import *
from abipy.electrons import ElectronDosPlotter, ElectronBandsPlotter, SigresPlotter
from abipy.electrons.gsr import GsrFile
from abipy.electrons.psps import PspsFile
from abipy.electrons.gw import SigresFile, SigresPlotter 
from abipy.electrons.bse import MdfFile
from abipy.electrons.scissors import ScissorsBuilder
from abipy.electrons.scr import ScrFile
from abipy.dfpt import PhbstFile, PhononBands, PhdosFile, PhdosReader
from abipy.dfpt.ddb import DdbFile
from abipy.dynamics.hist import HistFile
from abipy.core.mixins import AbinitInputFile, AbinitLogFile, AbinitOutputFile
from abipy.waves import WfkFile
from abipy.iotools import Visualizer

# Tools for unit conversion
import pymatgen.core.units as units
FloatWithUnit = units.FloatWithUnit
ArrayWithUnit = units.ArrayWithUnit

# Documentation.
from abipy.abio.abivars_db import get_abinit_variables, abinit_help, docvar

# Utils for notebooks.
from abipy.tools.notebooks import mpld3_enable_notebook


def _straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def abifile_subclass_from_filename(filename):
    """Returns the appropriate class associated to the given filename."""
    ext2ncfile = {
        "SIGRES.nc": SigresFile,
        "WFK-etsf.nc": WfkFile,
        "MDF.nc": MdfFile,
        "GSR.nc": GsrFile,
        "SCR-etsf.nc": ScrFile,
        "PHBST.nc": PhbstFile,
        "PHDOS.nc": PhdosFile,
        "DDB": DdbFile,
        "HIST.nc": HistFile,
        "PSPS.nc": PspsFile,
    }

    # Abinit text files.
    if filename.endswith(".abi"): return AbinitInputFile
    if filename.endswith(".abo"): return AbinitOutputFile
    if filename.endswith(".log"): return AbinitLogFile

    # CIF files.
    if filename.endswith(".cif"): return Structure

    ext = filename.split("_")[-1]
    try:
        return ext2ncfile[ext]
    except KeyError:
        #raise KeyError("No class has been registered for extension %s" % ext)
        for ext, cls in ext2ncfile.items():
            if filename.endswith(ext): return cls

        raise ValueError("No class has been registered for filename %s" % filename)


def abiopen(filepath):
    """
    Factory function that opens any file supported by abipy.
    File type is detected from the extension

    Args:
        filepath: string with the filename. 
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
    # Executables must be in $PATH. Unfortunately we cannot test the version of the binaries.
    # A possible approach would be to execute "exe -v" but supporting argv in Fortran is not trivial.
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
    This decorator is used to decorate main functions producing `Flows`.
    It adds the initialization of the logger and an argument parser that allows one to select 
    the loglevel, the workdir of the flow as well as the YAML file with the parameters of the `TaskManager`.
    The main function shall have the signature:

        main(options)

    where options in the container with the command line options generated by `ArgumentParser`.

    Args:
        main: main function.
    """
    from functools import wraps

    @wraps(main)
    def wrapper(*args, **kwargs):
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('--loglevel', default="ERROR", type=str,
                            help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

        parser.add_argument("-w", '--workdir', default="", type=str, help="Working directory of the flow.")

        parser.add_argument("-m", '--manager', default=None, 
                            help="YAML file with the parameters of the task manager. " 
                                 "Default None i.e. the manager is read from standard locations: "
                                 "working directory first then ~/.abinit/abipy/manager.yml.")

        parser.add_argument("-s", '--scheduler', action="store_true", default=False, 
                            help="Run the flow with the scheduler")

        parser.add_argument("-b", '--batch', action="store_true", default=False, 
                            help="Run the flow in batch mode")

        #parser.add_argument("-r", '--remove', action="store_true", default=False, help="Run the flow with the scheduler")

        parser.add_argument("--prof", action="store_true", default=False, help="Profile code wth cProfile ")

        options = parser.parse_args()

        # loglevel is bound to the string value obtained from the command line argument. 
        # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
        import logging
        numeric_level = getattr(logging, options.loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % options.loglevel)
        logging.basicConfig(level=numeric_level)

        # Istantiate the manager.
        options.manager = TaskManager.as_manager(options.manager)

        def execute():
            """This is the function that performs the work depending on options."""
            flow = main(options)

            if options.scheduler:
                flow.rmtree()
                return flow.make_scheduler().start()

            elif options.batch:
                flow.rmtree()
                flow.build_and_pickle_dump()
                return flow.batch()

            return 0

        if options.prof:
            # Profile execute
            import pstats, cProfile
            cProfile.runctx("execute()", globals(), locals(), "Profile.prof")
            s = pstats.Stats("Profile.prof")
            s.strip_dirs().sort_stats("time").print_stats()
            return 0
        else:
            return execute()

    return wrapper
