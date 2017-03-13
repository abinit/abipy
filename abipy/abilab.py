"""
This module gathers the most important classes and helper functions used for scripting.
"""
from __future__ import print_function, division, unicode_literals

import os
import collections

from monty.os.path import which
from pymatgen.core.units import *
from pymatgen.io.abinit.eos import EOS
from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable
from pymatgen.io.abinit.wrappers import Mrgscr, Mrgddb, Mrggkk
from pymatgen.io.abinit.tasks import *
from pymatgen.io.abinit.works import *
from pymatgen.io.abinit.flows import (Flow, G0W0WithQptdmFlow, bandstructure_flow,
    g0w0_flow, phonon_flow, phonon_conv_flow, NonLinearCoeffFlow)
# Need new version of pymatgen.
try:
    from pymatgen.io.abinit.flows import PhononFlow
except ImportError:
    pass

from pymatgen.io.abinit.launcher import PyFlowScheduler #, BatchLauncher
from pymatgen.io.abinit.pseudos import Pseudo

from abipy.core.release import __version__, min_abinit_version
from abipy.core.structure import Lattice, Structure, StructureModifier, frames_from_structures
from abipy.core.mixins import AbinitLogFile, AbinitOutputFile, OutNcFile
from abipy.core.kpoints import set_atol_kdiff
from abipy.htc.input import AbiInput, LdauParams, LexxParams, input_gen
from abipy.iotools import Visualizer
from abipy.iotools.cube import CubeFile
from abipy.abio.timer import AbinitTimerParser
from abipy.abio.robots import Robot, GsrRobot, SigresRobot, MdfRobot, DdbRobot, abirobot
from abipy.abio.inputs import AbinitInput, MultiDataset, AnaddbInput, OpticInput
from abipy.abio.abivars import AbinitInputFile
from abipy.abio.factories import *
from abipy.electrons.ebands import ElectronBands, ElectronDosPlotter, ElectronBandsPlotter, frame_from_ebands
from abipy.electrons.gsr import GsrFile
from abipy.electrons.psps import PspsFile
from abipy.electrons.gw import SigresFile, SigresPlotter
from abipy.electrons.bse import MdfFile
from abipy.electrons.scissors import ScissorsBuilder
from abipy.electrons.scr import ScrFile
from abipy.electrons.ncdenpot import DensityNcFile
from abipy.electrons.fatbands import FatBandsFile
from abipy.dfpt.phonons import (PhbstFile, PhononBands, PhdosFile, PhdosReader,
                                phbands_gridplot)
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.anaddbnc import AnaddbNcFile
from abipy.dynamics.hist import HistFile
from abipy.waves import WfkFile

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

# Abinit text files. Use OrderedDict for nice output in show_abiopen_exc2class.
ext2file = collections.OrderedDict([
    (".abi", AbinitInputFile),
    (".in", AbinitInputFile),
    (".abo", AbinitOutputFile),
    (".out", AbinitOutputFile),
    (".log", AbinitLogFile),
    (".cif", Structure),
    ("POSCAR", Structure),
    (".cssr", Structure),
    (".cube", CubeFile),
    ("anaddb.nc", AnaddbNcFile),
    (".psp8", Pseudo),
    (".xml", Pseudo),
])

# Abinit files require a special treatment.
abiext2ncfile = collections.OrderedDict([
    ("GSR.nc", GsrFile),
    ("DEN.nc", DensityNcFile),
    ("OUT.nc", OutNcFile),
    ("WFK.nc", WfkFile),
    ("HIST.nc", HistFile),
    ("PSPS.nc", PspsFile),
    ("DDB", DdbFile),
    ("PHBST.nc", PhbstFile),
    ("PHDOS.nc", PhdosFile),
    ("SCR.nc", ScrFile),
    ("SIGRES.nc", SigresFile),
    ("MDF.nc", MdfFile),
    ("FATBANDS.nc", FatBandsFile),
])


def abiopen_ext2class_table():
    """
    Print the association table between file extensions and File classes.
    """
    from itertools import chain
    from tabulate import tabulate
    table = []

    for ext, cls in chain(ext2file.items(), abiext2ncfile.items()):
        table.append((ext, str(cls)))

    return tabulate(table, headers=["Extension", "Class"])


def abifile_subclass_from_filename(filename):
    """Returns the appropriate class associated to the given filename."""
    for ext, cls in ext2file.items():
        if filename.endswith(ext): return cls

    ext = filename.split("_")[-1]
    try:
        return abiext2ncfile[ext]
    except KeyError:
        for ext, cls in abiext2ncfile.items():
            if filename.endswith(ext): return cls

    msg = ("No class has been registered for file:\n\t%s\n\nFile extensions supported:\n%s" %
        (filename, abiopen_ext2class_table()))
    raise ValueError(msg)


def dir2abifiles(top, recurse=True):
    """
    Analyze the filesystem starting from directory `top` and
    return an ordered dictionary mapping the directory name to the list
    of files supported by `abiopen` contained within that directory.
    If not `recurse`, children directories are not analyzed.
    """
    dl = collections.defaultdict(list)

    if recurse:
        for dirpath, dirnames, filenames in os.walk(top):
            for f in filenames:
                path = os.path.join(dirpath, f)
                if not isabifile(path): continue
                dl[dirpath].append(path)
    else:
        for f in os.listdir(top):
            path = os.path.join(top, f)
            if not isabifile(path): continue
            dl[dirpath].append(path)

    return collections.OrderedDict([(k, dl[k]) for k in sorted(dl.keys())])


def isabifile(filepath):
    """
    Return True if `filepath` can be opened with `abiopen`.
    """
    try:
        abifile_subclass_from_filename(filepath)
        return True
    except ValueError:
        return False


def abiopen(filepath):
    """
    Factory function that opens any file supported by abipy.
    File type is detected from the extension

    Args:
        filepath: string with the filename.
    """
    if os.path.basename(filepath) == "__AbinitFlow__.pickle":
        return Flow.pickle_load(filepath)

    # Handle old output files produced by Abinit.
    import re
    outnum = re.compile(".+\.out[\d]+")
    abonum = re.compile(".+\.abo[\d]+")
    if outnum.match(filepath) or abonum.match(filepath):
        return AbinitOutputFile.from_file(filepath)

    cls = abifile_subclass_from_filename(filepath)
    return cls.from_file(filepath)


def print_frame(frame, title=None):
    """
    Print entire pandas DataFrame.

    Args:
        frame: pandas DataFrame.
        title: Optional string to print as initial title.
    """
    if title is not None: print(title)
    import pandas as pd
    with pd.option_context('display.max_rows', len(frame),
                           'display.max_columns', len(list(frame.keys()))):
        print(frame)
    print()


def display_structure(obj, **kwargs):
    """
    Use Jsmol to display a structure in the jupyter notebook.
    Requires `nbjsmol` notebook extension installed on the local machine.
    Install it with `pip install nbjsmol`. See also https://github.com/gmatteo/nbjsmol.

    Args:
        obj: Structure object or file with a structure or python object with a `structure` attribute.
        kwargs: Keyword arguments passed to `nbjsmol_display`
    """
    try:
        from nbjsmol import nbjsmol_display
    except ImportError as exc:
        raise ImportError(str(exc) +
                          "\ndisplay structure requires nbjsmol package\n."
                          "Install it with `pip install nbjsmol.`\n"
                          "See also https://github.com/gmatteo/nbjsmol.")

    # Cast to structure, get string with cif data and pass it to nbjsmol.
    structure = Structure.as_structure(obj)
    return nbjsmol_display(structure.to(fmt="cif"), ext=".cif", **kwargs)


def software_stack():
    """
    Import all the hard dependencies. Returns ordered dict: package --> string with version info.
    """
    # Mandatory
    import numpy, scipy, netCDF4, pymatgen, apscheduler, pydispatch, yaml

    d = collections.OrderedDict([
        ("numpy", numpy.version.version),
        ("scipy", scipy.version.version),
        ("netCDF4", netCDF4.__version__),
        ("apscheduler", apscheduler.version),
        ("pydispatch", pydispatch.__version__),
        ("yaml", yaml.__version__),
        ("pymatgen", pymatgen.__version__),
    ])

    # Optional but strongly suggested.
    try:
        import matplotlib
        d["matplotlib"] = "%s (backend: %s)" % (matplotlib.__version__, matplotlib.get_backend())
    except ImportError:
        pass

    return d


def abicheck(verbose=0):
    """
    This function tests if the most important ABINIT executables
    can be found in $PATH and whether the python modules needed
    at run-time can be imported. Return string with error messages, empty if success.
    """
    from monty.termcolor import cprint
    err_lines = []
    app = err_lines.append

    try:
        manager = TaskManager.from_user_config()
    except Exception:
        manager = None
        app(_straceback())

    # Get info on the Abinit build.
    from abipy.core.testing import cmp_version
    if manager is not None:
        cprint("AbiPy Manager:\n%s\n" % str(manager), color="green")
        build = AbinitBuild(manager=manager)
        if not build.has_netcdf: app("Abinit executable does not support netcdf")
        cprint("Abinitbuild:\n%s" % str(build), color="magenta")
        if verbose: cprint(str(build.info), color="magenta")
        print()
        if not cmp_version(build.version, min_abinit_version, op=">="):
            app("Abipy requires Abinit version >= %s but got %s" % (min_abinit_version, build.version))

    # Get info on the scheduler.
    from pymatgen.io.abinit.launcher import PyFlowScheduler
    try:
        scheduler = PyFlowScheduler.from_user_config()
        cprint("Abipy Scheduler:\n%s\n" % str(scheduler), color="yellow")
    except Exception as exc:
        app(_straceback())

    from tabulate import tabulate
    try:
        d = software_stack()
        cprint("Installed packages:", color="blue")
        cprint(tabulate(list(d.items()), headers=["Package", "Version"]), color="blue")
        print()
    except ImportError:
        app(_straceback())

    return "\n".join(err_lines)


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

        parser.add_argument("-r", "--remove", default=False, action="store_true", help="Remove old flow workdir")

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


def abipy_logo1():
    """http://www.text-image.com/convert/pic2ascii.cgi"""
    return """\

                 `:-                                                               -:`
         --`  .+/`                              `                                  `/+.  .-.
   `.  :+.   /s-                   `yy         .yo                                   -s/   :+. .`
 ./.  +o`   /s/           `-::-`   `yy.-::-`   `:-    .:::-`   -:`     .:`            /s/   :s- ./.
.o.  /o:   .oo.         .oyo++syo. `yyys++oys. -ys  -syo++sy+` sy-     +y:            .oo-   oo` `o.
++   oo.   /oo          yy-    -yy `yy:    .yy`-ys .ys`    /yo sy-     +y:             oo/   /o:  ++
+/   oo`   /oo         `yy.    .yy` yy.    `yy`-ys :ys     :yo oy/     oy:             +o/   :o:  /o
-/   :+.   -++`         -sy+::+yyy` .sy+::+yy- -ys :yys/::oys. `oyo::/syy:            `++-   /+.  /:
 --  `//    /+-           -/++/-//    -/++/-   `+: :yo:/++/.     .:++/:oy:            -+/   `+-  --
  `.`  -:    :/`                                   :yo                 +y:           `/:`  `:. `.`
        `..   .:.                                   .`                 `.           .:.  `..
                ...                                                               ...
"""

def abipy_logo2():
    """http://www.text-image.com/convert/pic2ascii.cgi"""
    return """\
MMMMMMMMMMMMMMMMNhdMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMdhmMMMMMMMMMMMMMMM
MMMMMMMMMddNMMmoyNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNyomMMmhmMMMMMMMM
MMMmmMMhomMMMy/hMMMMMMMMMMMMMMMMMMMN::MMMMMMMMMm:oMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMd+yMMMhomMmmMMM
MmsmMMs+NMMMy+yMMMMMMMMMMMNhyyhmMMMN::mhyyhmMMMNhdMMMMmhyydNMMMdyNMMMMMmyNMMMMMMMMMMMMy+yMMMh+dMmsmM
m+mMMy+hMMMd++mMMMMMMMMMm+:+ss+:+mMN:::/ss+:+mMd:/MMd/:+so/:oNM/:dMMMMMo:yMMMMMMMMMMMMm++dMMMo+NMN+m
osMMMo+mMMMy+oMMMMMMMMMM::dMMMMd:/MN::hMMMMd::Nd:/Mm:/NMMMMy:oM/:dMMMMMo:yMMMMMMMMMMMMMo+yMMMy+hMMso
oyMMMooNMMMyooMMMMMMMMMN::mMMMMm::NM::mMMMMN::Nd:/Mh:/MMMMMh:+Mo:yMMMMM+:yMMMMMMMMMMMMMooyMMMyohMMyo
dyMMMysmMMMdooNMMMMMMMMMd/:oyys:::NMd/:oyys::dMd:/Mh::/shyo:/mMN+:+yhs/::yMMMMMMMMMMMMNoodMMMysmMMyh
MddMMNyhMMMMysdMMMMMMMMMMMdyooydssMMMMdysoydMMMNsyMh:+hsosydMMMMMmysosho:yMMMMMMMMMMMMdsyMMMNydMMddM
MMNmNMMddMMMMhyNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMh:oMMMMMMMMMMMMMMMMMo:yMMMMMMMMMMMNyhMMMNhmMNmNMM
MMMMMMMMNmmMMMmhmMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMmNMMMMMMMMMMMMMMMMMNdMMMMMMMMMMMmhmMMNmmMMMMMMMM
MMMMMMMMMMMMMMMMmmNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNmmMMMMMMMMMMMMMMM
"""


def abipy_logo3():
    """http://www.text-image.com/convert/pic2ascii.cgi"""
    return """\
             `-.                                                  `--`
      -:. `//`              `/.       ::                           `+: `-:`
 --``+:  `o+        `.--.`  -y/.--.`  ::   `---`  `-     ..         `o+  `o-`--
:/  /o   /o-       -oo/:+s/ -yy+:/os- ss .oo/:/s+`:y.    yo          :o:  :o. /:
o-  o/   oo`       ss    /y..y/    ss ss +y`   -y::y-    yo          -o/  .o- -o
:-  //   /+.       :s+--/sy- +s/--+s: ss oyo:-:so``os:-:oyo          -+:  -+. -/
 -` `/.  `/:        `-:::-:`  `-::-`  -- oy-:::.    .:::-yo          //`  :- `-
   `  ..` `:-                            :+              /:         --` `-` `
            `.`                                                   ..`
"""
