"""
This module gathers the most important classes and helper functions used for scripting.
"""
import sys
import os
import collections

from itertools import chain
from tabulate import tabulate

####################
### Monty import ###
####################
from monty.os.path import which
from monty.termcolor import cprint

#######################
### Pymatgen import ###
#######################
# Tools for unit conversion
import pymatgen.core.units as units
FloatWithUnit = units.FloatWithUnit
ArrayWithUnit = units.ArrayWithUnit

####################
### Abipy import ###
####################
from abipy.flowtk import Pseudo, PseudoTable, Mrgscr, Mrgddb, Flow, Work, TaskManager, AbinitBuild, flow_main
from abipy.core.release import __version__, min_abinit_version
from abipy.core.globals import enable_notebook, in_notebook, disable_notebook
from abipy.core import restapi
from abipy.core.structure import (Lattice, Structure, StructureModifier, dataframes_from_structures,
  mp_match_structure, mp_search, cod_search)
from abipy.core.mixins import CubeFile
from abipy.core.func1d import Function1D
from abipy.core.kpoints import set_atol_kdiff
from abipy.abio.robots import Robot
from abipy.abio.inputs import AbinitInput, MultiDataset, AnaddbInput, OpticInput
from abipy.abio.abivars import AbinitInputFile
from abipy.abio.outputs import AbinitLogFile, AbinitOutputFile, OutNcFile, AboRobot
from abipy.tools.printing import print_dataframe
from abipy.tools.notebooks import print_source, print_doc
from abipy.tools.plotting import get_ax_fig_plt, get_axarray_fig_plt, get_ax3d_fig_plt
from abipy.abio.factories import *
from abipy.electrons.ebands import (ElectronBands, ElectronBandsPlotter, ElectronDos, ElectronDosPlotter,
    dataframe_from_ebands, EdosFile)
from abipy.electrons.gsr import GsrFile, GsrRobot
from abipy.electrons.eskw import EskwFile
from abipy.electrons.psps import PspsFile
from abipy.electrons.gw import SigresFile, SigresRobot
from abipy.electrons.bse import MdfFile, MdfRobot
from abipy.electrons.scissors import ScissorsBuilder
from abipy.electrons.scr import ScrFile
from abipy.electrons.denpot import (DensityNcFile, VhartreeNcFile, VxcNcFile, VhxcNcFile, PotNcFile,
    DensityFortranFile, Cut3dDenPotNcFile)
from abipy.electrons.fatbands import FatBandsFile
from abipy.electrons.optic import OpticNcFile, OpticRobot
from abipy.electrons.fold2bloch import Fold2BlochNcfile
from abipy.dfpt.phonons import (PhbstFile, PhbstRobot, PhononBands, PhononBandsPlotter, PhdosFile, PhononDosPlotter,
    PhdosReader, phbands_gridplot)
from abipy.dfpt.ddb import DdbFile, DdbRobot
from abipy.dfpt.anaddbnc import AnaddbNcFile, AnaddbNcRobot
from abipy.dfpt.gruneisen import GrunsNcFile
from abipy.dynamics.hist import HistFile, HistRobot
from abipy.waves import WfkFile
from abipy.eph.a2f import A2fFile, A2fRobot
from abipy.eph.sigeph import SigEPhFile, SigEPhRobot
from abipy.eph.eph_plotter import EphPlotter
from abipy.eph.v1sym import V1symFile
from abipy.eph.gkq import GkqFile, GkqRobot
from abipy.eph.v1qnu import V1qnuFile
from abipy.eph.v1qavg import V1qAvgFile
from abipy.eph.rta import RtaFile, RtaRobot
from abipy.eph.transportfile import TransportFile
from abipy.wannier90 import WoutFile, AbiwanFile, AbiwanRobot
from abipy.electrons.lobster import CoxpFile, ICoxpFile, LobsterDoscarFile, LobsterInput, LobsterAnalyzer
#from abipy.electrons.abitk import ZinvConvFile, TetraTestFile

# Abinit Documentation.
from abipy.abio.abivars_db import get_abinit_variables, abinit_help, docvar


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
    (".abivars", Structure),
    (".ucell", Structure),
    ("POSCAR", Structure),
    (".cssr", Structure),
    (".cube", CubeFile),
    ("anaddb.nc", AnaddbNcFile),
    ("DEN", DensityFortranFile),
    (".psp8", Pseudo),
    (".pspnc", Pseudo),
    (".fhi", Pseudo),
    ("JTH.xml", Pseudo),
    (".wout", WoutFile),
    # Lobster files.
    ("COHPCAR.lobster", CoxpFile),
    ("COOPCAR.lobster", CoxpFile),
    ("ICOHPLIST.lobster", ICoxpFile),
    ("DOSCAR.lobster", LobsterDoscarFile),
    #("ZINVCONV.nc", ZinvConvFile),
    #("TETRATEST.nc", TetraTestFile),
    ("EDOS", EdosFile),
])

# Abinit files require a special treatment.
abiext2ncfile = collections.OrderedDict([
    ("GSR.nc", GsrFile),
    ("ESKW.nc", EskwFile),
    ("DEN.nc", DensityNcFile),
    ("OUT.nc", OutNcFile),
    ("VHA.nc", VhartreeNcFile),
    ("VXC.nc", VxcNcFile),
    ("VHXC.nc", VhxcNcFile),
    ("POT.nc", PotNcFile),
    ("WFK.nc", WfkFile),
    ("HIST.nc", HistFile),
    ("PSPS.nc", PspsFile),
    ("DDB", DdbFile),
    ("PHBST.nc", PhbstFile),
    ("PHDOS.nc", PhdosFile),
    ("SCR.nc", ScrFile),
    ("SIGRES.nc", SigresFile),
    ("GRUNS.nc", GrunsNcFile),
    ("MDF.nc", MdfFile),
    ("FATBANDS.nc", FatBandsFile),
    ("FOLD2BLOCH.nc", Fold2BlochNcfile),
    ("CUT3DDENPOT.nc", Cut3dDenPotNcFile),
    ("OPTIC.nc", OpticNcFile),
    ("A2F.nc", A2fFile),
    ("SIGEPH.nc", SigEPhFile),
    ("TRANSPORT.nc",TransportFile),
    ("RTA.nc",RtaFile),
    ("V1SYM.nc", V1symFile),
    ("GKQ.nc", GkqFile),
    ("V1QNU.nc", V1qnuFile),
    ("V1QAVG.nc", V1qAvgFile),
    ("ABIWAN.nc", AbiwanFile),
])


def abiopen_ext2class_table():
    """
    Print the association table between file extensions and File classes.
    """
    table = []

    for ext, cls in chain(ext2file.items(), abiext2ncfile.items()):
        table.append((ext, str(cls)))

    return tabulate(table, headers=["Extension", "Class"])


def extcls_supporting_panel(as_table=True, **tabulate_kwargs):
    """
    Return list of (ext, cls) tuples where ext is the file extension supporting the `get_panel` method
    i.e files that can be integrated with the panel dashboard and `cls` the python class associated to it.
    """
    items = []

    for ext, cls in chain(ext2file.items(), abiext2ncfile.items()):
        if hasattr(cls, "get_panel"): items.append((ext, str(cls)))

    if not as_table: return items

    return tabulate(items, headers=["Extension", "Class"], **tabulate_kwargs)


def abipanel(**kwargs):
    """
    Activate panel extensions used by AbiPy. Return panel module.
    """
    from abipy.panels.core import abipanel
    return abipanel(**kwargs)


def abifile_subclass_from_filename(filename):
    """
    Returns the appropriate class associated to the given filename.
    """
    if os.path.basename(filename) == Flow.PICKLE_FNAME:
        return Flow

    from abipy.tools.text import rreplace
    for ext, cls in ext2file.items():
        # This to support gzipped files.
        if filename.endswith(".gz"): filename = rreplace(filename, ".gz", "", occurrence=1)
        if filename.endswith(ext): return cls

    ext = filename.split("_")[-1]
    try:
        return abiext2ncfile[ext]
    except KeyError:
        for ext, cls in abiext2ncfile.items():
            if filename.endswith(ext): return cls

    msg = ("No class has been registered for file:\n\t%s\n\nFile extensions supported:\n\n%s" %
        (filename, abiopen_ext2class_table()))
    raise ValueError(msg)


def dir2abifiles(top, recurse=True):
    """
    Analyze the filesystem starting from directory `top` and
    return an ordered dictionary mapping the directory name to the list
    of files supported by ``abiopen`` contained within that directory.
    If not ``recurse``, children directories are not analyzed.
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
            dl[top].append(path)

    return collections.OrderedDict([(k, dl[k]) for k in sorted(dl.keys())])


def isabifile(filepath):
    """
    Return True if `filepath` can be opened with ``abiopen``.
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
    # Handle ~ in filepath.
    filepath = os.path.expanduser(filepath)

    # Handle zipped files by creating temporary file with correct extension.
    root, ext = os.path.splitext(filepath)
    if ext in (".bz2", ".gz", ".z"):
        from monty.io import zopen
        with zopen(filepath, "rt") as f:
            import tempfile
            _, tmp_path = tempfile.mkstemp(suffix=os.path.basename(root), text=True)
            cprint("Creating temporary file: %s" % tmp_path, "yellow")
            with open(tmp_path, "wt") as t:
                t.write(f.read())
            filepath = tmp_path

    if os.path.basename(filepath) == "__AbinitFlow__.pickle":
        return Flow.pickle_load(filepath)

    # Handle old output files produced by Abinit.
    import re
    outnum = re.compile(r".+\.out[\d]+")
    abonum = re.compile(r".+\.abo[\d]+")
    if outnum.match(filepath) or abonum.match(filepath):
        return AbinitOutputFile.from_file(filepath)

    if os.path.basename(filepath) == "log":
        # Assume Abinit log file.
        return AbinitLogFile.from_file(filepath)

    cls = abifile_subclass_from_filename(filepath)
    return cls.from_file(filepath)


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


def mjson_load(filepath, **kwargs):
    """
    Read JSON file in MSONable format with MontyDecoder. Return dict with python objects.
    """
    import json
    from monty.json import MontyDecoder
    with open(filepath, "rt") as fh:
        return json.load(fh, cls=MontyDecoder, **kwargs)


def mjson_loads(string, **kwargs):
    """
    Read JSON string in MSONable format with MontyDecoder. Return dict with python objects.
    """
    import json
    from monty.json import MontyDecoder
    return json.loads(string, cls=MontyDecoder, **kwargs)


def mjson_write(d, filepath, **kwargs):
    """
    Write dictionary d to filepath in JSON format using MontyDecoder
    """
    import json
    from monty.json import MontyEncoder
    with open(filepath, "wt") as fh:
        json.dump(d, fh, cls=MontyEncoder, **kwargs)


def software_stack(as_dataframe=False):
    """
    Import all the hard dependencies and some optional packages.
    Returns ordered dict: package --> string with version info or pandas dataframe if as_dataframe.
    """
    import platform
    system, node, release, version, machine, processor = platform.uname()
    # These packages are required
    import numpy, scipy, netCDF4, pymatgen, apscheduler, pydispatch, yaml, plotly

    from importlib import import_module
    def get_version(pkg_name):
        """Return version of package from string."""
        try:
            mod = import_module(pkg_name)
            return mod.__version__
        except Exception:
            return None

    try:
        from pymatgen.core import __version__ as pmg_version
        #from pymatgen.settings import __version__ as pmg_version
    except AttributeError:
        pmg_version = pymatgen.__version__

    d = collections.OrderedDict([
        ("system", system),
        ("python_version", platform.python_version()),
        ("numpy", numpy.version.version),
        ("scipy", scipy.version.version),
        ("netCDF4", netCDF4.__version__),
        ("apscheduler", apscheduler.version),
        ("pydispatch", pydispatch.__version__),
        ("yaml", yaml.__version__),
        ("boken", get_version("bokeh")),
        ("panel", get_version("panel")),
        ("plotly", get_version("plotly")),
        ("ase", get_version("ase")),
        ("phonopy", get_version("phonopy")),
        ("monty", get_version("monty")),
        ("pymatgen", pmg_version),
        ("abipy", __version__),
    ])

    if not as_dataframe: return d
    import pandas as pd
    return pd.Series(data=d, name="version").to_frame().rename_axis("Package")


def abicheck(verbose=0):
    """
    This function tests if the most important ABINIT executables
    can be found in $PATH and whether the python modules needed
    at run-time can be imported. Return string with error messages, empty if success.
    """

    err_lines = []
    app = err_lines.append

    try:
        manager = TaskManager.from_user_config()
    except Exception:
        manager = None
        app(_straceback())

    # Get info on the Abinit build.
    from abipy.core.testing import cmp_version
    from abipy.flowtk import PyFlowScheduler

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
    try:
        scheduler = PyFlowScheduler.from_user_config()
        cprint("Abipy Scheduler:\n%s\n" % str(scheduler), color="yellow")
    except Exception as exc:
        app(_straceback())

    try:
        d = software_stack()
        cprint("Installed packages:", color="blue")
        cprint(tabulate(list(d.items()), headers=["Package", "Version"]), color="blue")
        print()
    except ImportError:
        app(_straceback())

    return "\n".join(err_lines)


def install_config_files(workdir=None, force_reinstall=False):
    """
    Install pre-defined configuration files for the TaskManager and the Scheduler
    in the workdir directory.

    Args:
        workdir: Directory when configuration files should be produced. Use ~/abinit/abipy/ if None
        force_reinstall: Allow overwrite pre-existent configuration files. By default, the function
            raises RuntimeError if configuration files are already present.
    """
    workdir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy") if workdir is None else workdir
    print("Installing configuration files in directory:", workdir)
    from monty.os import makedirs_p
    makedirs_p(workdir)

    scheduler_path = os.path.join(workdir, "scheduler.yaml")
    scheduler_yaml = """
# The launcher will stop submitting jobs when the
# number of jobs in the queue is >= Max number of jobs
max_njobs_inqueue: 2

# Maximum number of cores that can be used by the scheduler.
max_ncores_used: 2

# number of hours to wait.
#hours: 0

# number of minutes to wait.
#minutes: 0

# number of seconds to wait.
seconds: 2

# Send mail to the specified address (accepts string or list of strings).
# PRO TIP: the scheduler WILL try to send and email after a default time of 4 days. If you
#          comment out the mailto address, this will cause the scheduler to terminate, with
#          potentially nefarious effects on your running jobs. If you do not wish to receive
#          emails, a work around is to set the variable `remindme_s` below to something very
#          large (say, 100 days).
#mailto: nobody@nowhere.com

# verbosity level (int, default 0)
#verbose: 0

# The scheduler will shutdown when the number of python exceptions is > max_num_pyexcs
#max_num_pyexcs: 2

# The scheduler will shutdown when the number of Abinit errors is > max_num_abierrs
#max_num_abierrs: 0

# The scheduler will shutdow when the total number of tasks launched is > safety_ratio * tot_num_tasks.
#safety_ratio: 5

# Send an e-mail to mailto every remindme_s seconds.
#remindme_s: 345600
"""

    manager_path = os.path.join(workdir, "manager.yaml")
    manager_yaml = """
qadapters:
    -
      priority: 1
      queue:
        qname: abipy
        qtype: shell
      job:
        mpi_runner: mpirun
        pre_run:
         - export OMP_NUM_THREADS=1
         # IMPORTANT: Change the below line so that the abinit executable is in PATH
         #- export PATH=$HOME/git_repos/abinit/_build/src/98_main:$PATH
         #- ulimit -s unlimited; ulimit -n 2048

      limits:
         min_cores: 1
         max_cores: 2
         timelimit: 0:10:0
      hardware:
         num_nodes: 1
         sockets_per_node: 1
         cores_per_socket: 2
         mem_per_node: 4 Gb
"""

    # Write configuration files.
    if not os.path.isfile(scheduler_path) or force_reinstall:
        with open(scheduler_path, "wt") as fh:
            fh.write(scheduler_yaml)
        print("Scheduler configuration file written to:", scheduler_path)
    else:
        raise RuntimeError("Configuration file: %s already exists.\nUse force_reinstall option to overwrite it" % scheduler_path)

    if not os.path.isfile(manager_path) or force_reinstall:
        with open(manager_path, "wt") as fh:
            fh.write(manager_yaml)
        print("Manager configuration file written to:", manager_path)
    else:
        raise RuntimeError("Configuration file: %s already exists.\nUse force_reinstall option to overwrite it" % manager_path)

    print("""
Configuration files installed successfully.
Please edit the configuration options according to your installation.
In particular, edit the `pre_run` section in manager.yml
so that the abinit executable is in $PATH.
""")

    return 0


def abipy_logo1():
    """http://www.text-image.com/convert/pic2ascii.cgi"""
    return r"""
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
    return r"""
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
    return r"""\
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
