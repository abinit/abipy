# coding: utf-8
# flake8: noqa
"""
Common test support for all AbiPy test scripts.

This single module should provide all the common functionality for abipy tests
in a single location, so that test scripts can just import it and work right away.
"""
from __future__ import annotations

import os
import numpy
import json
import tempfile
import unittest
#import subprocess
#import time
#import atexit
#import shutil
try:
    import numpy.testing as nptu
except ImportError:
    import numpy.testing.utils as nptu
import abipy.data as abidata

from typing import Optional
from functools import wraps
from monty.os.path import which
from monty.string import is_string
from pymatgen.util.testing import PymatgenTest
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput

root = os.path.dirname(__file__)

__all__ = [
    "AbipyTest",
]


def cmp_version(this: str, other: str, op: str = ">=") -> bool:
    """
    Compare two version strings with the given operator ``op``
    >>> assert cmp_version("1.1.1", "1.1.0") and not cmp_version("1.1.1", "1.1.0", op="==")
    """
    from pkg_resources import parse_version
    from monty.operator import operator_from_str
    op = operator_from_str(op)
    return op(parse_version(this), parse_version(other))


def has_abinit(version: Optional[str] = None, op: str = ">=", manager=None) -> bool:
    """
    True if abinit is available via TaskManager configuration options.
    If version is not None, `abinit_version op version` is evaluated and the result is returned.
    """
    from abipy.flowtk import TaskManager, AbinitBuild
    manager = TaskManager.from_user_config() if manager is None else manager
    build = AbinitBuild(manager=manager)
    if version is None:
        return build.version != "0.0.0"
    else:
        return cmp_version(build.version, version, op=op)


_HAS_MATPLOTLIB_CALLS = 0


def has_matplotlib(version: Optional[str] = None, op: str = ">=") -> bool:
    """
    True if matplotlib_ is installed.
    If version is None, the result of matplotlib.__version__ `op` version is returned.
    """
    try:
        import matplotlib
        # have_display = "DISPLAY" in os.environ
    except ImportError:
        print("Skipping matplotlib test")
        return False

    global _HAS_MATPLOTLIB_CALLS
    _HAS_MATPLOTLIB_CALLS += 1

    if _HAS_MATPLOTLIB_CALLS == 1:
        matplotlib.use("Agg")
        #matplotlib.use("Agg", force=True)  # Use non-graphical display backend during test.

    import matplotlib.pyplot as plt
    # http://stackoverflow.com/questions/21884271/warning-about-too-many-open-figures
    plt.close("all")

    backend = matplotlib.get_backend()
    if backend.lower() != "agg":
        #raise RuntimeError("matplotlib backend now is %s" % backend)
        #matplotlib.use("Agg", warn=True, force=False)
        # Switch the default backend.
        # This feature is experimental, and is only expected to work switching to an image backend.
        plt.switch_backend("Agg")

    if version is None: return True
    return cmp_version(matplotlib.__version__, version, op=op)


def has_plotly(version: Optional[str] = None, op: str = ">=") -> bool:
    """
    True if plotly is installed.
    If version is None, the result of plotly.__version__ `op` version is returned.
    """
    try:
        import plotly
        # have_display = "DISPLAY" in os.environ
    except ImportError:
        print("Skipping plotlyt test")
        return False

    if version is None: return True
    return cmp_version(plotly.__version__, version, op=op)


def has_seaborn() -> bool:
    """True if seaborn_ is installed."""
    try:
        import seaborn as sns
        return True
    except ImportError:
        return False


def has_phonopy(version: Optional[str] = None, op: str = ">=") -> bool:
    """
    True if phonopy_ is installed.
    If version is None, the result of phonopy.__version__ `op` version is returned.
    """
    try:
        import phonopy
    except ImportError:
        print("Skipping phonopy test")
        return False

    if version is None: return True
    return cmp_version(phonopy.__version__, version, op=op)


def get_mock_module():
    """Return mock module for testing. Raises ImportError if not found."""
    try:
        # py > 3.3
        from unittest import mock
    except ImportError:
        try:
            import mock
        except ImportError:
            print("mock module required for unit tests")
            print("Use py > 3.3 or install it with `pip install mock` if py2.7")
            raise

    return mock


def json_read_abinit_input_from_path(json_path: str) -> AbinitInput:
    """
    Read a json file from the absolute path ``json_path``,
    returns: |AbinitInput| instance.
    """
    with open(json_path, "rt") as fh:
        d = json.load(fh)

    # Convert pseudo paths: extract basename and build path in abipy/data/pseudos.
    for pdict in d["pseudos"]:
        pdict["filepath"] = os.path.join(abidata.dirpath, "pseudos", os.path.basename(pdict["filepath"]))

    return AbinitInput.from_dict(d)


def input_equality_check(ref_file, input2, rtol=1e-05, atol=1e-08, equal_nan=False):
    """
    Function to compare two inputs
    ref_file takes the path to reference input in json: json.dump(input.as_dict(), fp, indent=2)
    input2 takes an AbinintInput object
    tol relative tolerance for floats
    we check if all vars are uniquely present in both inputs and if the values are equal (integers, strings)
    or almost equal (floats)
    """
    def check_int(i, j):
        return i != j

    def check_float(x, y):
        return not numpy.isclose(x, y, rtol=rtol, atol=atol, equal_nan=equal_nan)

    def check_str(s, t):
        return s != t

    def check_var(v, w):
        _error = False
        if isinstance(v, int):
            _error = check_int(v, w)
        elif isinstance(v, float):
            _error = check_float(v, w)
        elif is_string(v):
            _error = check_str(v, w)
        return _error

    def flatten_var(o, tree_types=(list, tuple, numpy.ndarray)):
        flat_var = []
        if isinstance(o, tree_types):
            for value in o:
                for sub_value in flatten_var(value, tree_types):
                    flat_var.append(sub_value)
        else:
            flat_var.append(o)
        return flat_var

    input_ref = json_read_abinit_input_from_path(os.path.join(root, '..', 'test_files', ref_file))

    errors = []
    diff_in_ref = [var for var in input_ref.vars if var not in input2.vars]
    diff_in_actual = [var for var in input2.vars if var not in input_ref.vars]
    if len(diff_in_ref) > 0 or len(diff_in_actual) > 0:
        error_description = 'not the same input parameters:\n' \
                            '     %s were found in ref but not in actual\n' \
                            '     %s were found in actual but not in ref\n' % \
                            (diff_in_ref, diff_in_actual)
        errors.append(error_description)

    for var, val_r in input_ref.vars.items():
        try:
            val_t = input2.vars[var]
        except KeyError:
            errors.append('variable %s from the reference is not in the actual input\n' % str(var))
            continue
        val_list_t = flatten_var(val_t)
        val_list_r = flatten_var(val_r)
        error = False
        #print(var)
        #print(val_list_r, type(val_list_r[0]))
        #print(val_list_t, type(val_list_t[0]))
        for k, var_item in enumerate(val_list_r):
            try:
                error = error or check_var(val_list_t[k], val_list_r[k])
            except IndexError:
                #print(val_list_t, type(val_list_t[0]))
                #print(val_list_r, type(val_list_r[0]))
                raise RuntimeError('two value lists were not flattened in the same way, try to add the collection'
                                   'type to the tree_types tuple in flatten_var')

        if error:
            error_description = 'var %s differs: %s (reference) != %s (actual)' % \
                                (var, val_r, val_t)
            errors.append(error_description)

    if input2.structure != input_ref.structure:
        errors.append('Structures are not the same.\n')
        print(input2.structure, input_ref.structure)

    if len(errors) > 0:
        msg = 'Two inputs were found to be not equal:\n'
        for err in errors:
            msg += '   ' + err + '\n'
        raise AssertionError(msg)


def get_gsinput_si(usepaw=0, as_task=False):
    """
    Build and return a GS input file for silicon or a Task if `as_task`
    """
    pseudos = abidata.pseudos("14si.pspnc") if usepaw == 0 else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    silicon = abidata.cif_file("si.cif")

    from abipy.abio.inputs import AbinitInput
    scf_input = AbinitInput(silicon, pseudos)
    ecut = 6
    scf_input.set_vars(
        ecut=ecut,
        nband=6,
        paral_kgb=0,
        iomode=3,
        toldfe=1e-9,
    )
    if usepaw:
        scf_input.set_vars(pawecutdg=4 * ecut)

    # K-point sampling (shifted)
    scf_input.set_autokmesh(nksmall=4)

    if not as_task:
        return scf_input
    else:
        from abipy.flowtk.tasks import ScfTask
        return ScfTask(scf_input)


def get_gsinput_alas_ngkpt(ngkpt, usepaw=0, as_task=False):
    """
    Build and return a GS input file for AlAs or a Task if `as_task`
    """
    if usepaw != 0: raise NotImplementedError("PAW")
    pseudos = abidata.pseudos("13al.981214.fhi", "33as.pspnc")
    structure = abidata.structure_from_ucell("AlAs")

    scf_input = AbinitInput(structure, pseudos=pseudos)

    scf_input.set_vars(
        nband=5,
        ecut=8.0,
        ngkpt=ngkpt,
        nshiftk=1,
        shiftk=[0, 0, 0],
        tolvrs=1.0e-6,
        diemac=12.0,
    )

    if not as_task:
        return scf_input
    else:
        from abipy.flowtk.tasks import ScfTask
        return ScfTask(scf_input)


class AbipyTest(PymatgenTest):
    """
    Extends PymatgenTest with Abinit-specific methods.
    Several helper functions are implemented as static methods so that we
    can easily reuse the code in the pytest integration tests.
    """

    SkipTest = unittest.SkipTest

    @staticmethod
    def which(program: str) -> bool:
        """Returns full path to a executable. None if not found or not executable."""
        return which(program)

    @staticmethod
    def has_abinit(version: Optional[str] = None, op: str = ">=") -> bool:
        """Return True if abinit is in $PATH and version is op min_version."""
        return has_abinit(version=version, op=op)

    def skip_if_abinit_not_ge(self, version: str) -> None:
        """Skip test if Abinit version is not >= `version`"""
        op = ">="
        if not self.has_abinit(version, op=op):
            raise unittest.SkipTest("This test requires Abinit version %s %s" % (op, version))

    @staticmethod
    def has_matplotlib(version: Optional[str] = None, op: str = ">=") -> bool:
        return has_matplotlib(version=version, op=op)

    @staticmethod
    def has_plotly(version: Optional[str] = None, op: str = ">=") -> bool:
        return has_plotly(version=version, op=op)

    @staticmethod
    def has_seaborn() -> bool:
        return has_seaborn()

    @staticmethod
    def has_ase(version: Optional[str] = None, op: str = ">=") -> bool:
        """True if ASE_ package is available."""
        try:
            import ase
        except ImportError:
            return False

        if version is None: return True
        return cmp_version(ase.__version__, version, op=op)

    @staticmethod
    def has_ifermi() -> bool:
        """True if ifermi package is available."""
        try:
            from ifermi.interpolate import FourierInterpolator
            return True
        except ImportError:
            return False

    @staticmethod
    def has_skimage() -> bool:
        """True if skimage package is available."""
        try:
            from skimage import measure
            return True
        except ImportError:
            return False

    @staticmethod
    def has_python_graphviz(need_dotexec: bool = True) -> bool:
        """
        True if python-graphviz package is installed and dot executable in path.
        """
        try:
            from graphviz import Digraph
        except ImportError:
            return False

        return which("dot") is not None if need_dotexec else True

    @staticmethod
    def has_mayavi() -> bool:
        """
        True if mayavi_ is available. Set also offscreen to True
        """
        # Disable mayavi for the time being.
        #return False
        # This to run mayavi tests only on Travis
        if not os.environ.get("TRAVIS"): return False
        try:
            from mayavi import mlab
        except ImportError:
            return False

        #mlab.clf()
        mlab.options.offscreen = True
        mlab.options.backend = "test"
        return True

    def has_panel(self) -> bool:
        """False if Panel library is not installed."""
        try:
            import param
            import panel as pn
            import bokeh
            return pn
        except ImportError:
            return False

    def has_networkx(self) -> bool:
        """False if networkx library is not installed."""
        try:
            import networkx as nx
            return nx
        except ImportError:
            return False

    def has_graphviz(self) -> bool:
        """True if graphviz library is installed and `dot` in $PATH"""
        try:
            from graphviz import Digraph
            import graphviz
        except ImportError:
            return False

        if self.which("dot") is None: return False
        return graphviz

    def has_phonopy(self, version: Optional[str] = None, op: str = ">=") -> bool:
        """
        True if phonopy_ is installed.
        If version is None, the result of phonopy.__version__ `op` version is returned.
        """
        return has_phonopy(version=version, op=op)

    @staticmethod
    def get_abistructure_from_abiref(basename: str) -> Structure:
        """Return an Abipy |Structure| from the basename of one of the reference files."""
        from abipy.core.structure import Structure
        return Structure.as_structure(abidata.ref_file(basename))

    @staticmethod
    def mkdtemp(**kwargs) -> str:
        """Invoke mkdtep with kwargs, return the name of a temporary directory."""
        return tempfile.mkdtemp(**kwargs)

    @staticmethod
    def tmpfileindir(basename: str, **kwargs) -> str:
        """
        Return the absolute path of a temporary file with basename ``basename`` created in a temporary directory.
        """
        tmpdir = tempfile.mkdtemp(**kwargs)
        return os.path.join(tmpdir, basename)

    @staticmethod
    def get_tmpname(**kwargs) -> str:
        """Invoke mkstep with kwargs, return the name of a temporary file."""
        _, tmpname = tempfile.mkstemp(**kwargs)
        return tmpname

    def tmpfile_write(self, string: str) -> str:
        """
        Write string to a temporary file. Returns the name of the temporary file.
        """
        fd, tmpfile = tempfile.mkstemp(text=True)

        with open(tmpfile, "wt") as fh:
            fh.write(string)

        return tmpfile

    @staticmethod
    def has_nbformat() -> bool:
        """Return True if nbformat is available and we can test the generation of jupyter_ notebooks."""
        try:
            import nbformat
            return True
        except ImportError:
            return False

    #def run_nbpath(self, nbpath: str):
    #    """Test that the notebook in question runs all cells correctly."""
    #    nb, errors = notebook_run(nbpath)
    #    return nb, errors

    @staticmethod
    def has_ipywidgets():
        """Return True if ipywidgets_ package is available."""
        # Disabled due to:
        # AttributeError: 'NoneType' object has no attribute 'session'
        return False
        # Disable widget tests on TRAVIS
        #if os.environ.get("TRAVIS"): return False
        try:
            import ipywidgets as ipw
            return True
        except ImportError:
            return False

    @staticmethod
    def assert_almost_equal(actual, desired, decimal=7, err_msg='', verbose=True):
        """
        Alternative naming for assertArrayAlmostEqual.
        """
        return nptu.assert_almost_equal(actual, desired, decimal, err_msg, verbose)

    @staticmethod
    def assert_equal(actual, desired, err_msg='', verbose=True):
        """
        Alternative naming for assertArrayEqual.
        """
        return nptu.assert_equal(actual, desired, err_msg=err_msg, verbose=verbose)

    @staticmethod
    def json_read_abinit_input(json_basename: str):
        """Return an |AbinitInput| from the basename of the file in abipy/data/test_files."""
        return json_read_abinit_input_from_path(os.path.join(root, '..', 'test_files', json_basename))

    @staticmethod
    def assert_input_equality(ref_basename, input_to_test, rtol=1e-05, atol=1e-08, equal_nan=False):
        """
        Check equality between an input and a reference in test_files.
        only input variables and structure are compared.

        Args:
            ref_basename: base name of the reference file to test against in test_files
            input_to_test: |AbinitInput| object to test
            rtol: passed to numpy.isclose for float comparison
            atol: passed to numpy.isclose for float comparison
            equal_nan: passed to numpy.isclose for float comparison

        Returns:
            raises an assertion error if the two inputs are not the same
        """
        ref_file = os.path.join(root, '..', 'test_files', ref_basename)
        input_equality_check(ref_file, input_to_test, rtol=rtol, atol=atol, equal_nan=equal_nan)

    @staticmethod
    def straceback():
        """Returns a string with the traceback."""
        import traceback
        return traceback.format_exc()

    @staticmethod
    def skip_if_not_phonopy(version: Optional[str] = None, op: str = ">=") -> None:
        """
        Raise SkipTest if phonopy_ is not installed.
        Use ``version`` and ``op`` to ask for a specific version
        """
        if not has_phonopy(version=version, op=op):
            if version is None:
                msg = "This test requires phonopy"
            else:
                msg = "This test requires phonopy version %s %s" % (op, version)
            raise unittest.SkipTest(msg)

    @staticmethod
    def skip_if_not_bolztrap2(version: Optional[str] = None, op: str = ">=") -> None:
        """
        Raise SkipTest if bolztrap2 is not installed.
        Use ``version`` and ``op`` to ask for a specific version
        """
        try:
            import BoltzTraP2 as bzt
        except ImportError:
            raise unittest.SkipTest("This test requires bolztrap2")

        from BoltzTraP2.version import PROGRAM_VERSION
        if version is not None and not cmp_version(PROGRAM_VERSION, version, op=op):
            msg = "This test requires bolztrap2 version %s %s" % (op, version)
            raise unittest.SkipTest(msg)

    def skip_if_not_executable(self, executable: str) -> None:
        """
        Raise SkipTest if executable is not installed.
        """
        if self.which(executable) is None:
            raise unittest.SkipTest("This test requires `%s` in PATH" % str(executable))

    @staticmethod
    def skip_if_not_pseudodojo() -> None:
        """
        Raise SkipTest if pseudodojo package is not installed.
        """
        try:
            from pseudo_dojo import OfficialTables
        except ImportError:
            raise unittest.SkipTest("This test requires pseudodojo package.")

    @staticmethod
    def get_mock_module():
        """Return mock module for testing. Raises ImportError if not found."""
        return get_mock_module()

    def decode_with_MSON(self, obj):
        """
        Convert obj into JSON assuming MSONable protocol. Return new object decoded with MontyDecoder
        """
        from monty.json import MSONable, MontyDecoder
        self.assertIsInstance(obj, MSONable)
        return json.loads(obj.to_json(), cls=MontyDecoder)

    @staticmethod
    def abivalidate_input(abinput: AbinitInput, must_fail: bool = False) -> None:
        """
        Invoke Abinit to test validity of an |AbinitInput| object
        Print info to stdout if failure before raising AssertionError.
        """
        v = abinput.abivalidate()
        if must_fail:
            assert v.retcode != 0 and v.log_file.read()
        else:
            if v.retcode != 0:
                print("type abinput:", type(abinput))
                print("abinput:\n", abinput)
                lines = v.log_file.readlines()
                i = len(lines) - 50 if len(lines) >= 50 else 0
                print("Last 50 line from logfile:")
                print("".join(lines[i:]))

            assert v.retcode == 0

    @staticmethod
    def abivalidate_multi(multi) -> None:
        """
        Invoke Abinit to test validity of a |MultiDataset| or a list of |AbinitInput| objects.
        """
        if hasattr(multi, "split_datasets"):
            inputs = multi.split_datasets()
        else:
            inputs = multi

        errors = []
        for inp in inputs:
            try:
                AbipyTest.abivalidate_input(inp)
            except Exception as exc:
                errors.append(AbipyTest.straceback())
                errors.append(str(exc))

        if errors:
            for e in errors:
                print(90 * "=")
                print(e)
                print(90 * "=")

        assert not errors

    def abivalidate_work(self, work):
        """Invoke Abinit to test validity of the inputs of a |Work|"""
        from abipy.flowtk import Flow
        tmpdir = tempfile.mkdtemp()
        flow = Flow(workdir=tmpdir)
        flow.register_work(work)
        return self.abivalidate_flow(flow)

    @staticmethod
    def abivalidate_flow(flow):
        """
        Invoke Abinit to test validity of the inputs of a |Flow|
        """
        isok, errors = flow.abivalidate_inputs()
        if not isok:
            for e in errors:
                if e.retcode == 0: continue
                #print("type abinput:", type(abinput))
                #print("abinput:\n", abinput)
                lines = e.log_file.readlines()
                i = len(lines) - 50 if len(lines) >= 50 else 0
                print("Last 50 line from logfile:")
                print("".join(lines[i:]))
            raise RuntimeError("flow.abivalidate_input failed. See messages above.")

    @staticmethod
    @wraps(get_gsinput_si)
    def get_gsinput_si(*args, **kwargs):
        return get_gsinput_si(*args, **kwargs)

    @staticmethod
    @wraps(get_gsinput_alas_ngkpt)
    def get_gsinput_alas_ngkpt(*args, **kwargs):
        return get_gsinput_alas_ngkpt(*args, **kwargs)


ABIPY_TESTDB_NAME = "abipy_unit_tests"


def abipy_has_mongodb(host='localhost', port=27017, name=ABIPY_TESTDB_NAME, username=None, password=None):
    try:
        from pymongo import MongoClient
        connection = MongoClient(host, port, j=True)
        db = connection[name]
        if username:
            db.authenticate(username, password)
        return True
    except Exception:
        return False


class AbipyTestWithMongoDb(AbipyTest):
    """A suite of tests requiring a MongoDB database."""

    #def has_mongodb(self):
    #    """True if mongodb server is reachable."""
    #    return abipy_has_mongodb()

    #@classmethod
    #def setUpClass(cls):

    #@classmethod
    #def tearDownClass(cls):

    #@classmethod
    #def setup_mongodb(cls):
    #    try:
    #        cls._connection = connect(db=TESTDB_NAME)
    #        cls._connection.drop_database(TESTDB_NAME)
    #        cls.db = get_db()
    #    except Exception:
    #        cls.db = None
    #        cls._connection = None
    #
    #@classmethod
    #def teardown_mongodb(cls):
    #    if cls._connection:
    #        cls._connection.drop_database(TESTDB_NAME)

#class MongoTemporaryInstance:
#    """Singleton to manage a temporary MongoDB instance
#
#    Use this for testing purpose only. The instance is automatically destroyed
#    at the end of the program.
#
#    """
#    _instance = None
#
#    @classmethod
#    def get_instance(cls):
#        if cls._instance is None:
#            cls._instance = cls()
#            atexit.register(cls._instance.shutdown)
#        return cls._instance
#
#    def __init__(self):
#        self._tmpdir = tempfile.mkdtemp()
#        self._process = subprocess.Popen(['mongod', '--bind_ip', 'localhost',
#                                          '--port', str(MONGODB_TEST_PORT),
#                                          '--dbpath', self._tmpdir,
#                                          '--nojournal', '--nohttpinterface',
#                                          '--noauth', '--smallfiles',
#                                          '--syncdelay', '0',
#                                          '--maxConns', '10',
#                                          '--nssize', '1', ],
#                                         stdout=open(os.devnull, 'wb'),
#                                         stderr=subprocess.STDOUT)
#
#        # wait for the instance to be ready
#        # Mongo is ready in a glance, we just wait to be able to open a
#        # Connection.
#        import pymongo
#        for i in range(3):
#            time.sleep(0.1)
#            try:
#                self._conn = pymongo.Connection('localhost', MONGODB_TEST_PORT)
#            except pymongo.errors.ConnectionFailure:
#                continue
#            else:
#                break
#        else:
#            self.shutdown()
#            assert False, 'Cannot connect to the mongodb test instance'
#
#    @property
#    def conn(self):
#        return self._conn
#
#    def shutdown(self):
#        if self._process:
#            self._process.terminate()
#            self._process.wait()
#            self._process = None
#            shutil.rmtree(self._tmpdir, ignore_errors=True)


#def notebook_run(path):
#    """
#    Execute a notebook via nbconvert and collect output.
#
#    Taken from: https://blog.thedataincubator.com/2016/06/testing-jupyter-notebooks/
#
#    Args:
#        path (str): file path for the notebook object
#
#    Returns: (parsed nb object, execution errors)
#    """
#    import nbformat
#    dirname, _ = os.path.split(path)
#    os.chdir(dirname)
#    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
#        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
#                "--ExecutePreprocessor.timeout=300",
#                "--ExecutePreprocessor.allow_errors=True",
#                "--output", fout.name, path]
#        subprocess.check_call(args)
#
#        fout.seek(0)
#        nb = nbformat.read(fout, nbformat.current_nbformat)
#
#    errors = [output for cell in nb.cells if "outputs" in cell
#              for output in cell["outputs"] if output.output_type == "error"]
#
#    return nb, errors
