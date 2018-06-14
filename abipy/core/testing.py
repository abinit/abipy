# coding: utf-8
"""
Common test support for all AbiPy test scripts.

This single module should provide all the common functionality for abipy tests
in a single location, so that test scripts can just import it and work right away.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy
import subprocess
import json
import tempfile
import shutil
import unittest
import numpy.testing.utils as nptu
import abipy.data as abidata

from functools import wraps
from monty.os.path import which
from monty.string import is_string
from pymatgen.util.testing import PymatgenTest

import logging
logger = logging.getLogger(__file__)

root = os.path.dirname(__file__)

__all__ = [
    "AbipyTest",
]


def cmp_version(this, other, op=">="):
    """
    Compare two version strings with the given operator ``op``
    >>> assert cmp_version("1.1.1", "1.1.0") and not cmp_version("1.1.1", "1.1.0", op="==")
    """
    from pkg_resources import parse_version
    from monty.operator import operator_from_str
    op = operator_from_str(op)
    return op(parse_version(this), parse_version(other))


#TODO: Replace with abinit build and manager
def has_abinit(version=None, op=">="):
    """
    True if abinit is in $PATH.
    If version is not None, abinit version op version is evaluated and the result is returned.
    False if condition is not fulfilled or the execution of ``abinit -v`` raised CalledProcessError
    """
    abinit = which("abinit")
    if abinit is None: return False
    if version is None: return abinit is not None

    try:
        abinit_version = str(subprocess.check_output(["abinit", "-v"]))

    except subprocess.CalledProcessError:
        # Some MPI implementations require the mpirunner.
        try:
            abinit_version = subprocess.check_output(["mpirun", "-n", "1", "abinit", "-v"])
        except subprocess.CalledProcessError:
            try:
                abinit_version = subprocess.check_output(["mpiexec", "-n", "1", "abinit", "-v"])
            except subprocess.CalledProcessError as exc:
                logger.warning(exc.output)
                return False

    return cmp_version(abinit_version, version, op=op)


def has_matplotlib(version=None, op=">="):
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


def has_seaborn():
    """True if seaborn_ is installed."""
    try:
        import seaborn as sns
        return True
    except ImportError:
        return False


def has_phonopy(version=None, op=">="):
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


def json_read_abinit_input_from_path(json_path):
    """
    Read a json file from the absolute path ``json_path``, return |AbinitInput| instance.
    """
    from abipy.abio.inputs import AbinitInput

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
        print(var)
        print(val_list_r, type(val_list_r[0]))
        print(val_list_t, type(val_list_t[0]))
        for k, var_item in enumerate(val_list_r):
            try:
                error = error or check_var(val_list_t[k], val_list_r[k])
            except IndexError:
                print(val_list_t, type(val_list_t[0]))
                print(val_list_r, type(val_list_r[0]))
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
    # Build GS input file.
    pseudos = abidata.pseudos("14si.pspnc") if usepaw == 0 else data.pseudos("Si.GGA_PBE-JTH-paw.xml")
    #silicon = abilab.Structure.zincblende(5.431, ["Si", "Si"], units="ang")
    silicon = abidata.cif_file("si.cif")

    from abipy.abio.inputs import AbinitInput
    scf_input = AbinitInput(silicon, pseudos)
    ecut = 6
    scf_input.set_vars(
        ecut=ecut,
        pawecutdg=40,
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


class AbipyTest(PymatgenTest):
    """
    Extends PymatgenTest with Abinit-specific methods.
    Several helper functions are implemented as static methods so that we
    can easily reuse the code in the pytest integration tests.
    """

    SkipTest = unittest.SkipTest

    @staticmethod
    def which(program):
        """Returns full path to a executable. None if not found or not executable."""
        return which(program)

    @staticmethod
    def has_abinit(version=None, op=">="):
        """Return True if abinit is in $PATH and version is op min_version."""
        return has_abinit(version=version, op=op)

    @staticmethod
    def has_matplotlib(version=None, op=">="):
        return has_matplotlib(version=version, op=op)

    @staticmethod
    def has_seaborn():
        return has_seaborn()

    @staticmethod
    def has_ase(version=None, op=">="):
        """True if ASE_ package is available."""
        try:
            import ase
        except ImportError:
            return False

        if version is None: return True
        return cmp_version(ase.__version__, version, op=op)

    @staticmethod
    def has_skimage():
        """True if skimage package is available."""
        try:
            from skimage import measure
            return True
        except ImportError:
            return False

    @staticmethod
    def has_python_graphviz(need_dotexec=False):
        """
        True if python-graphviz package is installed and dot executable in path.
        """
        try:
            from graphviz import Digraph
        except ImportError:
            return False

        return which("dot") is not None if need_dotexec else True

    @staticmethod
    def has_mayavi():
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

    @staticmethod
    def get_abistructure_from_abiref(basename):
        """Return an Abipy |Structure| from the basename of one of the reference files."""
        from abipy.core.structure import Structure
        return Structure.as_structure(abidata.ref_file(basename))

    @staticmethod
    def mkdtemp(**kwargs):
        """Invoke mkdtep with kwargs, return the name of a temporary directory."""
        return tempfile.mkdtemp(**kwargs)

    @staticmethod
    def tmpfileindir(basename, **kwargs):
        """
        Return the absolute path of a temporary file with basename ``basename`` created in a temporary directory.
        """
        tmpdir = tempfile.mkdtemp(**kwargs)
        return os.path.join(tmpdir, basename)

    @staticmethod
    def get_tmpname(**kwargs):
        """Invoke mkstep with kwargs, return the name of a temporary file."""
        fd, tmpname = tempfile.mkstemp(**kwargs)
        return tmpname

    @staticmethod
    def has_nbformat():
        """Return True if nbformat is available and we can test the generation of jupyter_ notebooks."""
        try:
            import nbformat
            return True
        except ImportError:
            return False

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
    def json_read_abinit_input(json_basename):
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
    def skip_if_not_phonopy(version=None, op=">="):
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
    def skip_if_not_pseudodojo():
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

    @staticmethod
    def abivalidate_input(abinput, must_fail=False):
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
    def abivalidate_multi(multi):
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
                print(e)

        assert not errors

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
