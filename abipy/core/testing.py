# coding: utf-8
"""
Common test support for all abipy test scripts.

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

from monty.os.path import which
from monty.string import is_string
from pymatgen.util.testing import PymatgenTest

import logging
logger = logging.getLogger(__file__)

root = os.path.dirname(__file__)

__all__ = [
    "AbipyTest",
    "AbipyFileTest",
]


def cmp_version(this, other, op=">="):
    """
    Compare two version strings with the given operator `op`

    >>> assert cmp_version("1.1.1", "1.1.0") and not cmp_version("1.1.1", "1.1.0", op="==")
    """
    from pkg_resources import parse_version
    from monty.operator import operator_from_str
    op = operator_from_str(op)
    return op(parse_version(this), parse_version(other))


def has_abinit(version=None, op=">="):
    """
    True if abinit is in $PATH.
    If version is not None, abinit version op version is evaluated and the result is returned.
    False if condition is not fulfilled or the execution of `abinit -v` raised CalledProcessError
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
    True if matplotlib is installed.
    If version is None, the result of matplotlib.__version__ `op` version is returned.
    """
    try:
        # have_display = "DISPLAY" in os.environ
        import matplotlib
        matplotlib.use("Agg")  # Use non-graphical display backend during test.

    except ImportError:
        print("Skipping matplotlib test")
        return False

    # http://stackoverflow.com/questions/21884271/warning-about-too-many-open-figures
    import matplotlib.pyplot as plt
    plt.close("all")

    if version is None: return True
    return cmp_version(matplotlib.__version__, version, op=op)


def has_seaborn():
    """True if seaborn is installed."""
    try:
        import seaborn.apionly as sns
        return True
    except ImportError:
        return False


def has_phonopy(version=None, op=">="):
    """
    True if phonopy is installed.
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


def has_fireworks():
    """True if fireworks is installed."""
    try:
        import fireworks
        return True
    except ImportError:
        return False


def has_mongodb(host='localhost', port=27017, name='mongodb_test', username=None, password=None):
    try:
        from pymongo import MongoClient
        connection = MongoClient(host, port, j=True)
        db = connection[name]
        if username:
            db.authenticate(username, password)

        return True
    except:
        return False


def json_read_abinit_input_from_path(json_path):
    """
    Read a json file from the absolute path `json_path`, return AbinitInput instance.
    """
    from abipy.abio.inputs import AbinitInput
    import abipy.data as abidata

    with open(json_path, "rt") as fh:
        d = json.load(fh)

    # Convert pseudo paths: extract basename and build path in abipy/data/pseudos.
    for pdict in d["pseudos"]:
        pdict["filepath"] = os.path.join(abidata.dirpath, "pseudos", os.path.basename(pdict["filepath"]))

    return AbinitInput.from_dict(d)


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


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


class AbipyTest(PymatgenTest):
    """Extends PymatgenTest with Abinit-specific methods """

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
        """True if ASE package is available."""
        try:
            import ase
        except ImportError:
            return False

        if version is None: return True
        return cmp_version(ase.__version__, version, op=op)

    def has_mayavi(self):
        """
        True if Mayavi is available. Set also offscreen to True
        """
        try:
            from mayavi import mlab
        except ImportError:
            return False

        #mlab.clf()
        mlab.options.offscreen = True
        return True

    def assertFwSerializable(self, obj):
        self.assertTrue('_fw_name' in obj.to_dict())
        self.assertDictEqual(obj.to_dict(), obj.__class__.from_dict(obj.to_dict()).to_dict())

    @staticmethod
    def get_abistructure_from_abiref(basename):
        """Return an Abipy structure from the basename of one of the reference files."""
        import abipy.data as abidata
        from abipy.core.structure import Structure
        return Structure.as_structure(abidata.ref_file(basename))

    @staticmethod
    def mkdtemp(**kwargs):
        """Invoke mkdtep with kwargs, return the name of a temporary directory."""
        return tempfile.mkdtemp(**kwargs)

    @staticmethod
    def tmpfileindir(basename, **kwargs):
        """
        Return the absolute path of a temporary file with basename `basename` created in a temporary directory.
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
        """Return True if nbformat is available and we can test the generation of ipython notebooks."""
        try:
            import nbformat
            return True
        except ImportError:
            return False

    @staticmethod
    def has_ipywidgets():
        """Return True if ipywidgets is available."""
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
        """Return an AbinitInput from the basename of the file in abipy/data/test_files."""
        return json_read_abinit_input_from_path(os.path.join(root, '..', 'test_files', json_basename))

    @staticmethod
    def assert_input_equality(ref_basename, input_to_test, rtol=1e-05, atol=1e-08, equal_nan=False):
        """
        Check equality between an input and a reference in test_files.
        only input variables and structure are compared.
        Args:
            ref_basename: base name of the reference file to test against in test_files
            input_to_test: AbinitInput object to test
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
        return straceback()

    @staticmethod
    def skip_if_not_phonopy(version=None, op=">="):
        """
        Raise SkipTest if phonopy is not installed.
        Use `version` and `op` to ask for a specific version
        """
        if not has_phonopy(version=version, op=op):
            if version is None:
                msg = "This test requires phonopy"
            else:
                msg = "This test requires phonopy version %s %s" % (op, version)
            raise unittest.SkipTest(msg)

    @staticmethod
    def skip_if_not_abinit(version=None, op=">="):
        """
        Raise SkipTest if the specified version of abinit is not installed.
        Use `version` and `op` to ask for a specific version
        """
        if not has_abinit(version=version, op=op):
            if version is None:
                msg = "This test requires abinit"
            else:
                msg = "This test requires abinit version %s %s" % (op, version)
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

    def abivalidate_input(self, abinput, must_fail=False):
        """
        Invoke Abinit to test validity of an Input object
        Print info to stdout if failuer before raising AssertionError.
        """
        v = abinput.abivalidate()
        if must_fail:
            assert v.retcode != 0 and v.log_file.read()
        else:
            if v.retcode != 0:
                print(type(abinput))
                print(abinput)
                lines = v.log_file.readlines()
                i = len(lines) - 50 if len(lines) >= 50 else 0
                print("Last 50 line from logfile:")
                print("".join(lines[i:]))

            assert v.retcode == 0

    def abivalidate_multi(self, multi):
        """
        Invoke Abinit to test validity of a `MultiDataset` or a list of input objects.
        """
        if hasattr(multi, "split_datasets"):
            inputs = multi.split_datasets()
        else:
            inputs = multi

        errors = []
        for inp in inputs:
            try:
                self.abivalidate_input(inp)
            except Exception as exc:
                errors.append(str(exc))

        if errors:
            for e in errors:
                print(e)

        assert not errors


class AbipyFileTest(AbipyTest):
    """
    Test class for files with a __str__ attribute.
    At setup, must set the 'file' attribute of the AbipyFileTest.
    """
    file = None

    @staticmethod
    def normalize(string):
        string = string.replace('\t', '  ')
        string = string.replace("$", " CASH ")
        string = string.replace("(", " LP ")
        string = string.replace(")", " RP ")
        string = string.replace("*", " STAR ")
        string = string.strip()
        string = '\n'.join([line.strip() for line in string.splitlines()])

        return string

    def assertContains(self, expression):
        """
        Assert that the string representation of the file contains 'expression'
        'expression' is trimmed of leading new line.
        Each line of 'expression' is trimmed of blank spaces.
        Empty lines are ignored.
        """
        expression = self.normalize(expression)
        ref = self.normalize(str(self.file))

        return self.assertRegexpMatches(ref, expression)

    def assertContainsNot(self, expression):
        """
        Assert that the string representation of the file does not contain
        'expression'.
        """
        expression = self.normalize(expression)
        ref = self.normalize(str(self.file))

        return self.assertNotRegexpMatches(ref, expression)

    def assertEmpty(self):
        """Assert the string representation is empty."""
        s = str(self.file).strip()
        self.assertFalse(bool(s))

    def assertOrder(self, expression1, expression2):
        """
        Assert that the string representation of the file
        contains 'expression1' before 'expression2'.
        """
        expression1 = self.normalize(expression1)
        expression2 = self.normalize(expression2)
        ref = self.normalize(str(self.file))
        self.assertRegexpMatches(ref, expression1)
        self.assertRegexpMatches(ref, expression2)
        last = ref.split(expression1)[-1]

        return self.assertRegexpMatches(last, expression2)


CONF_FILE = None
BKP_FILE = None


def change_matplotlib_backend(new_backend=""):
    """Change the backend by modifying the matplotlib configuration file."""
    global CONF_FILE, BKP_FILE

    if not new_backend:
        return

    home = os.environ["HOME"]
    CONF_FILE = conf_file = os.path.join(home, ".matplotlib", "matplotlibrc")

    BKP_FILE = conf_file + ".bkp"

    if os.path.exists(conf_file):
        shutil.copy(conf_file, BKP_FILE)

        with open(conf_file, "rt") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line.strip().startswith("backend"):
                lines[i] = "backend : " + new_backend + "\n"

        with open(conf_file, "w") as f:
            f.writelines(lines)


def revert_matplotlib_backend():
    """Revert matplotlib backend to the previous value."""
    global CONF_FILE, BKP_FILE
    # print("reverting: BKP_FILE %s --> CONF %s" % (BKP_FILE, CONF_FILE))
    if BKP_FILE is not None:
        shutil.move(BKP_FILE, CONF_FILE)
