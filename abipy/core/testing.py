# coding: utf-8
"""
Common test support for all abipy test scripts.

This single module should provide all the common functionality for abipy tests
in a single location, so that test scripts can just import it and work right away.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import subprocess
import json
import tempfile
import numpy.testing.utils as nptu

from monty.os.path import which
from pymatgen.util.testing import PymatgenTest

import logging
logger = logging.getLogger(__file__)


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
        #have_display = "DISPLAY" in os.environ
        import matplotlib
        matplotlib.use("Agg")  # Use non-graphical display backend during test.
        if version is None: return True

    except ImportError:
        print("Skipping matplotlib test")
        return False

    return cmp_version(matplotlib.__version__, version, op=op)


def has_seaborn():
    """True if seaborn is installed."""
    try:
        import seaborn.apionly as sns
        return True
    except ImportError:
        return False


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


class AbipyTest(PymatgenTest):
    """Extends PymatgenTest with Abinit-specific methods """

    @staticmethod
    def which(program):
        """Returns full path to a executable. None if not found or not executable."""
        return which(program)

    @staticmethod
    def has_abinit(version=None, op=">="):
        """Return True if abinit is in $PATH and version is op min_version."""
        return has_abinit(version=None, op=op)

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

    def assertFwSerializable(self, obj):
        self.assertTrue('_fw_name' in obj.to_dict())
        self.assertDictEqual(obj.to_dict(), obj.__class__.from_dict(obj.to_dict()).to_dict())

    def get_abistructure_from_abiref(self, basename):
        """Return an Abipy structure from the basename of one of the reference files."""
        import abipy.data as abidata
        from abipy.core.structure import Structure
        return Structure.as_structure(abidata.ref_file(basename))

    def get_tmpname(self, **kwargs):
        """Invoke mkstep with kwargs, return the name of a temporary file."""
        fd, tmpname = tempfile.mkstemp(**kwargs)
        return tmpname

    def has_nbformat(self):
        """Return True if nbformat is available and we can test the generation of ipython notebooks."""
        try:
            import nbformat
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
