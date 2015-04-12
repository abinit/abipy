# coding: utf-8
"""
Common test support for all abipy test scripts.

This single module should provide all the common functionality for abipy tests
in a single location, so that test scripts can just import it and work right away.
"""
from __future__ import print_function, division, unicode_literals

import subprocess
import json

from monty.os.path import which
from monty.json import MontyDecoder
from pymatgen.util.testing import PymatgenTest

import logging
logger = logging.getLogger(__file__)


__all__ = [
    "AbipyTest",
    "AbipyFileTest",
]


def has_abinit(version, cmp=">="):
    """
    Return True if abinit is in $PATH and version is cmp version.
    False if condition is not fulfilled or the execution of `abinit -v`
    raised CalledProcessError
    """
    if which("abinit") is None:
        return False

    try:
        abiver = str(subprocess.check_output(["abinit", "-v"]))

    except subprocess.CalledProcessError:
        # Some MPI implementations require the mpirunner.
        try:
            abiver = subprocess.check_output(["mpirun", "-n", "1", "abinit", "-v"])
        except subprocess.CalledProcessError:
            try:
                abiver = subprocess.check_output(["mpiexec", "-n", "1", "abinit", "-v"])
            except subprocess.CalledProcessError as exc:
                logger.warning(exc.output)
                return False

    return {">=": abiver.strip() >= version.strip(),
            "==": abiver.strip() == version.strip()}[cmp]


def has_matplotlib():
    try:
        import matplotlib.pyplot as plt
        return True
    except ImportError:
        return False


def has_fireworks():
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
    """Extend TestCase with functions from numpy.testing.utils that support ndarrays."""

    @staticmethod
    def which(program):
        """Returns full path to a executable. None if not found or not executable."""
        return which(program)

    @staticmethod
    def has_abinit(version, cmp=">="):
        """Return True if abinit is in $PATH and version is cmp min_version."""
        return has_abinit(version, cmp=cmp)

    def assertFwSerializable(self, obj):
        self.assertTrue('_fw_name' in obj.to_dict())
        self.assertDictEqual(obj.to_dict(), obj.__class__.from_dict(obj.to_dict()).to_dict())


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
