"""
Common test support for all abipy test scripts.

This single module should provide all the common functionality for abipy tests
in a single location, so that test scripts can just import it and work right away.
"""
from __future__ import print_function, division

from os.path import join as pj, dirname, abspath
from unittest import TestCase

import os
import numpy.testing.utils as nptu

__all__ = [
    "WFK_NCFILES",
    "DEN_NCFILES",
    "ALL_NCFILES",
    "get_ncfile",
    "get_ncfiles_with_ext",
    "get_reference_file",
    "AbipyTest",
    "AbipyFileTest",
]

_DATA_ABSPATH = abspath( pj(dirname(__file__), 'data') )

try:
    _DATA_NCFILES = [pj(_DATA_ABSPATH, f) for f in os.listdir(_DATA_ABSPATH) if f.endswith(".nc")]
except:
    _DATA_NCFILES = list()

def get_ncfile(filename):
    """Return the absolute path of file filename locate in data. None if not found"""
    for f in _DATA_NCFILES:
        if os.path.basename(f) == filename:
            return f
    else:
        return None

def get_ncfiles_with_ext(ext):
    """Return a list with the absolute path of the files with extension ext"""
    ncfiles = []
    for filename in _DATA_NCFILES:
        f = filename.rstrip(".nc").rstrip("-etsf")
        if f.endswith("_"+ext):
            ncfiles.append(filename)
    return ncfiles

def get_reference_file(filename):
    """Returns the absolute path of filename in the tests/data directory."""
    return os.path.join(_DATA_ABSPATH, filename)

def get_datadir():
    return _DATA_ABSPATH

##########################################################################################

WFK_NCFILES = get_ncfiles_with_ext("WFK")

DEN_NCFILES = get_ncfiles_with_ext("DEN")

ALL_NCFILES = WFK_NCFILES + DEN_NCFILES

##########################################################################################

class AbipyTest(TestCase):
    """Extend TestCase with functions from numpy.testing.utils that support ndarrays."""

    @staticmethod
    def assert_almost_equal(actual, desired, decimal=7, err_msg='', verbose=True):
        return nptu.assert_almost_equal(actual, desired, decimal, err_msg, verbose)

    @staticmethod
    def assert_equal(actual, desired, err_msg='', verbose=True):
        return nptu.assert_equal(actual, desired, err_msg=err_msg, verbose=verbose)

class AbipyFileTest(AbipyTest):
    """
    Test class for files with a __str__ attribute.
    At setup, must set the 'file' attribute of the AbipyFileTest.
    """
    file = None

    def assertContains(self, expression):
        """
        Assert that the string representation of the file contains 'expression'
        'expression' is trimmed of leading new line.
        Each line of 'expression' is trimmed of blank spaces.
        Empty lines are ignored.
        """
        def normalize(string):
            string = string.replace('\t', '  ')
            string = string.replace("$", " CASH ")
            string = string.replace("(", " LP ")
            string = string.replace(")", " RP ")
            string = '\n'.join([line.strip() for line in string.splitlines()])
            return string

        expression = normalize(expression)
        ref = normalize(str(self.file))
        return  self.assertRegexpMatches(ref, expression)
