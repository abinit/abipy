"""
Common test support for all abipy test scripts.

This single module should provide all the common functionality for abipy tests
in a single location, so that test scripts can just import it and work right away.
"""
from __future__ import print_function, division

from unittest import TestCase

import numpy.testing.utils as nptu

__all__ = [
    "AbipyTest",
    "AbipyFileTest",
]

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
