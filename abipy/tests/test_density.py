#!/usr/bin/env python
"""Tests for core.density module"""
from __future__ import print_function, division

from abipy.core import Density
from abipy.iotools import *
from abipy.tests import DEN_NCFILES, AbipyTest


class TestDensity(AbipyTest):
    """Unit tests for Density."""

    def test_ncread_density(self):
        """Read density from NC example data files"""
        assert DEN_NCFILES

        for filename in DEN_NCFILES:
            print("Reading file %s " % filename)

            # Read data directly from file.
            with ETSF_Reader(filename) as r:
                nelect_file = r.read_value("number_of_electrons")

            # Compute nelect from data.
            den = Density.from_file(filename)
            print(den)
            den.get_rhor_tot()
            den.get_rhog_tot()
            nelect_calc = den.get_nelect().sum()

            # Diff between nelect computed and the one written on file.
            self.assert_almost_equal(nelect_calc, nelect_file)
            visu = den.export(".xsf")


if __name__ == "__main__":
    import unittest
    unittest.main()

