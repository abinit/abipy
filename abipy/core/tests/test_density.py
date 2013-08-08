#!/usr/bin/env python
"""Tests for core.density module"""
from __future__ import print_function, division

import abipy.data as data 

from abipy.core import Density
from abipy.core.testing import *
from abipy.iotools import *


class TestDensity(AbipyTest):
    """Unit tests for Density."""

    def test_ncread_density(self):
        """Read density from NC example data files"""
        assert data.DEN_NCFILES

        for filename in data.DEN_NCFILES:
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

            # Export data in xsf format.
            visu = den.export(".xsf")
            self.assertTrue(callable(visu))


if __name__ == "__main__":
    import unittest
    unittest.main()
