#!/usr/bin/env python
"""Tests for core.density module"""
from __future__ import print_function, division


import tempfile
import abipy.data as abidata

from abipy.core import Density
from abipy.core.testing import *
from abipy.iotools import *


class TestDensity(AbipyTest):
    """Unit tests for Density."""

    def test_ncread_density(self):
        """Read density from NC example data files"""
        assert abidata.DEN_NCFILES

        for path in abidata.DEN_NCFILES:
            print("Reading DEN file %s " % path)

            # Read data directly from file.
            with ETSF_Reader(path) as r:
                nelect_file = r.read_value("number_of_electrons")

            # Compute nelect from data.
            den = Density.from_file(path)
            print(den)
            structure = den.structure
            rhor_tot = den.total_rhor
            rhog_tot = den.total_rhog
            nelect_calc = den.get_nelect().sum()

            # Diff between nelect computed and the one written on file.
            self.assert_almost_equal(nelect_calc, nelect_file)
            self.assert_almost_equal(rhog_tot[0,0,0] * structure.volume, nelect_file)

            if self.which("xcrysden") is not None:
                # Export data in xsf format.
                visu = den.export(".xsf")
                self.assertTrue(callable(visu))

            # Test CUBE files.
            _, tmp_cubefile = tempfile.mkstemp(text=True)
            den.export_to_cube(tmp_cubefile)
            total_den = Density.from_cube(tmp_cubefile)
            assert total_den.structure == den.structure
            assert abs(total_den.get_nelect().sum() - nelect_file) < 1e-3
