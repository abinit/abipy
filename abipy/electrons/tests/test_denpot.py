# coding: utf-8
"""Tests for DEN/POT files."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.electrons.denpot import DensityNcFile


class DensityNcFileTest(AbipyTest):

    def test_silicon_density_ncfile(self):
        """Testing si_DEN.nc"""
        with DensityNcFile(abidata.ref_file("si_DEN.nc")) as denc:
            repr(denc); str(denc)
            assert denc.structure.formula == "Si2"
            assert denc.ebands.structure == denc.structure
            assert str(denc.xc) == "LDA_XC_TETER93"

            # Kpoint sampling
            # kptopt1 1 ngkpt1 8 8 8 nshiftk1 1 shiftk1   0 0 0
            assert denc.ebands.kpoints.is_mpmesh
            assert not denc.ebands.kpoints.is_path

            ksamp = denc.ebands.kpoints.ksampling
            str(ksamp); repr(ksamp)

            self.assert_equal(ksamp.mpdivs, [8, 8, 8])
            #assert 0

            self.assert_equal(ksamp.kptrlatt_orig, 8 * np.eye(3, 3))
            self.assert_equal(ksamp.shifts_orig, 0)

            self.assert_equal(ksamp.kptrlatt, 8 * np.eye(3, 3))
            self.assert_almost_equal(ksamp.shifts.flatten(), [0.0, 0.0, 0.0])
            assert ksamp.kptopt == 1
            assert ksamp.is_mesh

            density = denc.density
            assert density.nspinor == 1 and density.nsppol == 1 and density.nspden == 1
            self.assert_almost_equal(density.get_nelect(), 8)

            # Test converters.
            denc.write_chgcar(filename=self.get_tmpname(text=True))
            denc.write_xsf(filename=self.get_tmpname(text=True, suffix=".xsf"))
            denc.write_cube(filename=self.get_tmpname(text=True), spin="total")

            # Test ipython notebooks.
            if self.has_nbformat():
                denc.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_nickel_density_ncfile(self):
        """Testing ni_666k_DEN.nc"""
        with abilab.abiopen(abidata.ref_file("ni_666k_DEN.nc")) as denc:
            repr(denc); str(denc)
            assert denc.structure.formula == "Ni1"
            assert denc.ebands.structure == denc.structure
            assert str(denc.xc) == "PBE"

            # Kpoint sampling
            # ngkpt 6 6 6 nshiftk  4
            # shiftk   1/2 1/2 1/2
            #          1/2 0.0 0.0
            #          0.0 1/2 0.0
            #          0.0 0.0 1/2
            assert not denc.ebands.kpoints.is_mpmesh
            assert not denc.ebands.kpoints.is_path
            ksamp = denc.ebands.kpoints.ksampling
            print(ksamp)

            assert ksamp.mpdivs is None
            #assert 0

            self.assert_equal(ksamp.kptrlatt_orig, 6 * np.eye(3, 3))
            self.assert_equal(ksamp.shifts_orig.flatten(), [
                1/2, 1/2, 1/2, 1/2, 0.0, 0.0, 0.0, 1/2, 0.0, 0.0, 0.0, 1/2])

            self.assert_equal(ksamp.kptrlatt.flatten(), [6, -6, 6, -6, 6, 6, -6, -6, 6])
            self.assert_almost_equal(ksamp.shifts.flatten(), [0.5, 0.5, 0.5])

            assert ksamp.kptopt == 1
            assert ksamp.is_mesh

            density = denc.density
            assert density.nspinor == 1 and density.nsppol == 2 and density.nspden == 2
            self.assert_almost_equal(density.get_nelect(), 18)

            # Test converters.
            denc.write_chgcar(filename=self.get_tmpname(text=True))
            denc.write_xsf(filename=self.get_tmpname(text=True, suffix=".xsf"))
            denc.write_cube(filename=self.get_tmpname(text=True), spin="total")

            # Test ipython notebooks.
            if self.has_nbformat():
                denc.write_notebook(nbpath=self.get_tmpname(text=True))
