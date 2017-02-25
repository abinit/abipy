# coding: utf-8
"""Tests for DEN/POT files."""
from __future__ import division, print_function, unicode_literals

import abipy.data as abidata

from abipy.core.testing import *
from abipy.electrons.denpot import DensityNcFile

#from abipy.core.testing import has_matplotlib


class DensityNcFileTest(AbipyTest):

    def test_silicon_density(self):
        """Test si_DEN.nc"""
        with DensityNcFile(abidata.ref_file("si_DEN.nc")) as denc:
            print(denc)
            assert denc.ebands.structure == denc.structure
            assert str(denc.xc) == "LDA_XC_TETER93"
            assert denc.ebands.kpoints.is_mpmesh
            assert not denc.ebands.kpoints.is_path

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
