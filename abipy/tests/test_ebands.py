"""Tests for electrons.ebands module"""
from __future__ import print_function, division

from abipy.electrons import ElectronBands
from abipy.tests import *


class ElectronBandsTest(AbipyTest):

    def setUp(self):
        assert WFK_NCFILES

    def test_ncread_ebands(self):
        """Read ElectronBands from WFK data files"""
        for filename in WFK_NCFILES:
            bands = ElectronBands.from_file(filename)

    def test_dos(self):
        """Test DOS methods."""
        gs_bands = ElectronBands.from_file(get_reference_file("si_WFK-etsf.nc"))
        dos = gs_bands.get_dos()
        mu = dos.find_mu(8, atol=1.e-4)
        self.assert_almost_equal(mu, 6.1443350264585996, decimal=4)

    def test_jdos(self):
        """Test JDOS methods."""
        bands = ElectronBands.from_file(get_reference_file("si_WFK-etsf.nc"))

        spin = 0
        conduction = [4,]
        for v in range(1,5):
            valence = range(0,v)
            jdos = bands.get_jdos(spin, valence, conduction)
            intg = jdos.integral()[-1][-1]
            self.assert_almost_equal(intg, len(conduction) * len(valence))

        # Need a homogeneous sampling.
        with self.assertRaises(ValueError):
            bands = ElectronBands.from_file(get_reference_file("si_nscf_WFK-etsf.nc"))
            bands.get_jdos(spin, 0, 4)

