"""Tests for electrons.ebands module"""
from __future__ import print_function, division

from abipy.core.kpoints import KpointList
from abipy.electrons.ebands import Ebands_Reader, ElectronBands
from abipy.tests import *

class EbandsReaderTest(AbipyTest):

    def test_reader(self):
        """Test Ebands_Reader."""
        filepath = get_reference_file("si_WFK-etsf.nc")

        with Ebands_Reader(filepath) as r:
            kpoints = r.read_kpoints()
            self.assertTrue(isinstance(kpoints, KpointList))
            #self.assertTrue(len(kpoints) == ??)
            #self.assert_all_equal(self.read_nband_sk == ??))

            eigens = r.read_eigens()
            occfacts = r.read_occfacts()
            fermie = r.read_fermie()

            #self.assertTrue(r.read_nelect() == 8)
            #smearing = rread_smearing()


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


if __name__ == "__main__":
    import unittest
    unittest.main()
