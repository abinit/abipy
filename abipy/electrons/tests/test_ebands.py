"""Tests for electrons.ebands module"""
from __future__ import print_function, division

import abipy.data as data

from abipy.core.kpoints import KpointList
from abipy.electrons.ebands import ElectronsReader, ElectronBands
from abipy.core.testing import *


class EbandsReaderTest(AbipyTest):

    def test_reader(self):
        """Test ElectronsReader."""
        filepath = data.ref_file("si_scf_WFK-etsf.nc")

        with ElectronsReader(filepath) as r:
            kpoints = r.read_kpoints()
            self.assertTrue(isinstance(kpoints, KpointList))
            #self.assertTrue(len(kpoints) == ??)
            #self.assert_all_equal(self.read_nband_sk == ??))

            eigens = r.read_eigenvalues()
            occfacts = r.read_occupations()
            fermie = r.read_fermie()
            self.assertTrue(r.read_nelect() == 8)
            #smearing = r.read_smearing()


class ElectronBandsTest(AbipyTest):

    def test_read_ebands_from_WFK(self):
        """Read ElectronBands from WFK files."""
        for filename in data.WFK_NCFILES:
            ebands = ElectronBands.from_file(filename)
            ebands.to_pymatgen()

            # FIXME
            #print(ebands.as_dict())
            #ElectronBands.from_dict(ebands.as_dict())
            #self.assertMSONable(ebands, test_if_subclass=False)

    def test_read_ebands_from_GSR(self):
        """Read ElectronBands from GSR files."""
        for filename in data.GSR_NCFILES:
            ebands = ElectronBands.from_file(filename)
            ebands.to_pymatgen()

            # FIXME
            #ElectronBands.from_dict(ebands.as_dict())
            #self.assertMSONable(ebands, test_if_subclass=False)

    def test_dos(self):
        """Test DOS methods."""
        gs_bands = ElectronBands.from_file(data.ref_file("si_scf_GSR.nc"))
        dos = gs_bands.get_edos()
        mu = dos.find_mu(8, atol=1.e-4)
        imu = dos.tot_idos.find_mesh_index(mu)
        self.assert_almost_equal(dos.tot_idos[imu][1], 8, decimal=2)

        self.serialize_with_pickle(dos, protocols=[-1], test_eq=False)

    def test_jdos(self):
        """Test JDOS methods."""
        bands = ElectronBands.from_file(data.ref_file("si_scf_GSR.nc"))

        spin = 0
        conduction = [4,]
        for v in range(1,5):
            valence = range(0,v)
            jdos = bands.get_ejdos(spin, valence, conduction)
            intg = jdos.integral()[-1][-1]
            self.assert_almost_equal(intg, len(conduction) * len(valence))

        self.serialize_with_pickle(jdos, protocols=[-1])

        nscf_bands = ElectronBands.from_file(data.ref_file("si_nscf_GSR.nc"))

        # Test the detection of denerate states.
        degs = nscf_bands.degeneracies(spin=0, kpoint=[0,0,0], bands_range=range(8))

        ref_degbands = [[0], [1,2,3], [4,5,6], [7]]

        for i, (e, deg_bands) in enumerate(degs):
            self.assertEqual(deg_bands, ref_degbands[i])

        # JDOS a homogeneous sampling.
        with self.assertRaises(ValueError):
            nscf_bands.get_ejdos(spin, 0, 4)

    def test_pymatgen_converter(self):
        """Test abipy-->pymatgen converter"""
        nscf_bands = ElectronBands.from_file(data.ref_file("si_nscf_GSR.nc"))
        pmg_bands = nscf_bands.to_pymatgen()


if __name__ == "__main__":
    import unittest
    unittest.main()
