"""Tests for electrons.ebands module"""
from __future__ import print_function, division

import numpy as np
import abipy.data as data

import pymatgen.core.units as units
from abipy.core.kpoints import KpointList
from abipy.electrons.ebands import ElectronsReader, ElectronBands, ElectronDos
from abipy.core.testing import *


class EbandsReaderTest(AbipyTest):

    def test_reader(self):
        """Test ElectronsReader."""
        filepath = data.ref_file("si_scf_WFK.nc")

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
            ebands.to_pdframe()
            assert ElectronBands.as_ebands(ebands) is ebands

            # FIXME
            #print(ebands.as_dict())
            #ElectronBands.from_dict(ebands.as_dict())
            #self.assertMSONable(ebands, test_if_subclass=False)

    def test_read_ebands_from_GSR(self):
        """Read ElectronBands from GSR files."""
        for filename in data.GSR_NCFILES:
            ebands = ElectronBands.from_file(filename)
            ebands.to_pymatgen()
            ebands.to_pdframe()

            # FIXME
            #ElectronBands.from_dict(ebands.as_dict())
            #self.assertMSONable(ebands, test_if_subclass=False)

    def test_dos(self):
        """Test DOS methods."""
        gs_bands = ElectronBands.from_file(data.ref_file("si_scf_GSR.nc"))
        dos = gs_bands.get_edos()
        print(dos)
        assert ElectronDos.as_edos(dos, {}) is dos
        edos_samevals = ElectronDos.as_edos(gs_bands, {})
        assert ElectronDos.as_edos(gs_bands, {}) == dos
        assert ElectronDos.as_edos(data.ref_file("si_scf_GSR.nc"), {}) == dos

        mu = dos.find_mu(8)
        imu = dos.tot_idos.find_mesh_index(mu)
        self.assert_almost_equal(dos.tot_idos[imu][1], 8, decimal=2)

        self.serialize_with_pickle(dos, protocols=[-1], test_eq=False)

        # Test plot methods
        #gs_bands.boxplot()

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
        """Testing abipy-->pymatgen converter"""
        nscf_bands = ElectronBands.from_file(data.ref_file("si_nscf_GSR.nc"))
        pmg_bands = nscf_bands.to_pymatgen()

    def test_derivatives(self):
        """Testing computation of effective masses."""
        ebands = ElectronBands.from_file(data.ref_file("si_nscf_GSR.nc"))

        # Hack eigens to simulate free-electron bands. This will produce all(effective masses == 1)
        new_eigens = np.empty(ebands.shape)
        branch = 0.5 * units.Ha_to_eV * np.array([(k.norm * units.bohr_to_ang)**2 for k in ebands.kpoints])
        for spin in ebands.spins:
            for band in range(ebands.mband):
                new_eigens[spin, :, band] = branch
        ebands._eigens = new_eigens

        effm_lines = ebands.effective_masses(spin=0, band=0, acc=2)

        # Flatten structure (.flatten does not work in this case)
        values = []
        for arr in effm_lines:
            values.extend(arr)
        self.assertArrayAlmostEqual(np.array(values), 1.0)
