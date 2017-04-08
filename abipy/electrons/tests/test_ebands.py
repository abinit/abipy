"""Tests for electrons.ebands module"""
from __future__ import print_function, division, unicode_literals

import sys
import numpy as np
import unittest
import abipy.data as abidata
import pymatgen.core.units as units

from abipy.core.kpoints import KpointList
from abipy.electrons.ebands import (ElectronBands, ElectronDos, ElectronBandsPlotter, ElectronDosPlotter,
    ElectronsReader, frame_from_ebands)
from abipy.core.testing import *


class EbandsReaderTest(AbipyTest):

    def test_reader(self):
        """Test ElectronsReader."""
        filepath = abidata.ref_file("si_scf_WFK.nc")

        with ElectronsReader(filepath) as r:
            kpoints = r.read_kpoints()
            assert isinstance(kpoints, KpointList)
            #assert len(kpoints) == ??
            #self.assert_all_equal(self.read_nband_sk == ??))

            eigens = r.read_eigenvalues()
            occfacts = r.read_occupations()
            fermie = r.read_fermie()
            self.assertTrue(r.read_nelect() == 8)
            #smearing = r.read_smearing()


class ElectronBandsTest(AbipyTest):

    def test_read_ebands_from_WFK(self):
        """Read ElectronBands from WFK files."""
        for ii, filename in enumerate(abidata.WFK_NCFILES):
            ebands = ElectronBands.from_file(filename)
            ebands.to_pymatgen()
            ebands.to_pdframe()
            assert ElectronBands.as_ebands(ebands) is ebands

            self.serialize_with_pickle(ebands, test_eq=False)
            ElectronBands.from_dict(ebands.as_dict())
            self.assertMSONable(ebands, test_if_subclass=False)

            if ii == 0:
                if self.has_matplotlib(): ebands.plot(show=False)
                ebands.to_xmgrace(self.get_tmpname(text=True))

    def test_read_ebands_from_GSR(self):
        """Read ElectronBands from GSR files."""
        for filename in abidata.GSR_NCFILES:
            ebands = ElectronBands.from_file(filename)
            ebands.to_pymatgen()
            ebands.to_pdframe()

            self.serialize_with_pickle(ebands, test_eq=False)
            ElectronBands.from_dict(ebands.as_dict())
            self.assertMSONable(ebands, test_if_subclass=False)

    def test_edos(self):
        """Test electron DOS methods."""
        gs_bands = ElectronBands.from_file(abidata.ref_file("si_scf_GSR.nc"))
        assert not gs_bands.has_metallic_scheme
        repr(gs_bands)
        str(gs_bands)
        assert gs_bands.to_string(title="Title", with_structure=False, with_kpoints=True, verbose=1)

        # Test get_e0
        assert gs_bands.get_e0("fermie") == gs_bands.fermie
        assert gs_bands.get_e0(None) == 0.0
        assert gs_bands.get_e0("None") == 0.0
        assert gs_bands.get_e0(1.0) == 1.0

        edos = gs_bands.get_edos()
        print(edos)
        assert ElectronDos.as_edos(edos, {}) is edos
        edos_samevals = ElectronDos.as_edos(gs_bands, {})
        assert ElectronDos.as_edos(gs_bands, {}) == edos
        assert ElectronDos.as_edos(abidata.ref_file("si_scf_GSR.nc"), {}) == edos

        mu = edos.find_mu(8)
        imu = edos.tot_idos.find_mesh_index(mu)
        self.assert_almost_equal(edos.tot_idos[imu][1], 8, decimal=2)

        d, i = edos.dos_idos(spin=0)
        tot_d, tot_i = edos.dos_idos()
        self.assert_almost_equal(2 * d.values, tot_d.values)
        self.assert_almost_equal(2 * i.values, tot_i.values)

        self.serialize_with_pickle(edos, protocols=[-1], test_eq=False)

        # Test plot methods
        if self.has_matplotlib():
            edos.plot(show=False)
            edos.plot_dos_idos(show=False)
            edos.plot_up_minus_down(show=False)
            gs_bands.plot_with_edos(dos=edos, show=False)
            if self.has_seaborn(): gs_bands.boxplot(show=False)

        if self.has_ipywidgets():
            assert gs_bands.ipw_edos_widget() is not None
        else:
            raise RuntimeError()

    def test_jdos(self):
        """Test JDOS methods."""
        bands = ElectronBands.from_file(abidata.ref_file("si_scf_GSR.nc"))

        spin = 0
        conduction = [4,]
        for v in range(1,5):
            valence = range(0, v)
            jdos = bands.get_ejdos(spin, valence, conduction)
            intg = jdos.integral()[-1][-1]
            self.assert_almost_equal(intg, len(conduction) * len(valence))

        self.serialize_with_pickle(jdos, protocols=[-1])

        nscf_bands = ElectronBands.from_file(abidata.ref_file("si_nscf_GSR.nc"))

        diffs = nscf_bands.statdiff(nscf_bands)
        assert diffs is not None
        print(diffs)

        # Test the detection of denerate states.
        degs = nscf_bands.degeneracies(spin=0, kpoint=[0,0,0], bands_range=range(8))

        ref_degbands = [[0], [1,2,3], [4,5,6], [7]]
        for i, (e, deg_bands) in enumerate(degs):
            self.assertEqual(deg_bands, ref_degbands[i])

        # Test Electron
        e1 = nscf_bands._electron_state(spin=0, kpoint=[0, 0, 0], band=0)
        str(e1)
        e1_copy = e1.copy()
        assert isinstance(e1.as_dict(), dict)
        assert isinstance(e1.to_strdict(), dict)
        assert e1.spin == 0
        assert e1.skb[0] == 0
        str(e1.tips)

        # JDOS requires a homogeneous sampling.
        with self.assertRaises(ValueError):
            nscf_bands.get_ejdos(spin, 0, 4)

    def test_ebands_skw_interpolation(self):
        if sys.version[0:3] >= '3.4':
            raise unittest.SkipTest(
                "SKW interpolation is not tested if Python version >= 3.4 (linalg.solve portability issue)"
             )

        bands = ElectronBands.from_file(abidata.ref_file("si_scf_GSR.nc"))

        # Test interpolate.
        vertices_names = [((0.0, 0.0, 0.0), "G"), ((0.5, 0.5, 0.0), "M")]
        r = bands.interpolate(lpratio=10, vertices_names=vertices_names, kmesh=[8, 8, 8], verbose=1)
        assert r.ebands_kpath is not None
        assert r.ebands_kpath.kpoints.is_path
        assert not r.ebands_kpath.kpoints.is_ibz
        mpdivs, shifts = r.ebands_kpath.kpoints.mpdivs_shifts
        assert mpdivs is None and shifts is None

        assert r.ebands_kmesh is not None
        assert r.ebands_kmesh.kpoints.is_ibz
        assert not r.ebands_kmesh.kpoints.is_path
        assert r.ebands_kmesh.kpoints.ksampling is not None
        assert r.ebands_kmesh.kpoints.is_mpmesh
        mpdivs, shifts = r.ebands_kmesh.kpoints.mpdivs_shifts
        self.assert_equal(mpdivs, [8, 8, 8])
        self.assert_equal(shifts.flatten(), [0, 0, 0])

        # Export it in BXSF format.
        r.ebands_kmesh.to_bxsf(self.get_tmpname(text=True))

    def test_pymatgen_converter(self):
        """Testing abipy-->pymatgen converter"""
        nscf_bands = ElectronBands.from_file(abidata.ref_file("si_nscf_GSR.nc"))
        pmg_bands = nscf_bands.to_pymatgen()

    def test_derivatives(self):
        """Testing computation of effective masses."""
        ebands = ElectronBands.from_file(abidata.ref_file("si_nscf_GSR.nc"))

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

    def test_to_bxsf(self):
        """Testing Fermi surface exporter."""
        from abipy.abilab import abiopen
        with abiopen(abidata.ref_file("mgb2_kmesh181818_FATBANDS.nc")) as fbnc_kmesh:
            fbnc_kmesh.ebands.to_bxsf(self.get_tmpname(text=True))


class FrameFromEbandsTest(AbipyTest):

    def test_frame_from_ebands(self):
        """Testing frame_from_ebands."""
        gsr_scf_path = abidata.ref_file("si_scf_GSR.nc")
        gs_ebands = ElectronBands.as_ebands(gsr_scf_path)
        gsr_nscf_path = abidata.ref_file("si_nscf_GSR.nc")
        index = ["foo", "bar", "hello"]
        df = frame_from_ebands([gsr_scf_path, gs_ebands, gsr_nscf_path], index=index, with_spglib=True)
        #print(df)
        assert all(f == "Si2" for f in df["formula"])
        assert all(num == 227 for num in df["abispg_num"])
        assert all(df["spglib_num"] == df["abispg_num"])


class ElectronBandsPlotterTest(AbipyTest):

    def test_api(self):
        """Test ElelectronBandsPlotter API."""
        plotter = ElectronBandsPlotter(key_ebands=[("Si1", abidata.ref_file("si_scf_GSR.nc"))])
        plotter.add_ebands("Si2", abidata.ref_file("si_scf_GSR.nc"))
        str(plotter)
        repr(plotter)

        assert len(plotter.ebands_list) == 2
        assert len(plotter.edoses_list) == 0
        with self.assertRaises(ValueError):
            plotter.add_ebands("Si2", abidata.ref_file("si_scf_GSR.nc"))

        print(plotter.bands_statdiff())
        df = plotter.get_ebands_frame()
        assert df is not None

        if self.has_matplotlib():
            plotter.combiplot(title="Silicon band structure", show=False)
            if self.has_seaborn():
                plotter.combiboxplot(title="Silicon band structure", show=False)
            plotter.gridplot(title="Silicon band structure", show=False)
            plotter.boxplot(title="Silicon band structure", swarm=True, show=False)
            plotter.animate(show=False)

        if self.has_ipywidgets():
            assert plotter.ipw_select_plot() is not None

        if self.has_nbformat():
            plotter.write_notebook(nbpath=self.get_tmpname(text=True))

        pickle_path = self.get_tmpname(text=True)
        plotter.pickle_dump(pickle_path)
        same = ElectronBandsPlotter.pickle_load(pickle_path)
        assert len(same.ebands_dict) == len(plotter.ebands_dict)
        assert list(same.ebands_dict.keys()) == list(plotter.ebands_dict.keys())


class ElectronDosPlotterTest(AbipyTest):

    def test_api(self):
        """Test ElelectronDosPlotter API."""
        gsr_path = abidata.ref_file("si_scf_GSR.nc")
        gs_bands = ElectronBands.from_file(gsr_path)
        edos = gs_bands.get_edos()

        plotter = ElectronDosPlotter()
        plotter.add_edos("edos1", edos)
        plotter.add_edos("edos2", gsr_path, edos_kwargs=dict(method="gaussian", step=0.2, width=0.4))
        assert len(plotter.edos_list) == 2
        assert not plotter._can_use_basenames_as_labels()

        if self.has_matplotlib():
            plotter.combiplot(show=False)
            plotter.gridplot(show=False)

        if self.has_ipywidgets():
            assert plotter.ipw_select_plot() is not None

        if self.has_nbformat():
            plotter.write_notebook(nbpath=self.get_tmpname(text=True))
