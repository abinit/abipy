"""Tests for electrons.ebands module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import numpy as np
import unittest
import pymatgen.core.units as units
import abipy.data as abidata
from abipy import abilab

from abipy.electrons.ebands import (ElectronBands, ElectronDos, ElectronBandsPlotter, ElectronDosPlotter,
    ElectronsReader, dataframe_from_ebands, Smearing)
from abipy.core.testing import AbipyTest


class SmearingTest(AbipyTest):
    def test_smearing_info(self):
        """Testing SmearingInfo."""
        with self.assertRaises(ValueError):
            Smearing(scheme=None, occopt=1)

        sm = Smearing(scheme=None, occopt=3, tsmear_ev=0.0)
        repr(sm); str(sm)
        self.assertMSONable(sm, test_if_subclass=False)
        assert sm.has_metallic_scheme
        assert Smearing.as_smearing(sm) is sm
        assert Smearing.as_smearing(sm).occopt == 3
        assert Smearing.as_smearing(None).occopt == 1


class EbandsReaderTest(AbipyTest):

    def test_reader(self):
        """Testing ElectronsReader with WFK file."""

        with ElectronsReader(abidata.ref_file("si_scf_WFK.nc")) as r:
            nsppol = r.read_nsppol()
            nspden = r.read_nspden()
            nspinor = r.read_nspinor()
            assert nsppol == 1 and nspden == 1 and nspinor == 1
            kpoints = r.read_kpoints()
            assert kpoints.is_ibz and not kpoints.is_path
            assert kpoints.weights.sum() == 1
            nkpt = len(kpoints)
            assert nkpt == 29
            assert kpoints.to_array().shape == (nkpt, 3)
            assert kpoints.ksampling.kptopt == 1

            mband = 8
            nband_sk = r.read_nband_sk()
            self.assert_equal(nband_sk, mband)
            assert nband_sk.shape == (1, nkpt)

            eigens = r.read_eigenvalues()
            assert eigens.shape == (nsppol, nkpt, mband)
            assert str(eigens.unit) == "eV"
            occfacts = r.read_occupations()
            assert occfacts.shape == (nsppol, nkpt, mband)
            fermie = r.read_fermie()
            self.assert_almost_equal(fermie.to("Ha"), 0.205739364929578)
            assert r.read_nelect() == 8
            smearing = r.read_smearing()
            repr(smearing); str(smearing)
            assert smearing.occopt == 1
            self.assert_almost_equal(smearing.tsmear_ev.to("Ha"), 0.01)
            assert not smearing.has_metallic_scheme
            assert smearing.scheme == "none"

            self.assertMSONable(smearing, test_if_subclass=False)
            assert len(smearing.to_json())


class ElectronBandsTest(AbipyTest):

    def test_nickel_ebands_spin(self):
        """Testing Nickel electron bands with nsppol == 2"""
        ref_nelect = 18
        ni_ebands_kmesh = ElectronBands.from_file(abidata.ref_file("ni_666k_GSR.nc"))
        assert ElectronBands.as_ebands(ni_ebands_kmesh) is ni_ebands_kmesh
        with self.assertRaises(TypeError):
            ElectronBands.as_ebands(1)

        repr(ni_ebands_kmesh); str(ni_ebands_kmesh)
        assert ni_ebands_kmesh.nsppol == 2 and ni_ebands_kmesh.nspinor == 1 and ni_ebands_kmesh.nspden == 2
        assert ni_ebands_kmesh.nelect == ref_nelect
        assert ni_ebands_kmesh.kpoints.is_ibz and ni_ebands_kmesh.has_bzmesh and not ni_ebands_kmesh.has_bzpath
        assert ni_ebands_kmesh.has_timrev
        assert ni_ebands_kmesh.has_metallic_scheme
        smearing = ni_ebands_kmesh.smearing
        assert smearing.has_metallic_scheme
        assert smearing.occopt == 7
        self.assert_almost_equal(smearing.tsmear_ev.to("Ha"), 0.0075)
        assert smearing.scheme == "gaussian"
        assert not ni_ebands_kmesh.get_gaps_string()

        #ni_ebands_kmesh.copy()
        #ni_ebands_kmesh.deepcopy()

        ni_edos = ni_ebands_kmesh.get_edos()
        repr(ni_edos); str(ni_edos)
        assert ni_edos.to_string(verbose=2)

        # Get ElectronDosPlotter with nsppol == 2 and test matplotlib methods.
        edos_plotter = ni_ebands_kmesh.compare_gauss_edos(widths=[0.2, 0.4], step=0.2)
        assert len(edos_plotter) == 2
        if self.has_matplotlib():
            # Combiplot.
            assert edos_plotter.combiplot(title="default values", show=False)
            assert edos_plotter.combiplot(what_list=("dos", "idos"), spin_mode="resolved", show=False)
            assert edos_plotter.combiplot(e0=0, what_list="dos", spin_mode="automatic", fontsize=12, show=False)
            # Gridplot
            assert edos_plotter.gridplot(title="default values", show=False)
            assert edos_plotter.gridplot(what="idos", spin_mode="resolved", xlims=(-10, 10), show=False)
            assert edos_plotter.gridplot(e0=0, what="dos", spin_mode="resolved", fontsize=12, show=False)

        ni_ebands_kpath = ElectronBands.from_file(abidata.ref_file("ni_kpath_GSR.nc"))
        assert not ni_ebands_kpath.has_linewidths
        ni_ebands_kpath.linewidths = np.ones(ni_ebands_kpath.shape)
        assert ni_ebands_kpath.has_linewidths

        repr(ni_ebands_kpath); str(ni_ebands_kpath)
        assert ni_ebands_kpath.nsppol == 2 and ni_ebands_kpath.nspinor == 1 and ni_ebands_kpath.nspden == 2
        assert ni_ebands_kpath.nelect == ref_nelect
        assert ni_ebands_kpath.kpoints.is_path and not ni_ebands_kpath.has_bzmesh and ni_ebands_kpath.has_bzpath
        assert ni_ebands_kpath.has_timrev
        assert ni_ebands_kpath.fermie == ni_ebands_kmesh.fermie

        # Serialization
        self.serialize_with_pickle(ni_ebands_kpath, test_eq=False)
        self.assertMSONable(ni_ebands_kpath, test_if_subclass=False)
        assert len(ni_ebands_kpath.to_json())

        od = ni_ebands_kmesh.get_dict4pandas(with_spglib=False)
        assert od["nsppol"] == 2 and od["nspinor"] == 1 and od["nspden"] == 2

        df = ni_ebands_kpath.get_dataframe()
        ni_ebands_kpath.to_xmgrace(self.get_tmpname(text=True))
        ni_ebands_kpath.to_xmgrace(sys.stdout)

        # BXSF cannot be produced because.
        #ngkpt    6 6 6
        #nshiftk  4
        #shiftk   1/2 1/2 1/2 1/2 0.0 0.0 0.0 1/2 0.0 0.0 0.0 1/2
        with self.assertRaises(ValueError):
            ni_ebands_kmesh.to_bxsf(self.get_tmpname(text=True))

        # Interpolation
        r = ni_ebands_kmesh.interpolate(lpratio=10, kmesh=[8, 8, 8], verbose=1)

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

        # Test __add__ and __radd__ (should return ElectronBandsPlotter)
        p = ni_ebands_kmesh + r.ebands_kmesh + r.ebands_kpath
        assert hasattr(p, "combiplot")

        # Test plot methods
        if self.has_matplotlib():
            elims = [-10, 2]
            assert ni_ebands_kmesh.plot(show=False)
            assert ni_ebands_kmesh.plot_bz(show=False)
            assert ni_ebands_kpath.plot(ylims=elims, with_gaps=True, show=False)
            assert ni_ebands_kpath.plot_with_edos(ni_edos, ylims=elims, show=False)
            assert ni_ebands_kpath.plot_bz(show=False)
            assert ni_ebands_kpath.plot_transitions(4.4, qpt=(0, 0, 0), atol_ev=0.1, atol_kdiff=1e-4, show=False)
            assert ni_ebands_kpath.plot_transitions(4.4, qpt=(0.03, 0, 0), atol_ev=0.5, atol_kdiff=0.2, show=False)
            assert ni_ebands_kpath.plot_scatter3d(spin=0, band=8, show=False)
            assert ni_edos.plot(xlims=elims, show=False)
            assert ni_edos.plot_dos_idos(xlims=elims, show=False)
            assert ni_edos.plot_up_minus_down(xlims=elims, show=False)

            # Test linewidths
            assert ni_ebands_kmesh.plot_lws_vs_e0(show=False) is None
            assert ni_ebands_kpath.plot_lws_vs_e0(show=False)

            # TODO Generaliza jdos to metals.
            #vrange, crange = range(0, 4), range(4, 5)
            #assert ni_ebands_kmesh.plot_ejdosvc(vrange, crange, cumulative=False, show=False)
            #assert ni_ebands_kmesh.plot_ejdosvc(vrange, crange, cumulative=True, show=False)

            if self.has_seaborn():
                assert ni_ebands_kmesh.boxplot(brange=[5, 10], show=False,
                    title="Boxplot for up and down spin and 10 > band >= 5")

        # Test Abipy --> Pymatgen converter.
        pmg_bands_kpath = ni_ebands_kpath.to_pymatgen()
        assert hasattr(pmg_bands_kpath, "get_branch")   # Should be BandStructureSymmLine
        assert pmg_bands_kpath.efermi == ni_ebands_kpath.fermie
        assert pmg_bands_kpath.is_spin_polarized
        assert pmg_bands_kpath.is_metal()

        # Test Pymatgen --> Abipy converter.
        same_ekpath = ElectronBands.from_pymatgen(pmg_bands_kpath, ni_ebands_kpath.nelect)
        repr(same_ekpath); str(same_ekpath)
        self.assert_equal(same_ekpath.eigens, ni_ebands_kpath.eigens)
        assert same_ekpath.fermie == ni_ebands_kpath.fermie

        pmg_bands_kmesh = ni_ebands_kmesh.to_pymatgen()
        #assert hasattr(pmg_bands_kmesh, "get_branch")   # Should be BandStructure
        assert pmg_bands_kmesh.efermi == ni_ebands_kmesh.fermie
        assert pmg_bands_kmesh.is_spin_polarized
        assert pmg_bands_kmesh.is_metal()

        # Test Pymatgen --> Abipy converter.
        same_ekmesh = ElectronBands.from_pymatgen(pmg_bands_kmesh, ni_ebands_kmesh.nelect)
        self.assert_equal(same_ekmesh.eigens, ni_ebands_kmesh.eigens)
        assert same_ekmesh.fermie == ni_ebands_kmesh.fermie

    def test_silicon_ebands(self):
        """Testing electron bands with nsppol == 1"""
        si_ebands_kmesh = ElectronBands.from_file(abidata.ref_file("si_scf_GSR.nc"))
        assert not si_ebands_kmesh.has_metallic_scheme
        repr(si_ebands_kmesh); str(si_ebands_kmesh)
        assert si_ebands_kmesh.to_string(title="Title",
                with_structure=False, with_kpoints=True, verbose=1)

        for spin, ik, band in si_ebands_kmesh.skb_iter():
            assert spin == 0
            assert si_ebands_kmesh.nkpt >= ik >= 0
            assert si_ebands_kmesh.nband_sk[spin, ik] >= band >= 0

        # Test ElectronBands get_e0
        assert si_ebands_kmesh.get_e0("fermie") == si_ebands_kmesh.fermie
        assert si_ebands_kmesh.get_e0(None) == 0.0
        assert si_ebands_kmesh.get_e0("None") == 0.0
        assert si_ebands_kmesh.get_e0(1.0) == 1.0
        with self.assertRaises(ValueError):
            si_ebands_kmesh.get_e0("foo")

        r = si_ebands_kmesh.with_points_along_path(knames=["G", "X", "L", "G"])
        assert r.ebands.kpoints.is_path
        assert r.ebands.kpoints[0].is_gamma
        assert r.ebands.kpoints[0] == r.ebands.kpoints[-1]
        assert si_ebands_kmesh.kpoints[r.ik_new2prev[0]].is_gamma

        # Serialization
        self.serialize_with_pickle(si_ebands_kmesh, test_eq=False)
        self.assertMSONable(si_ebands_kmesh, test_if_subclass=False)
        assert len(si_ebands_kmesh.to_json())

        dless_states = si_ebands_kmesh.dispersionless_states()
        assert not dless_states

        estats = si_ebands_kmesh.spacing()
        self.assert_almost_equal(estats.mean, 2.3100587301616917)
        self.assert_almost_equal(estats.stdev, 2.164400652355628)
        self.assert_almost_equal(estats.min, 0)
        self.assert_almost_equal(estats.max, 11.855874158768694)
        repr(estats); str(estats)

        with self.assertRaises(NotImplementedError):
            si_ebands_kmesh.get_edos(method="tetrahedron")

        si_edos = si_ebands_kmesh.get_edos()
        repr(si_edos); str(si_edos)
        assert ElectronDos.as_edos(si_edos, {}) is si_edos
        assert si_edos == si_edos and not (si_edos != si_edos)
        edos_samevals = ElectronDos.as_edos(si_ebands_kmesh, {})
        assert ElectronDos.as_edos(si_ebands_kmesh, {}) == si_edos
        assert ElectronDos.as_edos(abidata.ref_file("si_scf_GSR.nc"), {}) == si_edos
        with self.assertRaises(TypeError):
            ElectronDos.as_edos({}, {})

        mu = si_edos.find_mu(8)
        imu = si_edos.tot_idos.find_mesh_index(mu)
        self.assert_almost_equal(si_edos.tot_idos[imu][1], 8, decimal=2)

        d, i = si_edos.dos_idos(spin=0)
        tot_d, tot_i = si_edos.dos_idos()
        self.assert_almost_equal(2 * d.values, tot_d.values)
        self.assert_almost_equal(2 * i.values, tot_i.values)
        self.assert_equal(si_edos.up_minus_down.mesh, si_edos.tot_dos.mesh)
        self.assert_equal(si_edos.up_minus_down.values, 0)

        # Test ElectronDos get_e0
        assert si_edos.get_e0("fermie") == si_edos.fermie
        assert si_edos.get_e0(None) == 0.0
        assert si_edos.get_e0("None") == 0.0
        assert si_edos.get_e0(1.0) == 1.0
        with self.assertRaises(TypeError):
            si_edos.get_e0("foo")

        self.serialize_with_pickle(si_edos, protocols=[-1], test_eq=False)

        # Test plot methods
        if self.has_matplotlib():
            klabels = {
                (0,0,0): r"$\Gamma$",
                (0.375, 0.375, 0.7500): "K",
                (0.5, 0.5, 1.0): "X",
                (0.5, 0.5, 0.5): "L",
                (0.5, 0.0, 0.5): "X",
                (0.5, 0.25, 0.75): "W",
            }

            assert si_edos.plot(show=False)
            assert si_edos.plot_dos_idos(show=False)
            assert si_edos.plot_up_minus_down(show=False)
            assert si_ebands_kmesh.plot_with_edos(edos=si_edos, klabels=klabels, with_gaps=True, show=False)
            assert si_ebands_kmesh.kpoints.plot(show=False)

            vrange, crange = range(0, 4), range(4, 5)
            assert si_ebands_kmesh.plot_ejdosvc(vrange, crange, cumulative=False, show=False)
            assert si_ebands_kmesh.plot_ejdosvc(vrange, crange, cumulative=True, show=False)
            assert si_ebands_kmesh.kpoints.plot(show=False)
            if self.has_seaborn():
                assert si_ebands_kmesh.boxplot(swarm=True, show=False)

        if self.has_ipywidgets():
            assert si_ebands_kmesh.ipw_edos_widget()

        # Test Abipy --> Pymatgen converter.
        pmg_bands_kmesh = si_ebands_kmesh.to_pymatgen()
        assert pmg_bands_kmesh.efermi == si_ebands_kmesh.fermie
        assert not pmg_bands_kmesh.is_spin_polarized
        #assert not pmg_bands_kmesh.is_metal()

        # Test Pymatgen --> Abipy converter.
        same_ekmesh = ElectronBands.from_pymatgen(pmg_bands_kmesh, si_ebands_kmesh.nelect)
        repr(same_ekmesh); str(same_ekmesh)
        self.assert_equal(same_ekmesh.eigens, si_ebands_kmesh.eigens)
        assert same_ekmesh.fermie == si_ebands_kmesh.fermie
        assert len(same_ekmesh.kpoints) == len(pmg_bands_kmesh.kpoints)

        # Test JDOS methods.
        spin = 0
        conduction = [4,]
        for v in range(1, 5):
            valence = range(0, v)
            jdos = si_ebands_kmesh.get_ejdos(spin, valence, conduction)
            intg = jdos.integral()[-1][-1]
            self.assert_almost_equal(intg, len(conduction) * len(valence))

        self.serialize_with_pickle(jdos, protocols=[-1])

        si_ebands_kpath = ElectronBands.from_file(abidata.ref_file("si_nscf_GSR.nc"))

        diffs = si_ebands_kpath.statdiff(si_ebands_kpath)
        assert diffs is not None
        repr(diffs); str(diffs)

        homo = si_ebands_kpath.homos[0]
        repr(homo); str(homo)
        assert homo.spin == 0 and homo.occ == 2.0 and homo.band == 3
        assert homo.kpoint == [0, 0, 0]
        assert homo == homo
        assert si_ebands_kpath.kpoints[homo.kidx] == homo.kpoint
        self.assert_almost_equal(homo.eig, 5.5983129712050665)
        assert "eig" in homo.__class__.get_fields()

        lumo = si_ebands_kpath.lumos[0]
        assert homo != lumo
        assert lumo.spin == 0 and lumo.occ == 0.0 and lumo.band == 4
        self.assert_almost_equal(lumo.kpoint.frac_coords, [0.,  0.4285714, 0.4285714])
        assert si_ebands_kpath.kpoints[lumo.kidx] == lumo.kpoint
        self.assert_almost_equal(lumo.eig, 6.1226526474610843)

        dir_gap = si_ebands_kpath.direct_gaps[0]
        fun_gap = si_ebands_kpath.fundamental_gaps[0]
        assert fun_gap.energy <= dir_gap.energy
        assert dir_gap.qpoint == [0, 0, 0]
        assert dir_gap.is_direct
        assert dir_gap == dir_gap
        assert dir_gap != fun_gap
        assert si_ebands_kpath.get_gaps_string()
        #print("repr_fun_gap", repr(fun_gap), id(fun_gap), id(fun_gap.qpoint))
        #print("repr_dir_gap", repr(dir_gap), id(dir_gap), id(dir_gap.qpoint))
        self.assert_almost_equal(dir_gap.energy, 2.5318279814319133)
        self.assert_almost_equal(fun_gap.qpoint.frac_coords, [0.,  0.4285714, 0.4285714])
        self.assert_almost_equal(fun_gap.energy, 0.52433967625601774)
        assert not fun_gap.is_direct

        e1 = si_ebands_kpath.lomo_sk(spin=0, kpoint=0)
        e2 = si_ebands_kpath.lomo_sk(spin=0, kpoint=si_ebands_kpath.kpoints[0])
        assert e1.eig == e2.eig

        # Test abipy-->pymatgen converter
        pmg_bands_kpath = si_ebands_kpath.to_pymatgen()
        assert hasattr(pmg_bands_kpath, "get_branch")  # Should be BandStructureSymmLine
        assert pmg_bands_kpath.efermi == si_ebands_kpath.fermie
        assert not pmg_bands_kpath.is_spin_polarized
        assert not pmg_bands_kpath.is_metal()

        # Test the detection of denerate states.
        degs = si_ebands_kpath.degeneracies(spin=0, kpoint=[0, 0, 0], bands_range=range(8))
        ref_degbands = [[0], [1, 2, 3], [4, 5, 6], [7]]
        for i, (e, deg_bands) in enumerate(degs):
            self.assertEqual(deg_bands, ref_degbands[i])

        # Test Electron
        e1 = si_ebands_kpath._electron_state(spin=0, kpoint=[0, 0, 0], band=0)
        repr(e1); str(e1)
        e1_copy = e1.copy()
        assert isinstance(e1.as_dict(), dict)
        assert isinstance(e1.to_strdict(), dict)
        assert e1.spin == 0
        assert e1.skb[0] == 0
        str(e1.tips)

        # JDOS requires a homogeneous sampling.
        with self.assertRaises(ValueError):
            si_ebands_kpath.get_ejdos(spin, 0, 4)

        if self.has_matplotlib():
            assert si_ebands_kpath.plot(spin=0, e0=0, with_gaps=True, max_phfreq=0.1, show=False)

    def test_ebands_skw_interpolation(self):
        """Testing SKW interpolation."""
        si_ebands_kmesh = ElectronBands.from_file(abidata.ref_file("si_scf_GSR.nc"))

        # Test interpolation.
        vertices_names = [((0.0, 0.0, 0.0), "G"), ((0.5, 0.5, 0.0), "M")]
        r = si_ebands_kmesh.interpolate(lpratio=10, vertices_names=vertices_names,
                                        kmesh=[8, 8, 8], verbose=1)
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

    def test_derivatives(self):
        """Testing computation of effective masses."""
        ebands = ElectronBands.from_file(abidata.ref_file("si_nscf_GSR.nc"))

        # Hack eigens to simulate free-electron bands.
        # This should produce all(effective masses == 1)
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

        self.assert_almost_equal(np.array(values), 1.0)

        em = ebands.effmass_line(spin=0, kpoint=(0, 0, 0), band=0)
        repr(em); str(em)
        #self.assert_almost_equal(np.array(values), 1.0)

    def test_fermi_surface(self):
        """Testing Fermi surface tools."""
        with abilab.abiopen(abidata.ref_file("mgb2_kmesh181818_FATBANDS.nc")) as fbnc_kmesh:
            ebands = fbnc_kmesh.ebands
            str(ebands)
            ebands.to_bxsf(self.get_tmpname(text=True))

            # Test Ebands3d
            eb3d = ebands.get_ebands3d()
            repr(eb3d); str(eb3d)
            assert eb3d.to_string(verbose=2)

            if self.has_matplotlib():
                assert eb3d.plot_contour(band=4, spin=0, plane="xy", elevation=0, show=False)
                if self.has_skimage():
                    assert eb3d.plot_isosurfaces(e0="fermie", cmap="rainbow", verbose=1, show=False)

            # Test Mayavi
            if self.has_mayavi():
                assert eb3d.mvplot_isosurfaces(verbose=1, show=False)
                #assert eb3d.mvplot_cutplanes(band=4, spin=0, show=False)

    def test_frame_from_ebands(self):
        """Testing dataframe_from_ebands."""
        gsr_kmesh = abidata.ref_file("si_scf_GSR.nc")
        si_ebands_kmesh = ElectronBands.as_ebands(gsr_kmesh)
        gsr_nscf_path = abidata.ref_file("si_nscf_GSR.nc")
        index = ["foo", "bar", "hello"]
        df = dataframe_from_ebands([gsr_kmesh, si_ebands_kmesh, gsr_nscf_path], index=index, with_spglib=True)
        #str(df)
        assert all(f == "Si2" for f in df["formula"])
        assert all(num == 227 for num in df["abispg_num"])
        assert all(df["spglib_num"] == df["abispg_num"])


class ElectronBandsFromRestApi(AbipyTest):

    def test_from_mpid(self):
        """Testing interpolation of SnO2 band energies from MP database."""
        #mpid = "mp-149"
        #mpid = "mp-856"
        mpid = "mp-3079"
        ebands = abilab.ElectronBands.from_mpid(mpid)
        # Use prune_step to remove k-points (too many k-points on a k-path can cause numerical instabilities)
        ebands = ebands.new_with_irred_kpoints(prune_step=2)
        # Interpolate on k-path + kmesh.
        # Results are very sensitive to the value of lpratio. The default is not enough in this case!!
        r = ebands.interpolate(lpratio=50, kmesh=[10, 10, 10])

        new_fermie = r.ebands_kpath.set_fermie_to_vbm()
        assert new_fermie == r.ebands_kpath.fermie

        edos = r.ebands_kmesh.get_edos()
        new_fermie = r.ebands_kpath.set_fermie_from_edos(edos)
        assert new_fermie == edos.fermie

        if self.has_matplotlib():
            # Plot bands + dos using interpolated energies.
            assert r.ebands_kpath.plot_with_edos(edos, show=False)


class ElectronBandsPlotterTest(AbipyTest):

    def test_ebands_plotter(self):
        """Testing ElelectronBandsPlotter."""
        plotter = ElectronBandsPlotter(key_ebands=[("Si1", abidata.ref_file("si_scf_GSR.nc"))])
        plotter.add_ebands("Si2", abidata.ref_file("si_scf_GSR.nc"))
        repr(plotter); str(plotter)

        assert len(plotter.ebands_list) == 2
        assert len(plotter.edoses_list) == 0
        with self.assertRaises(ValueError):
            plotter.add_ebands("Si2", abidata.ref_file("si_scf_GSR.nc"))

        # __add__ should merge plotters
        p2 = plotter.add_plotter(plotter)
        assert len(p2.ebands_list) == 2
        assert len(p2.edoses_list) == 0
        assert hasattr(p2, "combiplot")

        print(plotter.bands_statdiff())
        df = plotter.get_ebands_frame()
        assert df is not None

        if self.has_matplotlib():
            assert plotter.combiplot(title="Silicon band structure", show=False)
            # Alias for combiplot
            assert plotter.plot(e0=2, width_ratios=(3, 1), fontsize=12, show=False)
            if self.has_seaborn():
                plotter.combiboxplot(title="Silicon band structure", swarm=True, show=False)
            assert plotter.gridplot(title="Silicon band structure", with_gaps=True, show=False)
            assert plotter.boxplot(title="Silicon band structure", swarm=True, show=False)
            assert plotter.plot_band_edges(epad_ev=2.0, show=False)
            assert plotter.animate(show=False)

        if self.has_ipywidgets():
            assert plotter.ipw_select_plot() is not None

        if self.has_nbformat():
            assert plotter.write_notebook(nbpath=self.get_tmpname(text=True))

        pickle_path = self.get_tmpname(text=True)
        plotter.pickle_dump(pickle_path)
        same = ElectronBandsPlotter.pickle_load(pickle_path)
        assert len(same.ebands_dict) == len(plotter.ebands_dict)
        assert list(same.ebands_dict.keys()) == list(plotter.ebands_dict.keys())


class ElectronDosPlotterTest(AbipyTest):

    def test_api(self):
        """Testing ElelectronDosPlotter API."""
        gsr_path = abidata.ref_file("si_scf_GSR.nc")
        gs_bands = ElectronBands.from_file(gsr_path)
        si_edos = gs_bands.get_edos()

        plotter = ElectronDosPlotter()
        plotter.add_edos("edos1", si_edos)
        with self.assertRaises(ValueError):
            plotter.add_edos("edos1", si_edos)
        plotter.add_edos("edos2", gsr_path, edos_kwargs=dict(method="gaussian", step=0.2, width=0.4))
        assert len(plotter.edos_list) == 2
        assert not plotter._can_use_basenames_as_labels()

        if self.has_matplotlib():
            assert plotter.combiplot(show=False)
            assert plotter.gridplot(show=False)

        if self.has_ipywidgets():
            assert plotter.ipw_select_plot() is not None

        if self.has_nbformat():
            assert plotter.write_notebook(nbpath=self.get_tmpname(text=True))
