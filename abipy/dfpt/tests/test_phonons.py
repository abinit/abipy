"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import unittest
import sys
import os
import pickle
import numpy as np
import abipy.data as abidata
import abipy.core.abinit_units as abu

from abipy import abilab
from abipy.dfpt.phonons import (PhononBands, PhononDos, PhdosFile, phbands_gridplot,
        PhononBandsPlotter, PhononDosPlotter, dataframe_from_phbands)
from abipy.dfpt.ddb import DdbFile
from abipy.core.testing import AbipyTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


class PhononBandsTest(AbipyTest):

    def test_base(self):
        """Base tests for PhononBands"""
        filename = abidata.ref_file("trf2_5.out_PHBST.nc")
        phbands = PhononBands.from_file(filename)
        repr(phbands); str(phbands)
        assert phbands.to_string(title="Title", with_structure=False, with_qpoints=True, verbose=1)

        assert PhononBands.as_phbands(phbands) is phbands
        with self.assertRaises(TypeError):
            PhononBands.as_phbands({})
        assert np.array_equal(PhononBands.as_phbands(filename).phfreqs, phbands.phfreqs)

        with abilab.abiopen(abidata.ref_file("trf2_5.out_PHBST.nc")) as nc:
            repr(nc); str(nc)
            assert nc.to_string(verbose=1)
            assert nc.params["nqpt"] == len(nc.qpoints)
            assert nc.qpoints.is_path
            assert any(q.name is not None for q in nc.qpoints)
            same_phbands_nc = PhononBands.as_phbands(nc)
            self.assert_equal(same_phbands_nc.phfreqs, phbands.phfreqs)
            assert phbands.phdispl_cart.shape == (phbands.nqpt, phbands.num_branches, phbands.num_branches)
            # a + b gives plotter
            assert hasattr(same_phbands_nc + phbands, "combiplot")
            assert phbands.epsinf is None and phbands.zcart is None

        self.serialize_with_pickle(phbands, protocols=[-1], test_eq=False)

        # From pickle file.
        tmp_path = self.get_tmpname(suffix=".pickle")
        with open(tmp_path, "wb") as fh:
            pickle.dump(phbands, fh)
        same_phbands = PhononBands.as_phbands(tmp_path)
        self.assert_equal(same_phbands.phfreqs, phbands.phfreqs)

        # a + b + c gives plotter
        p = phbands + same_phbands + same_phbands_nc
        assert hasattr(p, "combiplot")

        assert phbands.minfreq == 0.0
        #self.assertEqual(phbands.maxfreq, 30)
        assert phbands.phfactor_ev2units("eV") == abu.phfactor_ev2units("eV")

        # Test XYZ vib
        phbands.create_xyz_vib(iqpt=0, filename=self.get_tmpname(text=True), max_supercell=[4, 4, 4])
        # Test ascii file
        phbands.create_ascii_vib(iqpts=0, filename=self.get_tmpname(text=True), pre_factor=1)
        # Test phononwebsite file
        phbands.create_phononwebsite_json(filename=self.get_tmpname(text=True), name='test')
        assert phbands.view_phononwebsite(verbose=1, dryrun=True) == 0
        # Test xmgrace
        phbands.to_xmgrace(self.get_tmpname(text=True))
        #phbands.to_xmgrace(sys.stdout)

        df = phbands.get_dataframe()
        assert "freq" in df and "mode" in df
        self.assert_almost_equal(df["freq"].values.min(), 0)

        umodes = phbands.get_unstable_modes(below_mev=-1000)
        assert len(umodes) == 0

        acoustic_modes = phbands.acoustic_indices((0, 0, 0))
        self.assertArrayEqual(acoustic_modes, [0, 1, 2])
        asr_breaking = phbands.asr_breaking()
        assert asr_breaking.absmax_break == 0

        # Test convertion to eigenvectors. Verify that they are orthonormal
        # Allow relatively large tolerance due to possible mismatching in the atomic masses between abinit and pmg
        # (Note that amu is None here)
        assert phbands.amu is None
        eig = phbands.dyn_mat_eigenvect
        assert len(eig) == phbands.nqpt

        cidentity = np.eye(len(eig[0]), dtype=np.complex)
        for iq in range(len(eig)):
            #print("About to test iq", iq, np.dot(eig[iq], eig[iq].T))
            #assert np.allclose(np.dot(eig[iq], eig[iq].T), cidentity , atol=1e-5, rtol=1e-3)
            assert np.allclose(np.dot(eig[iq].conjugate().T, eig[iq]), cidentity , atol=1e-5, rtol=1e-3)
            #self.assert_almost_equal(np.dot(eig[iq].conjugate().T, eig[iq]), cidentity)

        # Mapping reduced coordinates -> labels
        qlabels = {
            (0,0,0): r"$\Gamma$",
            (0.375, 0.375, 0.75): "K",
            (0.5, 0.5, 1.0): "X",
            (0.5, 0.5, 0.5): "L",
            (0.5, 0.0, 0.5): "X",
            (0.5, 0.25, 0.75): "W",
        }

        if self.has_matplotlib():
            assert phbands.plot(units="Thz", show=False, temp=300)
            assert phbands.plot_fatbands(units="ha", qlabels=qlabels, show=False)
            assert phbands.plot_fatbands(phdos_file=abidata.ref_file("trf2_5.out_PHDOS.nc"), units="thz", show=False)
            assert phbands.plot_colored_matched(units="cm^-1", show=False)
            assert phbands.plot_phdispl(qpoint=(0, 0, 0), units="cm^-1", hatches=None, show=False)
            assert phbands.plot_phdispl(qpoint=(0, 0, 0), units="cm^-1", hatches=None, show=False, cart_dir="x+y")
            assert phbands.plot_phdispl(qpoint=(0, 0, 0), units="cm^-1", hatches=None, show=False, use_sqrt=True,
                                        normalize=False)
            with self.assertRaises(ValueError):
                # No LO-TO terms
                assert phbands.plot_phdispl(qpoint=1, is_non_analytical_direction=True, show=False)
            assert phbands.plot_phdispl_cartdirs(qpoint=0, units="cm^-1", show=False)
            assert phbands.boxplot(units="ev", mode_range=[2, 4], show=False)

        # Cannot compute PHDOS with q-path
        with self.assertRaises(ValueError):
            phdos = phbands.get_phdos()

        # convert to pymatgen object
        phbands.to_pymatgen()

        # get frozen phonons
        phbands.get_frozen_phonons((0.5, 0.5, 1.0), 1, eta=0.5, max_supercell=[5,5,5])

        assert not phbands.has_linewidths
        phbands.linewidths = np.ones(phbands.shape)
        assert phbands.has_linewidths


class PlotterTest(AbipyTest):

    def test_plot_functions(self):
        """Testing plotting tools for phonons."""
        phbst_filename = abidata.ref_file("trf2_5.out_PHBST.nc")
        with abilab.abiopen(phbst_filename) as nc:
            phbands = nc.phbands

        phb_objects = [
            phbands,
            phbst_filename,
        ]

        phdos_filename = abidata.ref_file("trf2_5.out_PHDOS.nc")
        phdos = PhdosFile(phdos_filename)
        phdos_objects = [
            phdos,
            phdos_filename,
        ]

        if self.has_matplotlib():
            assert phbands_gridplot(phb_objects, titles=["phonons1", "phonons2"],
                                    phdos_objects=phdos_objects, units="meV", show=False)

        phdos.close()


class PhbstFileTest(AbipyTest):

    def test_phbst_file(self):
        """Testing PHBST file."""
        with abilab.abiopen(abidata.ref_file("trf2_5.out_PHBST.nc")) as ncfile:
            assert ncfile.phbands is not None
            for iq, qpt in enumerate(ncfile.qpoints):
                assert hasattr(qpt, "frac_coords")
                assert ncfile.qpoints[ncfile.qindex(qpt)] == qpt
                ii = ncfile.qindex(qpt)
                #print("iq", iq, "qpt", qpt, "ii", ii, "qpoints[ii]", ncfile.qpoints[ii])
                same_ii, same_qpt = ncfile.qindex_qpoint(ii)
                assert same_ii == ii and qpt == same_qpt
                same_ii, same_qpt = ncfile.phbands.qindex_qpoint(qpt)
                assert same_ii == ii and qpt == same_qpt

            qpoint = ncfile.qpoints[0]
            frame = ncfile.get_phframe(qpoint)
            assert frame.qpoint == qpoint

            mode0 = ncfile.get_phmode(qpoint, 0)
            repr(mode0); str(mode0)
            mode0.to_string(with_displ=True)
            assert mode0.qpoint == qpoint

            self.serialize_with_pickle(mode0, test_eq=False)

            last_mode = ncfile.get_phmode(ncfile.qpoints[0], -1)
            assert last_mode > mode0

        # Test notebook
        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_phbst_file_with_loto(self):
        """Testing PHBST file with LOTO terms from anaddb.nc."""
        with abilab.abiopen(abidata.ref_file("ZnSe_hex_886.out_PHBST.nc")) as ncfile:
            phbands = ncfile.phbands

        # Phonon frequencies with non analytical contributions, if calculated, are saved
        # in the anaddb.nc file produced by anaddb. The results should be fetched from there
        # and added to the phonon bands.
        phbands.read_non_anal_from_file(abidata.ref_file("ZnSe_hex_886.anaddb.nc"))
        nana = phbands.non_anal_ph
        assert nana.structure == phbands.structure
        str(nana.structure.reciprocal_lattice)
        self.assert_almost_equal(nana.directions.ravel(),
                [0.1234510847, -0.071274517, 0, 0.1646014463, 0, 0, 0, 0, 0.0751616546])

        for i, cart_direc in enumerate(nana.directions):
            assert nana.has_direction(cart_direc, cartesian=True)

            red_direc = phbands.structure.reciprocal_lattice.get_fractional_coords(cart_direc)
            idir, qdir = phbands.qindex_qpoint(red_direc, is_non_analytical_direction=True)
            assert i == idir
            idir, qdir = phbands.qindex_qpoint(qdir, is_non_analytical_direction=True)
            assert i == idir

        if self.has_matplotlib():
            assert phbands.plot(title="ZnSe with LO-TO splitting", show=False)
            red_direc = phbands.structure.reciprocal_lattice.get_fractional_coords(nana.directions[1])
            assert phbands.plot_phdispl(qpoint=red_direc, is_non_analytical_direction=True, show=False)
            assert phbands.plot_phdispl(qpoint=red_direc, is_non_analytical_direction=True, use_eigvec=True, show=False)

    def test_phbst_robot(self):
        """Testing PHBST robot."""
        paths = abidata.ref_files("trf2_5.out_PHBST.nc", "ZnSe_hex_886.out_PHBST.nc")

        with abilab.PhbstRobot.from_files(paths) as robot:
            assert len(robot) == len(paths)
            repr(robot); str(robot)
            assert robot.to_string(verbose=2)
            df = robot.get_phbands_dataframe()
            for k in ("min_freq", "max_freq", "std_freq"):
                assert k in df

            if self.has_matplotlib():
                assert robot.plot_phdispl(qpoint=(0, 0, 0), show=False)


class PhononBandsPlotterTest(AbipyTest):

    def test_phbands_plotter(self):
        """Testing phbands plotter."""
        phbst_paths = 2 * [abidata.ref_file("trf2_5.out_PHBST.nc")]
        phdos_paths = 2 * [abidata.ref_file("trf2_5.out_PHDOS.nc")]

        plotter = PhononBandsPlotter()
        plotter.add_phbands("AlAs", phbst_paths[0], phdos=phdos_paths[0])
        plotter.add_phbands("Same-AlAs", phbst_paths[1], phdos=phdos_paths[1])
        repr(plotter); str(plotter)

        assert len(plotter.phbands_list) == 2
        assert len(plotter.phdoses_list) == 2
        assert plotter.has_same_formula()

        # __add__ merges two plotters:
        p2 = plotter.add_plotter(plotter)
        assert len(p2.phbands_list) == 2
        assert len(p2.phdoses_list) == 2

        df = dataframe_from_phbands(plotter.phbands_list)
        assert "nqpt" in df

        df = plotter.get_phbands_frame()
        assert df is not None

        if self.has_matplotlib():
            assert plotter.combiplot(units="eV", show=True)
            assert plotter.gridplot(units="meV", show=True)
            assert plotter.boxplot(units="cm-1", show=True)
            assert plotter.combiboxplot(units="Thz", show=True)
            assert plotter.animate(show=False)

        if self.has_nbformat():
            plotter.write_notebook(nbpath=self.get_tmpname(text=True))

        if self.has_ipywidgets():
            assert plotter.ipw_select_plot() is not None

        # Plotter without PHDOS
        plotter = PhononBandsPlotter(key_phbands=[("foo", phbst_paths[0]), ("bar", phbst_paths[1])])

        if self.has_matplotlib():
            assert plotter.combiplot(units="EV", show=True)
            assert plotter.gridplot(units="MEV", show=True)
            assert plotter.boxplot(units="cm^-1", mode_range=[0, 3], swarm=True, show=True)
            assert plotter.combiboxplot(units="THZ", mode_range=[0, 3], show=True)


class PhononDosTest(AbipyTest):

    def test_api(self):
        """Testing PhononDos API with fake data."""
        phdos = PhononDos(mesh=[1, 2, 3], values=[4, 5, 6])
        assert phdos.mesh.tolist() == [1, 2, 3] and phdos.h == 1 and phdos.values.tolist() == [4, 5, 6]
        repr(phdos); str(phdos)
        phdos.idos
        with self.assertRaises(TypeError):
            PhononDos.as_phdos({}, {})

        assert PhononDos.as_phdos(phdos, {}) is phdos
        assert phdos.iw0 == 0

        # From pickle file.
        tmp_path = self.get_tmpname(suffix=".pickle")
        with open(tmp_path, "wb") as fh:
            pickle.dump(phdos, fh)
        same_phdos = PhononDos.as_phdos(tmp_path)
        same_phdos == phdos

        phdos.to_pymatgen()

    def test_from_phdosfile(self):
        """Testing PHDOS from netcdf file."""
        ncfile = PhdosFile(abidata.ref_file("trf2_5.out_PHDOS.nc"))
        repr(ncfile); str(ncfile)
        assert ncfile.to_string(verbose=1)
        assert hasattr(ncfile, "structure")
        nw = len(ncfile.wmesh)
        assert nw == 461
        natom3 = len(ncfile.structure) * 3

        # Read PJDOSes
        assert list(ncfile.pjdos_symbol.keys()) == ["Al", "As"]
        od = ncfile.reader.read_pjdos_symbol_xyz_dict()
        assert list(od.keys()) == ["Al", "As"]
        assert all(v.shape == (3, nw) for v in od.values())

        arr = ncfile.reader.read_pjdos_atdir()
        assert arr.shape == (len(ncfile.structure), 3, nw)

        phdos = ncfile.phdos
        # Test integrals of DOS (the tolerance is a bit low likely due to too coarse meshes)
        self.assert_almost_equal(phdos.integral_value, natom3, decimal=1)

        # Summing projected DOSes over types should give the total DOS.
        pj_sum = sum(pjdos.integral_value for pjdos in ncfile.pjdos_symbol.values())
        self.assert_almost_equal(phdos.integral_value, pj_sum)

        # Summing projected DOSes over types and directions should give the total DOS.
        # phdos_rc_type[ntypat, 3, nomega]
        values = ncfile.reader.read_value("pjdos_rc_type").sum(axis=(0, 1))
        tot_dos = abilab.Function1D(phdos.mesh, values)
        self.assert_almost_equal(phdos.integral_value, tot_dos.integral_value)

        assert phdos == PhononDos.as_phdos(abidata.ref_file("trf2_5.out_PHDOS.nc"))

        # Test Thermodinamics in the Harmonic approximation
        self.assert_almost_equal(phdos.zero_point_energy.to("Ha"), 0.0030872835637731303)

        u = phdos.get_internal_energy()
        self.assert_almost_equal(u.values[0], 0.084009326574073395)
        self.assert_almost_equal(u.values[-1], 0.17270110252071791)

        s = phdos.get_entropy()
        self.assert_almost_equal(s.values[0], 1.6270193052583423e-08)
        self.assert_almost_equal(s.values[-1], 0.00058238779354717)

        cv = phdos.get_cv()
        self.assert_almost_equal(cv.values[0], 5.3715488328604253e-08)
        self.assert_almost_equal(cv.values[-1], 0.00045871909188672578)

        f = phdos.get_free_energy()
        self.assert_almost_equal(f.values, (u - s.mesh * s.values).values)

        self.assertAlmostEqual(phdos.debye_temp, 469.01524830328606)
        self.assertAlmostEqual(phdos.get_acoustic_debye_temp(len(ncfile.structure)), 372.2576492728813)

        assert ncfile.to_pymatgen()

        if self.has_matplotlib():
            assert ncfile.plot_pjdos_type(show=False)
            assert ncfile.plot_pjdos_type(units="cm-1", stacked=False, colormap="viridis", show=False)
            assert ncfile.plot_pjdos_type(units="eV", stacked=True, colormap="jet",
                exchange_xy=True, fontsize=8, show=False)
            assert ncfile.plot_pjdos_cartdirs_type(units="Thz", stacked=True, show=False)
            assert ncfile.plot_pjdos_cartdirs_type(units="meV", stacked=False, alpha=0.5, show=False)
            assert ncfile.plot_pjdos_cartdirs_site(units="meV", stacked=False, alpha=0.5, show=False)
            assert ncfile.plot_pjdos_cartdirs_site(units="meV", view="all", stacked=True, alpha=0.5, show=False)

            assert phdos.plot(units="cm-1", show=False)
            assert phdos.plot_harmonic_thermo(tstar=20, tstop=350, units="eV", formula_units=1, show=False)
            assert phdos.plot_harmonic_thermo(tstar=20, tstop=350, units="Jmol", formula_units=2, show=False)

        # Test notebook
        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

        ncfile.close()


class PhononDosPlotterTest(AbipyTest):

    def test_phdos_plotter(self):
        """Testing PhononDosPlotter."""
        phdos_paths = 2 * [abidata.ref_file("trf2_5.out_PHDOS.nc")]

        plotter = PhononDosPlotter()
        plotter.add_phdos("AlAs", phdos_paths[0])
        plotter.add_phdos("Same-AlAs", phdos_paths[1])
        repr(plotter); str(plotter)
        assert len(plotter.phdos_list) == 2

        if self.has_matplotlib():
            assert plotter.combiplot(show=True)
            assert plotter.gridplot(show=True)
            assert plotter.plot_harmonic_thermo()

        if self.has_ipywidgets():
            assert plotter.ipw_select_plot() is not None
            assert plotter.ipw_harmonic_thermo() is not None

        if self.has_nbformat():
            plotter.write_notebook(nbpath=self.get_tmpname(text=True))


class InteratomicForceConstantsTest(AbipyTest):

    @classmethod
    def setUpClass(cls):
        cls.ddb = DdbFile(os.path.join(test_dir, "AlAs_444_nobecs_DDB"))
        cls.ifc = cls.ddb.anaget_ifc(ifcout=40, ngqpt=[4, 4, 4], verbose=1)

    @classmethod
    def tearDownClass(cls):
        cls.ddb.close()

    def test_filtering(self):
        """Testing IFC filtering."""
        self.ifc.ifc_local_coord_ewald
        self.ifc.ifc_local_coord_short_range

        dist, data = self.ifc.get_ifc_cartesian(atom_indices=[0, 1])
        self.assert_equal(np.shape(data), (80, 3, 3))

        dist, data = self.ifc.get_ifc_local(atom_element="Al")
        self.assert_equal(np.shape(data), (40, 3, 3))

        dist, data = self.ifc.get_ifc_cartesian(min_dist=1, max_dist=10)
        self.assert_equal(np.shape(data), (56, 3, 3))

    def test_plot(self):
        """Testing IFC plot."""
        if not self.has_matplotlib():
            raise unittest.SkipTest("matplotlib missing")

        assert self.ifc.plot_longitudinal_ifc(show=False)
        assert self.ifc.plot_longitudinal_ifc_short_range(show=False)
        assert self.ifc.plot_longitudinal_ifc_ewald(show=False)


class NonAnalyticalPhTest(AbipyTest):

    def test_read_from_file(self):
        """Testing non-analytical terms."""
        # no becs, so no splitting. The test only checks the parsing
        with DdbFile(os.path.join(test_dir, "ZnO_gamma_becs_DDB")) as ddb:
            phbands = ddb.anaget_phmodes_at_qpoint(qpoint=[0, 0, 0], lo_to_splitting=True)

            assert phbands.amu is not None
            #print(phbands.amu)
            # FIXME: I don't like that we Z as key in amu. Should be the symbol
            self.assert_almost_equal(phbands.amu[30.0], 0.6539e+02)
            self.assert_almost_equal(phbands.amu[8.0], 0.159994e+02)
            self.assert_almost_equal(phbands.amu_symbol["Zn"], phbands.amu[30.0])
            self.assert_almost_equal(phbands.amu_symbol["O"], phbands.amu[8.0])

            assert phbands.non_anal_ph is not None
            repr(phbands.non_anal_ph); str(phbands.non_anal_ph)
            assert phbands.structure == phbands.non_anal_ph.structure
            #assert phbands.non_anal_ph.has_direction(direction=,  cartesian=False)

            # TODO should check numerical values (?)
            assert phbands.non_anal_phfreqs is not None
            assert phbands.non_anal_directions is not None
            assert phbands.non_anal_phdispl_cart is not None
            assert phbands.non_anal_dyn_mat_eigenvect is not None
