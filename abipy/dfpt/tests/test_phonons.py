"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import unittest
import os
import pickle
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.dfpt.phonons import (PhononBands, PhononDos, PhdosFile, InteratomicForceConstants, phbands_gridplot,
        PhononBandsPlotter, PhononDosPlotter, frame_from_phbands)
from abipy.dfpt.ddb import DdbFile
from abipy.core.testing import AbipyTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


class TestUnitTools(AbipyTest):

    def test_units_api(self):
        from abipy.dfpt.phonons import _factor_ev2units, _unit_tag, _dos_label_from_units
        for units in ["ev", "meV" ,"ha", "cm-1", "cm^-1", "Thz"]:
            _factor_ev2units(units)
            _unit_tag(units)
            _dos_label_from_units(units)

        for func in [_factor_ev2units, _unit_tag, _dos_label_from_units]:
            with self.assertRaises(KeyError):
                func("foo")


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
            same_phbands = PhononBands.as_phbands(nc)
            self.assert_equal(same_phbands.phfreqs, phbands.phfreqs)

        self.serialize_with_pickle(phbands, protocols=[-1], test_eq=False)

        # From pickle file.
        tmp_path = self.get_tmpname(suffix=".pickle")
        with open(tmp_path, "wb") as fh:
            pickle.dump(phbands, fh)
        same_phbands = PhononBands.as_phbands(tmp_path)
        self.assert_equal(same_phbands.phfreqs, phbands.phfreqs)

        assert phbands.minfreq == 0.0
        #self.assertEqual(phbands.maxfreq, 30)

        # Test XYZ vib
        phbands.create_xyz_vib(iqpt=0, filename=self.get_tmpname(text=True), max_supercell=[4,4,4])
        # Test ascii file
        phbands.create_ascii_vib(iqpts=0, filename=self.get_tmpname(text=True), pre_factor=1)
        # Test phononwebsite file
        phbands.create_phononwebsite_json(filename=self.get_tmpname(text=True), name='test')
        # Test xmgrace
        phbands.to_xmgrace(self.get_tmpname(text=True))

        df = phbands.to_dataframe()
        assert "freq" in df and "mode" in df
        self.assert_almost_equal(df["freq"].values.min(), 0)

        umodes = phbands.get_unstable_modes(below_mev=-1000)
        assert len(umodes) == 0

        # Test convertion to eigenvectors. Verify that they are orthonormal
        # Allow relatively large tolerance due to possible mismatching in the atomic masses between abinit and pmg
        eig = phbands.dyn_mat_eigenvect
        assert np.allclose(np.dot(eig[0], eig[0].T), np.eye(len(eig[0]), dtype=np.complex), atol=1e-5, rtol=1e-3)

        # Mapping reduced coordinates -> labels
        qlabels = {
            (0,0,0): "$\Gamma$",
            (0.375, 0.375, 0.75): "K",
            (0.5, 0.5, 1.0): "X",
            (0.5, 0.5, 0.5): "L",
            (0.5, 0.0, 0.5): "X",
            (0.5, 0.25, 0.75): "W",
        }

        if self.has_matplotlib():
            assert phbands.plot(units="Thz", show=False)
            assert phbands.plot_fatbands(units="ha", qlabels=qlabels, show=False)
            assert phbands.plot_fatbands(phdos_file=abidata.ref_file("trf2_5.out_PHDOS.nc"), units="thz", show=False)
            assert phbands.plot_colored_matched(units="cm^-1", show=False)
            assert phbands.boxplot(units="ev", mode_range=[2, 4], show=False)

        # Cannot compute PHDOS with q-path
        with self.assertRaises(ValueError):
            phdos = phbands.get_phdos()

        # convert to pymatgen object
        phbands.to_pymatgen()


class PlotterTest(AbipyTest):

    def test_plot_functions(self):
        """Testing plotting tools for phonons."""
        phbst_filename = abidata.ref_file("trf2_5.out_PHBST.nc")
        from abipy import abilab
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
        from abipy import abilab
        with abilab.abiopen(abidata.ref_file("trf2_5.out_PHBST.nc")) as ncfile:
            assert ncfile.phbands is not None
            for iq, qpt in enumerate(ncfile.qpoints):
                assert ncfile.qpoints[ncfile.qindex(qpt)] == qpt
                ii = ncfile.qindex(qpt)
                #print("iq", iq, "qpt", qpt, "ii", ii, "qpoints[ii]", ncfile.qpoints[ii])
                #assert ii == iq

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

        df = frame_from_phbands(plotter.phbands_list)
        assert "nqpt" in df

        df = plotter.get_phbands_frame()
        assert df is not None

        if self.has_matplotlib():
            assert plotter.combiplot(units="eV", show=True)
            assert plotter.gridplot(units="meV", show=True)
            assert plotter.boxplot(units="cm-1", show=True)
            assert plotter.combiboxplot(units="Thz", show=True)

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

    def test_from_phdosfile(self):
        """Testing PHDOS from netcdf file."""
        ncfile = PhdosFile(abidata.ref_file("trf2_5.out_PHDOS.nc"))
        repr(ncfile); str(ncfile)
        assert hasattr(ncfile, "structure")
        nw = len(ncfile.wmesh)
        assert nw == 461

        # Read PJDOSes
        assert list(ncfile.pjdos_symbol.keys()) == ["Al", "As"]
        od = ncfile.reader.read_pjdos_symbol_rc_dict()
        assert list(od.keys()) == ["Al", "As"]
        assert all(v.shape == (3, nw) for v in od.values())

        phdos = ncfile.phdos
        assert phdos == PhononDos.as_phdos(abidata.ref_file("trf2_5.out_PHDOS.nc"))

        # Test Thermodinamics in the Harmonic approximation
        self.assert_almost_equal(phdos.zero_point_energy.to("Ha"), 0.0030872835637731303)

        u = phdos.get_internal_energy()
        self.assert_almost_equal(u.values[0], 0.16801859138305908)
        self.assert_almost_equal(u.values[-1], 0.25671036732970359)

        s = phdos.get_entropy()
        self.assert_almost_equal(s.values[0], 1.6270193052583423e-08)
        self.assert_almost_equal(s.values[-1], 0.00058238779354717)

        cv = phdos.get_cv()
        self.assert_almost_equal(cv.values[0], 5.3715488328604253e-08)
        self.assert_almost_equal(cv.values[-1], 0.00045871909188672578)

        f = phdos.get_free_energy()
        self.assert_almost_equal(f.values, (u - s.mesh * s.values).values)

        if self.has_matplotlib():
            assert ncfile.plot_pjdos_type(show=False)
            assert ncfile.plot_pjdos_type(units="cm-1", stacked=False, colormap="viridis", show=False)
            assert ncfile.plot_pjdos_redirs_type(units="Thz", stacked=True, show=False)
            assert ncfile.plot_pjdos_redirs_type(units="meV", stacked=False, alpha=0.5, show=False)
            assert ncfile.plot_pjdos_redirs_site(units="meV", stacked=False, alpha=0.5, show=False)
            assert ncfile.plot_pjdos_redirs_site(units="meV", stacked=True, alpha=0.5, show=False)

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

            assert phbands.non_anal_ph is not None
            str(phbands.non_anal_ph)
            assert phbands.structure == phbands.non_anal_ph.structure
            #assert phbands.non_anal_ph.has_direction(direction= cartesian=False)

            # TODO should check numerical values (?)
            assert phbands.non_anal_phfreqs is not None
            assert phbands.non_anal_directions is not None
            assert phbands.non_anal_phdispl_cart is not None
            assert phbands.non_anal_dyn_mat_eigenvect is not None
