"""Tests for phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import unittest
import os
import numpy as np
import abipy.data as abidata

from abipy.dfpt.phonons import PhononBands, PhononDos, PhdosFile, InteratomicForceConstants, phbands_gridplot
from abipy.dfpt.ddb import DdbFile
from abipy.core.testing import *

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


class PhononBandsTest(AbipyTest):

    def test_base(self):
        """Base tests for PhononBands"""
        filename = abidata.ref_file("trf2_5.out_PHBST.nc")
        phbands = PhononBands.from_file(filename)
        print(phbands)
        assert PhononBands.as_phbands(phbands) is phbands
        assert np.array_equal(PhononBands.as_phbands(filename).phfreqs, phbands.phfreqs)

        self.serialize_with_pickle(phbands, protocols=[-1], test_eq=False)

        self.assertEqual(phbands.minfreq, 0.0)
        #self.assertEqual(phbands.maxfreq, 30)

        # Test XYZ vib
        phbands.create_xyz_vib(iqpt=0, filename=self.get_tmpname(text=True), max_supercell=[4,4,4])

        # Test ascii file
        phbands.create_ascii_vib(iqpts=0, filename=self.get_tmpname(text=True), pre_factor=1)

        # Test phononwebsite file
        phbands.create_phononwebsite_json(filename=self.get_tmpname(text=True), name='test')

        # Test xmgrace
        phbands.to_xmgrace(self.get_tmpname(text=True))

        # Test convertion to eigenvectors. Verify that they are orthonormal
        # Allow relatively large tolerance due to possible mismatching in the atomic masses between abinit and pmg
        eig = phbands.dyn_mat_eigenvect
        assert np.allclose(np.dot(eig[0], eig[0].T), np.eye(len(eig[0]), dtype=np.complex), atol=1e-5, rtol=1e-3)

        if self.has_matplotlib():
            phbands.plot(show=False)
            phbands.plot_fatbands(show=False)
            phbands.plot_colored_matched(show=False)

        # Cannot compute PHDOS with q-path
        with self.assertRaises(ValueError):
            phdos = phbands.get_phdos()

        # convert to pymatgen object
        phbands.to_pymatgen()


class PlotterTest(AbipyTest):

    def test_plot_functions(self):
        """Testing plotting tools for phonons."""
        if not self.has_matplotlib():
            raise unittest.SkipTest("matplotlib missing")

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

        fig = phbands_gridplot(phb_objects, titles=["phonons1", "phonons2"],
                               phdos_objects=phdos_objects, show=False)
        assert fig is not None
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
            print(mode0)
            mode0.to_string(with_displ=True)
            assert mode0.qpoint == qpoint

            self.serialize_with_pickle(mode0, test_eq=False)

            last_mode = ncfile.get_phmode(ncfile.qpoints[0], -1)
            assert last_mode > mode0

        # Test notebook
        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))


class PhononDosTest(AbipyTest):

    def test_api(self):
        """Testing PhononDos API with fake data."""
        dos = PhononDos(mesh=[1, 2, 3], values=[4, 5, 6])
        assert dos.mesh.tolist() == [1,2,3] and dos.h == 1 and dos.values.tolist() == [4,5,6]
        print(dos)
        dos.idos
        assert PhononDos.as_phdos(dos, {}) is dos
        assert dos.iw0 == 0

    def test_from_phdosfile(self):
        ncfile = PhdosFile(abidata.ref_file("trf2_5.out_PHDOS.nc"))
        assert hasattr(ncfile, "structure")

        phdos = ncfile.phdos

        # Thermodinamics in the Harmonic approximation
        #self.assert_almost_equal(phdos.zero_point_energy, )
        f = phdos.get_free_energy()
        #self.assert_almost_equal(f.values, )
        u = phdos.get_internal_energy()
        #self.assert_almost_equal(u.values, )
        s = phdos.get_entropy()
        #self.assert_almost_equal(s.values, )
        cv = phdos.get_cv()
        #self.assert_almost_equal(cv.values, )

        if self.has_matplotlib():
            ncfile.plot_pjdos_type(show=False)
            phdos.plot(show=False)
            phdos.plot_harmonic_thermo(tstar=20, tstop=350)

        # Test notebook
        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

        ncfile.close()


class InteratomicForceConstantsTest(AbipyTest):

    @classmethod
    def setUpClass(cls):
        cls.ddb = DdbFile(os.path.join(test_dir, "AlAs_444_nobecs_DDB"))
        cls.ifc = cls.ddb.anaget_ifc(ifcout=40, ngqpt=[4,4,4], verbose=1)

    @classmethod
    def tearDownClass(cls):
        cls.ddb.close()

    def test_filtering(self):

        self.ifc.ifc_local_coord_ewald
        self.ifc.ifc_local_coord_short_range

        dist, data = self.ifc.get_ifc_cartesian(atom_indices=[0,1])
        self.assert_equal(np.shape(data), (80, 3, 3))

        dist, data = self.ifc.get_ifc_local(atom_element="Al")
        self.assert_equal(np.shape(data), (40, 3, 3))

        dist, data = self.ifc.get_ifc_cartesian(min_dist=1, max_dist=10)
        self.assert_equal(np.shape(data), (56, 3, 3))

    def test_plot(self):
        if not self.has_matplotlib():
            raise unittest.SkipTest("matplotlib missing")

        self.ifc.plot_longitudinal_ifc(show=False)
        self.ifc.plot_longitudinal_ifc_short_range(show=False)
        self.ifc.plot_longitudinal_ifc_ewald(show=False)


class NonAnalyticalPhTest(AbipyTest):

    def test_read_from_file(self):
        # no becs, so no splitting. The test only checks the parsing
        ddb = DdbFile(os.path.join(test_dir, "ZnO_gamma_becs_DDB"))

        phbands = ddb.anaget_phmodes_at_qpoint(qpoint=[0, 0, 0], lo_to_splitting=True)

        phbands.non_anal_phfreqs
        ddb.close()
