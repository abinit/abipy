"""Tests for phonons"""
from __future__ import print_function, division

import unittest
import os
import numpy as np
import abipy.data as abidata

from abipy.dfpt.phonons import PhononBands, PhononDos, PhdosFile, InteratomicForceConstants
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

        # Test xmgrace
        phbands.to_xmgrace(self.get_tmpname(text=True))

        # Test convertion to eigenvectors. Verify that they are orthonormal
        # Allow relatively large tolerance due to possible mismatching in the atomic masses between abinit and pmg
        eig = phbands.dyn_mat_eigenvect
        self.assertTrue(np.allclose(np.dot(eig[0], eig[0].T), np.eye(len(eig[0]), dtype=np.complex), atol=1e-5, rtol=1e-3))

        #dos = phbands.get_phdos()
        #print(dos)


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

        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))


class InteratomicForceConstantsTest(AbipyTest):

    @classmethod
    def setUpClass(cls):
        cls.ddb = DdbFile(os.path.join(test_dir, "AlAs_444_nobecs_DDB"))

        cls.ifc = cls.ddb.anaget_ifc(ifcout=40, ngqpt=[4,4,4], verbose=1)

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
