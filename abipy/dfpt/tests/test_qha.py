"""Tests for frozen_phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import warnings
import abipy.data as abidata

from abipy.dfpt.qha import QHA, QHA3PF, QHA3P
from abipy.dfpt.phonons import PhononBands
from abipy.core.testing import AbipyTest


class QhaTest(AbipyTest):

    @classmethod
    def setUpClass(cls):
        cls.strains = [-4, -2, 0, 2, 4, 6]
        path = os.path.join(abidata.dirpath, "refs", "si_qha")
        cls.gsr_paths = [os.path.join(path, "mp-149_{:+d}_GSR.nc".format(s)) for s in cls.strains]
        cls.dos_paths = [os.path.join(path, "mp-149_{:+d}_PHDOS.nc".format(s)) for s in cls.strains]
        cls.phbs = [PhononBands.from_file(os.path.join(path, "mp-149_{:+d}_PHBST.nc".format(s))) for s in
                    cls.strains[2:4]]

    def test_qha(self):
        """Base tests for QHA"""

        qha = QHA.from_files(self.gsr_paths, self.dos_paths)

        self.assertEqual(qha.nvols, len(self.strains))
        self.assertEqual(qha.natoms, 2)

        f = qha.fit_energies(tstart=0, tstop=300, num=3)
        self.assertArrayEqual(f.tot_en.shape, (len(self.strains), 3))
        self.assertAlmostEqual(f.min_en[0], -230.15471148501612, places=6)

        self.assertEqual(qha.eos._eos_name, "vinet")
        qha.set_eos("murnaghan")
        self.assertEqual(qha.eos._eos_name, "murnaghan")

        te = qha.get_thermal_expansion_coeff(num=4)
        self.assertAlmostEqual(te.values[1], 1.4676820862386381e-05)

        self.assertAlmostEqual(qha.get_vol_at_t(200), 41.07441539803265, places=4)
        temps = qha.get_t_for_vols([40, 41.065, 41.1])
        self.assertEqual(len(temps[0]), 0)
        self.assertEqual(len(temps[1]), 2)
        self.assertEqual(len(temps[2]), 1)

        tmpdir = self.mkdtemp()
        qha.write_phonopy_qha_inputs(num=3, path=tmpdir)

        if self.has_matplotlib():
            assert qha.plot_energies(num=6, show=False)
            assert qha.plot_thermal_expansion_coeff(num=6, show=False)
            assert qha.plot_vol_vs_t(num=6, show=False)
            # fake temperatures to test the plotting function.
            assert qha.plot_phbs(self.phbs, temperatures=[10, 20], show=False)

    def test_phonopy_object(self):
        self.skip_if_not_phonopy()

        qha = QHA.from_files(self.gsr_paths, self.dos_paths)

        from phonopy.qha import QHA as QHA_phonopy
        qha_ph = qha.get_phonopy_qha(tstop=500, num=11)
        self.assertIsInstance(qha_ph, QHA_phonopy)
        qha_ph.run()
        if self.has_matplotlib():
            qha_ph.plot_thermal_expansion()


class Qha3pfTest(AbipyTest):

    @classmethod
    def setUpClass(cls):
        cls.strains = [-4, -2, 0, 2, 4, 6]
        path = os.path.join(abidata.dirpath, "refs", "si_qha")
        cls.gsr_paths = [os.path.join(path, "mp-149_{:+d}_GSR.nc".format(s)) for s in cls.strains]
        cls.dos_paths = [os.path.join(path, "mp-149_{:+d}_PHDOS.nc".format(s)) for s in cls.strains]

    def test_qha(self):
        """Base tests for QHA3PF"""

        qha = QHA3PF.from_files(self.gsr_paths, self.dos_paths[1:4], ind_doses=[1, 2, 3])

        self.assertEqual(qha.nvols, len(self.strains))
        self.assertEqual(qha.natoms, 2)

        f = qha.fit_energies(tstart=0, tstop=300, num=3)
        self.assertArrayEqual(f.tot_en.shape, (len(self.strains), 3))
        self.assertAlmostEqual(f.min_en[0], -230.15433223031306, places=6)

        self.assertEqual(qha.eos._eos_name, "vinet")
        qha.set_eos("murnaghan")
        self.assertEqual(qha.eos._eos_name, "murnaghan")

        te = qha.get_thermal_expansion_coeff(num=4)
        self.assertAlmostEqual(te.values[1], 1.2773693323408941e-05)

        self.assertAlmostEqual(qha.get_vol_at_t(200), 41.10212044734946, places=4)

        tmpdir = self.mkdtemp()
        qha.write_phonopy_qha_inputs(num=3, path=tmpdir)

        if self.has_matplotlib():
            assert qha.plot_energies(num=6, show=False)
            assert qha.plot_thermal_expansion_coeff(num=6, show=False)
            assert qha.plot_vol_vs_t(num=6, show=False)

    def test_phonopy_object(self):
        self.skip_if_not_phonopy()

        qha = QHA3PF.from_files(self.gsr_paths, self.dos_paths[1:4], ind_doses=[1, 2, 3])

        from phonopy.qha import QHA as QHA_phonopy
        qha_ph = qha.get_phonopy_qha(tstop=500, num=11)
        self.assertIsInstance(qha_ph, QHA_phonopy)
        qha_ph.run()
        if self.has_matplotlib():
            qha_ph.plot_thermal_expansion()


class Qha3pTest(AbipyTest):

    @classmethod
    def setUpClass(cls):
        cls.strains = [-4, -2, 0, 2, 4, 6]
        path = os.path.join(abidata.dirpath, "refs", "si_qha")
        cls.gsr_paths = [os.path.join(path, "mp-149_{:+d}_GSR.nc".format(s)) for s in cls.strains]
        cls.gruns_path = os.path.join(path, "mp-149_GRUNS.nc")
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_qha(self):
        """Base tests for QHA3PF"""

        qha = QHA3P.from_files(self.gsr_paths, self.gruns_path, ind_doses=[1, 2, 3])

        self.assertEqual(qha.nvols, len(self.strains))
        self.assertEqual(qha.natoms, 2)

        f = qha.fit_energies(tstart=0, tstop=300, num=3)
        self.assertArrayEqual(f.tot_en.shape, (len(self.strains), 3))
        self.assertAlmostEqual(f.min_en[0], -230.15458036894498, places=6)

        self.assertEqual(qha.eos._eos_name, "vinet")
        qha.set_eos("murnaghan")
        self.assertEqual(qha.eos._eos_name, "murnaghan")

        te = qha.get_thermal_expansion_coeff(num=4)
        self.assertAlmostEqual(te.values[1], 1.2725767394824783e-05)

        self.assertAlmostEqual(qha.get_vol_at_t(200), 41.1083743159003, places=4)

        tmpdir = self.mkdtemp()
        qha.write_phonopy_qha_inputs(num=3, path=tmpdir)

        if self.has_matplotlib():
            assert qha.plot_energies(num=6, show=False)
            assert qha.plot_thermal_expansion_coeff(num=6, show=False)
            assert qha.plot_vol_vs_t(num=6, show=False)

    def test_phonopy_object(self):
        self.skip_if_not_phonopy()

        qha = QHA3P.from_files(self.gsr_paths, self.gruns_path, ind_doses=[1, 2, 3])

        from phonopy.qha import QHA as QHA_phonopy
        qha_ph = qha.get_phonopy_qha(tstop=500, num=11)
        self.assertIsInstance(qha_ph, QHA_phonopy)
        qha_ph.run()
        if self.has_matplotlib():
            qha_ph.plot_thermal_expansion()