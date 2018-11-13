"""Tests for phtk module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.dfpt.vsound import SoundVelocity


class SoundVelocityTest(AbipyTest):

    def test_from_ddb(self):

        ddb_path = os.path.join(abidata.dirpath, "refs", "si_sound_vel", "Si_DDB")
        sv = SoundVelocity.from_ddb(ddb_path, num_points=5)

        self.assertAlmostEqual(sv.sound_velocities[0,0], 5573.880942238931, places=4)
        self.assertEqual(sv.labels[0], "X")

        df = sv.get_dataframe()

        if self.has_matplotlib():
            assert sv.plot_fit_freqs_dir(0, show=False)
            assert sv.plot(show=False)

        if self.has_nbformat():
            assert sv.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_from_phb(self):
        phb_path = os.path.join(abidata.dirpath, "refs", "si_sound_vel", "Si_sound_PHBST.nc")
        sv = SoundVelocity.from_phbst(phb_path, ignore_neg_freqs=False)

        self.assertAlmostEqual(sv.sound_velocities[1, 1], 5624.764713502484, places=4)
        self.assertIsNone(sv.labels)
        df = sv.get_dataframe()
