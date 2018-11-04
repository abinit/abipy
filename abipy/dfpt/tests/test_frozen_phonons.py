"""Tests for frozen_phonons"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata

from abipy.dfpt.frozen_phonons import FrozenPhonon
from abipy.dfpt.phonons import PhononBands
from abipy.core.testing import AbipyTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


class FrozenPhononTest(AbipyTest):

    def test_frozen(self):
        """Base tests for FrozenPhonon"""
        filename = abidata.ref_file("trf2_5.out_PHBST.nc")
        phbands = PhononBands.from_file(filename)

        qpt_frac_coords = [0.5, 0.5, 0.5]
        fp = FrozenPhonon.from_phbands(phbands, qpt_frac_coords, 0 ,
                                       etas=[-0.2, -0.1, 0, 0.1, 0.2], max_supercell=[5, 5, 5])

        self.assertArrayEqual(fp.scale_matrix, [[-1,  0,  1], [-1,  1,  0], [-1, -1,  0]])

        w = phbands.phfreqs[phbands.qindex(qpt_frac_coords), 0]

        energies = [0.0704, 0.0176, 0. , 0.0175, 0.0703]

        with self.assertRaises(ValueError):
            fp.energies = energies[:3]

        fp.energies = energies

        assert fp.ieta0 == 2
        assert fp.n_displ == 5

        fit_data = fp.fit_to_frequency()
        self.assertAlmostEqual(w, fit_data.freq, places=5)

        fit_data = fp.fit_to_frequency(min_fit_eta=-0.15, max_fit_eta=0.15)
        self.assertNotAlmostEqual(w, fit_data.freq, places=5)

        self.assertArrayAlmostEqual(fp.qpt_cart_coords, [0.55954285702497808, 0.55954285702497808, 0.55954285702497808])

        if self.has_matplotlib():
            assert fp.plot_fit_energies(freq=w, show=False)
            assert fp.plot_anharmonic_contribution(freq=w, show=False)
            assert fp.plot_anharmonic_contribution(freq=w, relative=True, show=False)