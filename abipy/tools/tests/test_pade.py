"""Tests for frozen_phonons"""
import os
import numpy as np
#import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.tools.pade import SigmaPade


class PadeTest(AbipyTest):

    def test_sigma_pade(self):
        """Testing SigmaPade interface."""

        # List of points and values for the Pade'
        zs = np.array([0.3204639033968435j, 1.159257347625358j, 2.7851769489568112j,
                       6.6074592214004655j, 16.793387237006716j, 50.739714260170196j],
                       dtype=complex)

        f_zs = np.array([(-3.4591316381250774-0.09890936491537991j), (-3.4538645067311458-0.24793469903457657j),
                        (-3.4173978321863663-0.6621397224074206j), (-3.257949585597619-1.324699039193591j),
                        (-2.523177170083569-2.4583423787866163j), (-0.8792031880464662-2.2566436943461583j)],
                        dtype=complex)

        pade = SigmaPade(zs, f_zs)

        # Evaluate Pade' as w_vals.
        w_vals = np.array([-0.49705069294512505, -0.397050692945125, -0.297050692945125, -0.197050692945125])
        this_f, this_fp = pade.eval(w_vals)

        ref_f = [-3.36245542+0.0371258j, -3.38196076+0.03677179j, -3.40149899+0.036502j, -3.42108152+0.03632633j]
        self.assert_almost_equal(this_f, ref_f, decimal=7)
        ref_fp = [-0.19492091-0.00392504j, -0.19520083-0.00313668j, -0.19558294-0.00224237j, -0.1960902-0.00125779j]
        self.assert_almost_equal(this_fp, ref_fp, decimal=7)
