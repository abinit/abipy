"""Tests for wannier90 module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import os
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest


class TestWoutFile(AbipyTest):

    def test_example01_gaas(self):
        """Parsing example01_gaas.wout"""
        filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "example01_gaas.wout")
        with abilab.abiopen(filepath) as wout:
            repr(wout); str(wout)
            assert wout.to_string(verbose=2)
            assert wout.version == "2.1.0+git"
            assert wout.structure.formula == "Ga1 As1"
            self.assert_almost_equal(wout.structure.frac_coords.ravel(),
                    [0.00000, 0.00000, 0.00000, 0.25000, 0.25000, 0.25000])
            assert not wout.warnings
            for k in ("MAIN", "WANNIERISE"):
                assert k in wout.params_section
            assert not wout.use_disentangle
            assert wout.nwan == 4
            assert np.all(wout.grid_size == 2)
            assert len(wout.conv_df) == 20 + 1
            assert wout.conv_df.O_D[0] == 0.0083198 and wout.conv_df.O_OD[0] == 0.5036294
            assert wout.conv_df.O_D[20] == 0.0080300 and wout.conv_df.O_OD[20] == 0.5019880

            # numpy array (nwan, nstep, ...)
            #natom = len(wout.structure)
            nstep = 21
            assert wout.wf_centers.shape == (wout.nwan, nstep, 3)
            self.assert_equal(wout.wf_centers[1, 20], [-0.866253,  0.866253,  0.866253])
            assert wout.wf_spreads.shape == (wout.nwan, nstep)
            self.assert_equal(wout.wf_spreads[1, 20], 1.11672024)

            if self.has_matplotlib():
                assert wout.plot(show=False)
                assert wout.plot_centers_spread(show=False)

            if self.has_nbformat():
                assert wout.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_example03_silicon(self):
        """Parsing example02_silicon.wout with DISENTANGLE"""
        filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "example03_silicon.wout")
        with abilab.abiopen(filepath) as wout:
            repr(wout); str(wout)
            assert wout.to_string(verbose=2)
            assert wout.version == "2.1.0+git"
            assert wout.structure.formula == "Si2"
            self.assert_almost_equal(wout.structure.frac_coords.ravel(),
                    [-0.25000, 0.75000, -0.25000, 0.00000, 0.00000, 0.00000])
            assert not wout.warnings
            for k in ("MAIN", "WANNIERISE", "DISENTANGLE"):
                assert k in wout.params_section
            assert wout.use_disentangle
            assert wout.nwan == 8
            assert np.all(wout.grid_size == 4)
            assert len(wout.conv_df) == 6 + 1
            assert wout.conv_df.O_D[1] == 0.1213986 and wout.conv_df.O_OD[1] == 2.7017701
            assert wout.conv_df.O_D.values[-1] == 0.1054702 and wout.conv_df.O_OD.values[-1] == 2.5449106

            # numpy array (nwan, nstep, ...)
            nstep = len(wout.conv_df.O_D)
            assert wout.wf_centers.shape == (wout.nwan, nstep, 3)
            self.assert_equal(wout.wf_centers[7, -1], [0.888643,  0.888652,  1.810090 ])
            assert wout.wf_spreads.shape == (wout.nwan, nstep)
            self.assert_equal(wout.wf_spreads[7, -1], 1.81245236)

            if self.has_matplotlib():
                assert wout.plot(show=False)
                assert wout.plot_centers_spread(show=False)

            if self.has_nbformat():
                assert wout.write_notebook(nbpath=self.get_tmpname(text=True))
