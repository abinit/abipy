"""Tests for wannier90 module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import os
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest


class TestAbiwanFile(AbipyTest):

    def test_abiwan_nodis(self):
        """Testing abiwan file without DISENTANGLE."""
        #filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "foo_ABIWAN.nc")
        with abilab.abiopen(filepath) as abiwan:
            repr(abiwan); str(abiwan)
            assert abiwan.to_string(verbose=2)
            #assert abiwan.version == "2.1.0+git"
            #assert abiwan.structure.formula == "Ga1 As1"
            #self.assert_almost_equal(abiwan.structure.frac_coords.ravel(),
            #        [0.00000, 0.00000, 0.00000, 0.25000, 0.25000, 0.25000])
            #assert not abiwan.warnings
            #for k in ("MAIN", "WANNIERISE"):
            #    assert k in abiwan.params_section
            #assert not abiwan.use_disentangle
            #assert abiwan.nwan == 4
            #assert np.all(abiwan.grid_size == 2)
            #assert len(abiwan.conv_df) == 20 + 1
            #assert abiwan.conv_df.O_D[0] == 0.0083198 and abiwan.conv_df.O_OD[0] == 0.5036294
            #assert abiwan.conv_df.O_D[20] == 0.0080300 and abiwan.conv_df.O_OD[20] == 0.5019880

            ## numpy array (nwan, nstep, ...)
            ##natom = len(abiwan.structure)
            #nstep = 21
            #assert abiwan.wf_centers.shape == (abiwan.nwan, nstep, 3)
            #self.assert_equal(abiwan.wf_centers[1, 20], [-0.866253,  0.866253,  0.866253])
            #assert abiwan.wf_spreads.shape == (abiwan.nwan, nstep)
            #self.assert_equal(abiwan.wf_spreads[1, 20], 1.11672024)

            #if self.has_matplotlib():
            #    assert abiwan.plot(show=False)
            #    assert abiwan.plot_centers_spread(show=False)

            #if self.has_nbformat():
            #    assert abiwan.write_notebook(nbpath=self.get_tmpname(text=True))

    #def test_example03_silicon(self):
    #    """Parsing example02_silicon.wout with DISENTANGLE"""
    #    with abilab.abiopen(filepath) as abiwan:
    #        repr(abiwan); str(abiwan)
    #        assert abiwan.to_string(verbose=2)
    #        assert abiwan.version == "2.1.0+git"
    #        assert abiwan.structure.formula == "Si2"
    #        self.assert_almost_equal(abiwan.structure.frac_coords.ravel(),
    #                [-0.25000, 0.75000, -0.25000, 0.00000, 0.00000, 0.00000])
    #        assert not abiwan.warnings
    #        for k in ("MAIN", "WANNIERISE", "DISENTANGLE"):
    #            assert k in abiwan.params_section
    #        assert abiwan.use_disentangle
    #        assert abiwan.nwan == 8
    #        assert np.all(abiwan.grid_size == 4)
    #        assert len(abiwan.conv_df) == 6 + 1
    #        assert abiwan.conv_df.O_D[1] == 0.1213986 and abiwan.conv_df.O_OD[1] == 2.7017701
    #        assert abiwan.conv_df.O_D.values[-1] == 0.1054702 and abiwan.conv_df.O_OD.values[-1] == 2.5449106

    #        # numpy array (nwan, nstep, ...)
    #        nstep = len(abiwan.conv_df.O_D)
    #        assert abiwan.wf_centers.shape == (abiwan.nwan, nstep, 3)
    #        self.assert_equal(abiwan.wf_centers[7, -1], [0.888643,  0.888652,  1.810090 ])
    #        assert abiwan.wf_spreads.shape == (abiwan.nwan, nstep)
    #        self.assert_equal(abiwan.wf_spreads[7, -1], 1.81245236)

    #        if self.has_matplotlib():
    #            assert abiwan.plot(show=False)
    #            assert abiwan.plot_centers_spread(show=False)

    #        if self.has_nbformat():
    #            assert abiwan.write_notebook(nbpath=self.get_tmpname(text=True))
