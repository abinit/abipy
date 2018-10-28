"""Tests for boltztrap module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import collections
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.boltztrap.boltztrap import AbipyBoltztrap, BoltztrapResults
from abipy import abilab


class AbipyBoltztrapTest(AbipyTest):

    # TODO: Need new files with IBZ.
    def test_sigeph_boltztrap(self):
        """Test boltztrap interpolation"""
        self.skip_if_not_bolztrap2()

        sigeph = abilab.abiopen(abidata.ref_file("diamond_444q_full_SIGEPH.nc"))

        bt = AbipyBoltztrap.from_sigeph(sigeph)

        # get equivalences
        assert bt.rmesh == (5,5,5)
        assert bt.nequivalences == 5

        # get coefficients
        assert bt.ncoefficients == 53

        #get results
        btr = bt.run()

        #boltztrap_results
        btr.plot_dos()

        btr.pickle('diamond.npy')
        #fig = plt.figure
        #ax = fig.add_subplot(1,1,1)
        #btr.plot_sigma(ax)
        #btr.plot_seebeck(ax)
