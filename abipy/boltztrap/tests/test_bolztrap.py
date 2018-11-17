"""Tests for boltztrap module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import collections
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.boltztrap import AbipyBoltztrap, BoltztrapResult
from abipy import abilab


class AbipyBoltztrapTest(AbipyTest):

    # TODO: Need new files with IBZ.
    def test_sigeph_boltztrap(self):
        """Test boltztrap interpolation"""
        self.skip_if_not_bolztrap2()

        with abilab.abiopen(abidata.ref_file("diamond_444q_full_SIGEPH.nc")) as sigeph:
            bt = AbipyBoltztrap.from_sigeph(sigeph)
            repr(bt); str(bt)
            assert bt.to_string(verbose=2)

        # get equivalences
        assert bt.rmesh == (17, 17, 17)
        assert bt.nequivalences == 67

        # get coefficients
        assert bt.ncoefficients == 53
        bt.dump_rsphere(self.get_tmpname(text=True))

        # get ebands using boltztrap
        bt_ebands = bt.get_ebands()

        # Get boltztrap results using different DOS methods
        btr = bt.run(dos_method="histogram")
        btr = bt.run(npts=500,dos_method="gaussian:0.5 eV")
        btr = bt.run(npts=500,dos_method="lorentzian:0.5 eV")
        repr(btr); str(btr)
        assert btr.to_string(verbose=2)

        # Test pickle
        pickle_file = self.get_tmpname(suffix="diamond.npy")
        btr.pickle(pickle_file)
        same_result = BoltztrapResult.from_pickle(pickle_file)
        self.assert_equal(btr.tmesh, same_result.tmesh)

        if self.has_matplotlib():
            # Plot the density of states and VVDOS for multiple temperatures
            assert btr.plot_dos_vvdos(show=False)

            # Plot transport related quantities for different combinations of
            # tau temperature and boltztrap temperature
            assert btr.plot('sigma', itemp_list=None, itau_list=[3], show=False)
            assert btr.plot('seebeck', itemp_list=[3], itau_list=[1,2], show=False)
            assert btr.plot('powerfactor', itemp_list=[3], itau_list=None, show=False)
            assert btr.plot_transport(show=False)
