"""Tests for boltztrap module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import collections
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.boltztrap.boltztrap import AbipyBoltztrap
from abipy import abilab


class AbipyBoltztrapTest(AbipyTest):

    # TODO: Need new files with IBZ.
    def test_sigeph_boltztrap(self):
        """Test boltztrap interpolation"""
        sigeph = abilab.abiopen(abidata.ref_file("diamond_444q_full_SIGEPH.nc"))
        
        boltztrap = AbipyBoltztrap.from_sigeph(sigeph)

        # get equivalences
        assert boltztrap.rmesh == (5,5,5)
        assert boltztrap.nequivalences == 5
        
        # get coefficients
        assert boltztrap.ncoefficients == 53

        #get results
        boltztrap_results = boltztrap.run()
