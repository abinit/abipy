"""Tests for psps module."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata
import abipy.core

from abipy.core.testing import AbipyTest
from abipy.electrons.psps import PspsFile


class PspsFileTestCase(AbipyTest):

    def test_psps_nc_silicon(self):
        """Test PSPS.nc file with Ga.oncvpsp"""
        pseudo = abidata.pseudo("Ga.oncvpsp")

        with pseudo.open_pspsfile(ecut=10) as psps:
            repr(psps); print(psps)
            r = psps.reader
            assert r.usepaw == 0 and r.ntypat == 1
            assert not psps.params

            if self.has_matplotlib():
                psps.plot(what="all", with_qn=True, show=False)
                psps.compare(psps, show=False)
