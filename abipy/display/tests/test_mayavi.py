"""Tests for psps module."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata
import abipy.core

from abipy.core.testing import AbipyTest


class MayaviTest(AbipyTest):

    def test_psps_nc_silicon(self):
        """Test mayavi toolkit."""
        #if self.has_matplotlib():
        #    psps.plot(what="all", with_qn=True, show=False)
        #    psps.compare(psps, show=False)
