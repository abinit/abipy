"""Tests for the transport module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import collections
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy import abilab


class TransportFileTest(AbipyTest):

    def test_transportfile(self):
        """Test abinit transport file"""

        si_transport = abilab.abiopen('sio_DS1_TRANSPORT.nc')
        si_transport.get_mobility_mu(0,0)
