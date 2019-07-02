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

        with abilab.abiopen(abidata.ref_file('sio_DS1_TRANSPORT.nc')) as si_transport:
            assert repr(si_transport); assert str(si_transport); assert si_transport.to_string(verbose=2)
            si_transport.get_mobility_mu(0, 0)

            if self.has_matplotlib():
                assert si_transport.plot_dos(title="default values", show=False)
                assert si_transport.plot_vvdos(colormap="viridis", component="yy", show=False)
                assert si_transport.plot_mobility(colormap="viridis", component="yy", show=False)

            if self.has_nbformat():
                assert si_transport.write_notebook(nbpath=self.get_tmpname(text=True))
