"""Tests for the transport module."""

import abipy.data as abidata
from abipy import abilab
from abipy.core.testing import AbipyTest


class TransportFileTest(AbipyTest):
    def test_transportfile(self):
        """Test abinit transport file"""
        with abilab.abiopen(abidata.ref_file("sio_DS1_TRANSPORT.nc")) as si_transport:
            assert repr(si_transport)
            assert str(si_transport)
            assert si_transport.to_string(verbose=2)
            si_transport.get_mobility_mu(0, 0)

            if self.has_matplotlib():
                assert si_transport.plot_edos(title="default values", show=False)
                assert si_transport.plot_vvtau_dos(colormap="viridis", component="yy", show=False)
                assert si_transport.plot_mobility(colormap="viridis", component="yy", show=False)

            if self.has_nbformat():
                assert si_transport.write_notebook(nbpath=self.get_tmpname(text=True))
