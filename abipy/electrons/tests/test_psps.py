"""Tests for psps module."""
import numpy as np
import abipy.data as abidata
import abipy.core

from abipy.core.testing import AbipyTest
from abipy.electrons.psps import PspsFile, PspsRobot


class PspsFileTestCase(AbipyTest):

    def test_psps_nc_silicon(self):
        """Test PSPS.nc file with Ga.oncvpsp"""
        pseudo = abidata.pseudo("Ga.oncvpsp")

        with pseudo.open_pspsfile(ecut=10) as psps:
            repr(psps); print(psps)
            r = psps.reader
            assert r.usepaw == 0 and r.ntypat == 1
            assert not psps.params


            robot = PspsRobot.from_files([psps.filepath])
            repr(psps); print(psps)

            if self.has_matplotlib():
                # psps plots.
                assert psps.plot(what="all", with_qn=True, show=False)

                # robot plots.
                assert robot.plot_tcore_rspace(ders=(0, 1, 2, 3), with_qn=0, scale=None, fontsize=8, show=False)
                assert robot.plot_tcore_qspace(ders=(0, 1), with_qn=0, scale=None, fontsize=8, show=False)
                assert robot.plot_vlocq(ders=(0, 1), with_qn=0, with_fact=True, scale=None, fontsize=8, show=False)
                assert robot.plot_ffspl(ecut_ffnl=None, ders=(0, 1), with_qn=0, scale=None, fontsize=8, show=False)
