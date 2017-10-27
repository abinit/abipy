"""Tests for eph module."""
from __future__ import print_function, division, unicode_literals, absolute_import

#import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest


class EphFileTest(AbipyTest):

    def test_eph_file(self):
        """Tests for EphFile."""
        ncfile = abilab.abiopen(abidata.ref_file("al_888k_161616q_EPH.nc"))
        repr(ncfile); str(ncfile)
        assert ncfile.to_string(verbose=2)
        assert ncfile.nsppol == 1 and ncfile.nspden == 1 and ncfile.nspinor == 1
        assert ncfile.ebands.kpoints.is_ibz
        self.assert_equal(ncfile.ebands.kpoints.ksampling.mpdivs, [8, 8, 8])
        assert ncfile.phbands.qpoints.is_path
        assert ncfile.phbands.qpoints.ksampling is None

        # Test A2f(w) function.
        #a2f = ncfile.a2f
        #repr(a2f); str(a2f)
        #a2f.to_string(verbose=2)
        #assert a2f.nsppol == 1 and a2f.nmodes == 3
        #assert a2f.mesh
        #assert a2f.values_spin
        #assert a2f.get_momentum(order=1) ==
        #assert a2f.get_mcmillan_Tc(mustar=0.8) ==

        if self.has_matplotlib():
            # Test A2f plot methods
            #assert a2f.plot(show=False)

            # Test Eph plot methods.
            assert ncfile.plot(show=False)
            assert ncfile.plot_eph_strength(show=False)
            assert ncfile.plot_with_a2f(show=False)

        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))


from abipy.eph.eph import EphRobot

class EphRobotTest(AbipyTest):

    def test_eph_robot(self):
        """Test EphRobot."""
        return
        files = abidata.ref_files(
                "al_888k_161616q_EPH.nc",
                "al_888k_161616q_EPH.nc",
        )
        with EphRobot.from_files(files) as robot:
            assert len(robot) == 2
            repr(robot); str(robot)
            robot.to_string(verbose=2)
            #assert [t[2] for t in robot.sortby("nkpt")] == [10, 60, 182]

            # Test plot methods
            if self.has_matplotlib():
                assert robot.plot_lambda_convergence(show=False)
                #assert robot.plot_a2f_convergence(show=False)

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))
