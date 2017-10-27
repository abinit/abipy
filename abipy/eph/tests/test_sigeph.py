"""Tests for eph module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest


class SigEPhFileTest(AbipyTest):

    def test_sigeph_file(self):
        """Tests for SigEPhFile."""
        return
        ncfile = abilab.abiopen(abidata.ref_file("al_888k_161616q_EPH.nc"))
        repr(ncfile); str(ncfile)
        assert ncfile.to_string(verbose=2)
        assert ncfile.nsppol == 1 and ncfile.nspden == 1 and ncfile.nspinor == 1
        assert ncfile.ebands.kpoints.is_ibz
        self.assert_equal(ncfile.ebands.kpoints.ksampling.mpdivs, [8, 8, 8])
        assert ncfile.phbands.qpoints.is_path
        assert ncfile.phbands.qpoints.ksampling is None

        if self.has_matplotlib():
            # Test A2f plot methods
            #assert a2f.plot(show=False)

            # Test Eph plot methods.
            assert ncfile.plot(show=False)
            assert ncfile.plot_eph_strength(show=False)
            assert ncfile.plot_with_a2f(show=False)

        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))


    def test_sigeph_robot(self):
        """Tests for SigEPhRobot."""
        return

        with abilab.SigEPhRobot(filepaths) as robot:
            robot.to_string(verbose=2)

            # Test plot methods
            if self.has_matplotlib():
                assert robot.plot_linopt_convergence(show=False)
                assert robot.plot_shg_convergence(show=False)
                assert robot.plot_leo_convergence(show=False)

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))
