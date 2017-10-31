"""Tests for eph module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy import abilab


class SigEPhFileTest(AbipyTest):

    def test_sigeph_file(self):
        """Tests for SigEPhFile."""
        #return
        ncfile = abilab.abiopen(abidata.ref_file("test_SIGEPH.nc"))
        repr(ncfile); str(ncfile)
        assert ncfile.to_string(verbose=2)
        assert ncfile.nsppol == 1 and ncfile.nspden == 1 and ncfile.nspinor == 1
        assert ncfile.structure.formula == "C2"
        assert ncfile.ebands.kpoints.is_ibz
        self.assert_equal(ncfile.ebands.kpoints.ksampling.mpdivs, [4, 4, 4])

        assert ncfile.nkcalc == 2 #and ncfile.nqbz == 64
        self.assert_equal(ncfile.ngqpt.flatten(), [4, 4, 4])
        assert ncfile.symsigma == 0
        assert ncfile.ntemp ==  6
        assert ncfile.nband == 54
        #self.assert_almost_equal(ncfile.eta, 0.001)
        #assert len(self.mu_e) == self.ntemp

        assert ncfile.ks_dirgaps.shape == (ncfile.nsppol, ncfile.nkcalc)
        assert ncfile.qp_dirgaps_t.shape == (ncfile.nsppol, ncfile.nkcalc, ncfile.ntemp)

        assert len(ncfile.sigma_kpoints) == 2
        assert not ncfile.sigma_kpoints.is_ibz and not ncfile.sigma_kpoints.is_path
        assert ncfile.sigma_kpoints[0] == [0, 0, 0]
        assert ncfile.sigma_kpoints[1] == [0.5, 0, 0]
        assert ncfile.bstart_sk.shape == (ncfile.nsppol, ncfile.nkcalc)
        self.assert_equal(ncfile.bstart_sk.flatten(), [0, 0])
        self.assert_equal(ncfile.nbcalc_sk.flatten(), [8, 8])
        assert ncfile.bstart_sk.shape == ncfile.nbcalc_sk.shape
        self.assert_equal(ncfile.bstop_sk.flatten(), [8, 8])
        self.assert_almost_equal(ncfile.tmesh[0], 5)
        # FIXME
        #self.assert_almost_equal(ncfile.zcut, )
        #assert self.nbsum ==

        assert ncfile.sigkpt2index(0) == 0
        assert ncfile.sigkpt2index([0.5, 0.0, 0.0]) == 1

        data_sk = ncfile.get_dataframe_sk(spin=0, sig_kpoint=[0.5, 0.0, 0.0])
        #assert

        data = ncfile.get_dataframe()
        #assert

        if self.has_matplotlib():
            # Test SigEph plot methods.
            assert ncfile.plot_dirgaps_t(show=False)

        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_sigeph_robot(self):
        """Tests for SigEPhRobot."""
        #return
        filepaths = [
            abidata.ref_file("test_SIGEPH.nc"),
        ]
        with abilab.SigEPhRobot.from_files(filepaths) as robot:
            robot.add_file("same_file", filepaths[0])
            repr(robot); str(robot)
            robot.to_string(verbose=2)
            assert len(robot) == 2

            # Test plot methods
            #if self.has_matplotlib():
            #    assert robot.plot_qp_convergence(show=False)

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))
