"""Tests for wannier90 module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import os
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest


class TestAbiwanFile(AbipyTest):

    def test_abiwan_without_dis(self):
        """Testing abiwan file without DISENTANGLE."""
        filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "foo_ABIWAN.nc")
        with abilab.abiopen(filepath) as abiwan:
            repr(abiwan); str(abiwan)
            assert abiwan.to_string(verbose=2)
            assert abiwan.structure.formula == "Ga1 As1"
            self.assert_equal(abiwan.nwan, [4])
            assert abiwan.mwan == 4 and abiwan.nntot == 4
            self.assert_equal(abiwan.num_bands_spin, [4])
            self.assert_equal(abiwan.have_disentangled, [False])

            #abiwan.lwindow
            #natom = len(abiwan.structure)
            assert abiwan.wann_centers.shape == (abiwan.nsppol, abiwan.mwan, 3)
            #self.assert_equal(abiwan.wf_centers[1, 20], [-0.866253,  0.866253,  0.866253])
            assert abiwan.wan_spreads.shape == (abiwan.nsppol, abiwan.mwan)
            #self.assert_equal(abiwan.wf_spreads[1, 20], 1.11672024)
            assert len(abiwan.irvec) == len(abiwan.nedeg)
            #abiwan.params

            # Compare input eigenvalues with interpolated values.
            in_eigens = abiwan.ebands.eigens
            for spin in range(abiwan.nsppol):
                for ik, kpt in enumerate(abiwan.kpoints):
                    ews = abiwan.hwan.eval_sk(spin, kpt.frac_coords)
                    self.assert_almost_equal(ews, in_eigens[spin, ik, :self.nwan_spin[spin]])

            ebands_wan = abiwan.interpolate_ebands(line_density=10)
            plotter = abiwan.get_plotter_from_ebands(ebands_wan)

            if self.has_matplotlib():
                assert abiwan.hwan.plot(show=False)
                assert abiwan.plot_centers_spread(show=False)

            if self.has_nbformat():
                assert abiwan.write_notebook(nbpath=self.get_tmpname(text=True))

    #def test_abiwan_with_dis(self):
    #    """Testing abiwan file with DISENTANGLE."""
    #    filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "foo_ABIWAN.nc")
    #    with abilab.abiopen(filepath) as abiwan:
    #        repr(abiwan); str(abiwan)
    #        assert abiwan.to_string(verbose=2)
    #        #assert abiwan.structure.formula == "Ga1 As1"
    #        self.assert_equal(abiwan.nwan, [4])
    #        assert abiwan.mwan == 4
    #        assert abiwan.nntot == 4
    #        self.assert_equal(abiwan.num_bands_spin, [4]
    #        self.assert_equal(abiwan.have_disentangled, [False])

    #        #abiwan.lwindow
    #        ## numpy array (nwan, nstep, ...)
    #        ##natom = len(abiwan.structure)
    #        #nstep = 21
    #        #assert abiwan.wann_centers.shape == (abiwan.nsppol, abiwan.mwan, 3)
    #        #self.assert_equal(abiwan.wf_centers[1, 20], [-0.866253,  0.866253,  0.866253])
    #        #assert abiwan.wan_spreads.shape == (abiwan.nsppol, abiwan.mwan)
    #        assert len(abiwan.irvec) == len(abiwan.nedeg)
    #        #self.assert_equal(abiwan.wf_spreads[1, 20], 1.11672024)
    #        #abiwan.params

    #        Compare input eigenvalues with interpolated values
    #        in_eigens = abiwan.ebands.eigens
    #        for spin in range(abiwan.nsppol):
    #            for ik, kpt in enumerate(abiwan.kpoints):
    #                ews = abiwan.hwan.eval_sk(spin, kpt.frac_coords)
    #                self.assert_almost_equal(ews, in_eigens[spin, ik, :self.nwan_spin[spin]])

    #        #ebands_wan = awiwan.interpolate_ebands(line_density=10)
    #        #plotter = abiwan.get_plotter_from_ebands(ebands_wan)

    #        #if self.has_matplotlib():
    #        #    assert abiwan.hwan.plot(show=False)
    #        #    assert abiwan.plot_centers_spread(show=False)

    #        #if self.has_nbformat():
    #        #    assert abiwan.write_notebook(nbpath=self.get_tmpname(text=True))


    def test_abiwan_robot(self):
        """Testing abiwan file with DISENTANGLE."""
        filepaths = [os.path.join(abidata.dirpath, "refs", "wannier90", "foo_ABIWAN.nc")]

        robot = abilab.AbiwanRobot.from_files()
        assert repr(robot); assert str(robot)
        assert robot.to_string(verbose==2)
        assert len(robot.abifiles) == 1
        assert robot.EXT == "ABIWAN"

	# Get pandas dataframe.
        df = robot.get_dataframe()
        self.assert_equal(df["ecut"].values, 6.0)

        if self.has_matplotlib():
            assert ebands_plotter.gridplot(show=False)

        if self.has_nbformat():
            robot.write_notebook(nbpath=self.get_tmpname(text=True))

        robot.close()
