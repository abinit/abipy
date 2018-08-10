"""Tests for wannier90 module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import os
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest


class TestAbiwanFile(AbipyTest):

    def test_abiwan_without_dis(self):
        """Testing ABIWAN file without DISENTANGLE procedure."""
        filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "tutoplugs_tw90_1", "tw90_1o_DS2_ABIWAN.nc")
        with abilab.abiopen(filepath) as abiwan:
            repr(abiwan); str(abiwan)
            assert abiwan.to_string(verbose=2)
            assert abiwan.structure.formula == "Si2"
            self.assert_equal(abiwan.nwan_spin, [4])
            assert abiwan.mwan == 4 and abiwan.nntot == 8
            self.assert_equal(abiwan.num_bands_spin, [4])
            self.assert_equal(abiwan.have_disentangled_spin, [False])

            #abiwan.lwindow
            #natom = len(abiwan.structure)
            assert abiwan.wann_centers.shape == (abiwan.nsppol, abiwan.mwan, 3)
            #self.assert_equal(abiwan.wf_centers[1, 20], [-0.866253,  0.866253,  0.866253])
            #assert abiwan.wan_spreads.shape == (abiwan.nsppol, abiwan.mwan)
            #self.assert_equal(abiwan.wf_spreads[1, 20], 1.11672024)
            assert len(abiwan.irvec) == len(abiwan.ndegen)
            #abiwan.params

            # Compare input eigenvalues with interpolated values.
            in_eigens = abiwan.ebands.eigens
            for spin in range(abiwan.nsppol):
                n = abiwan.nwan_spin[spin]
                for ik, kpt in enumerate(abiwan.kpoints):
                    ews = abiwan.hwan.eval_sk(spin, kpt.frac_coords)
                    self.assert_almost_equal(ews[:n], in_eigens[spin, ik, :n])

            ebands_kmesh = abiwan.interpolate_ebands(ngkpt=(4, 4, 4))
            assert ebands_kmesh.kpoints.is_ibz

            ebands_wan = abiwan.interpolate_ebands(line_density=5)
            assert ebands_wan.kpoints.is_path
            plotter = abiwan.get_plotter_from_ebands(ebands_wan)

            if self.has_matplotlib():
                assert abiwan.hwan.plot(show=False)
                assert plotter.combiplot(show=False)
                #assert abiwan.plot_centers_spread(show=False)

            if self.has_nbformat():
                assert abiwan.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_abiwan_with_dis(self):
        """Testing ABIWAN file with DISENTANGLE."""
        filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "tutoplugs_tw90_4", "tw90_4o_DS3_ABIWAN.nc")
        with abilab.abiopen(filepath) as abiwan:
            repr(abiwan); str(abiwan)
            assert abiwan.to_string(verbose=2)
            #assert abiwan.structure.formula == "Si 2"
            self.assert_equal(abiwan.nwan_spin, [8])
            assert abiwan.mwan == 8 and abiwan.nntot == 8
            self.assert_equal(abiwan.num_bands_spin, [14])
            self.assert_equal(abiwan.have_disentangled_spin, [True])

            #abiwan.lwindow
            ## numpy array (nwan, nstep, ...)
            ##natom = len(abiwan.structure)
            #nstep = 21
            #assert abiwan.wann_centers.shape == (abiwan.nsppol, abiwan.mwan, 3)
            #self.assert_equal(abiwan.wf_centers[1, 20], [-0.866253,  0.866253,  0.866253])
            #assert abiwan.wan_spreads.shape == (abiwan.nsppol, abiwan.mwan)
            assert len(abiwan.irvec) == len(abiwan.ndegen)
            #self.assert_equal(abiwan.wf_spreads[1, 20], 1.11672024)
            #abiwan.params

            # Compare input eigenvalues with interpolated values
            in_eigens = abiwan.ebands.eigens
            for spin in range(abiwan.nsppol):
                # Remember that that we have disentangled so exclude last 2
                n = abiwan.nwan_spin[spin] - 2
                for ik, kpt in enumerate(abiwan.kpoints):
                    ews = abiwan.hwan.eval_sk(spin, kpt.frac_coords)
                    self.assert_almost_equal(ews[:n], in_eigens[spin, ik, :n])

            ebands_wan = abiwan.interpolate_ebands(line_density=3)
            plotter = abiwan.get_plotter_from_ebands(ebands_wan)

            if self.has_matplotlib():
                assert abiwan.hwan.plot(show=False)
                assert plotter.combiplot(show=False)

            if self.has_nbformat():
                assert abiwan.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_abiwan_robot(self):
        """Testing abiwan file with DISENTANGLE."""
        filepaths = [
            os.path.join(abidata.dirpath, "refs", "wannier90", "tutoplugs_tw90_1", "tw90_1o_DS2_ABIWAN.nc"),
            os.path.join(abidata.dirpath, "refs", "wannier90", "tutoplugs_tw90_4", "tw90_4o_DS3_ABIWAN.nc"),
        ]

        robot = abilab.AbiwanRobot.from_files(filepaths)
        assert repr(robot); assert str(robot)
        assert robot.to_string(verbose=2)
        assert len(robot.abifiles) == 2
        assert robot.EXT == "ABIWAN"

	# Get pandas dataframe.
        df = robot.get_dataframe()
        #self.assert_equal(df["ecut"].values, 6.0)

        plotter = robot.get_interpolated_ebands_plotter(line_density=3)

        if self.has_matplotlib():
            assert robot.plot_hwanr(show=False)
            assert plotter.gridplot(show=False)

        if self.has_nbformat():
            assert robot.write_notebook(nbpath=self.get_tmpname(text=True))

        robot.close()
