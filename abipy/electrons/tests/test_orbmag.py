# coding: utf-8
"""Tests for orbmag module."""
import pytest
import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.orbmag import OrbmagAnalyzer

root = "/Users/giantomassi/git_repos/ABIPY_WORKS/JOE_ORBMAG"

class OrbmagTest(AbipyTest):

    @pytest.mark.xfail(condition=not os.path.isdir(root), reason=f"{root=} does not exist")
    def test_orbmag_analyzer(self):
        """Testing OrbmagAnalyzer"""
        filepaths = [os.path.join(root, s) for s in ["gso_DS12_ORBMAG.nc", "gso_DS22_ORBMAG.nc", "gso_DS32_ORBMAG.nc"]]

        with OrbmagAnalyzer(filepaths) as orban:
            repr(orban); str(orban)
            #orban.to_string(verbose=1)
            assert orban.structure.formula == "Al1 P1"
            assert orban.mband == 4
            assert orban.nkpt == 64
            assert orban.nsppol == 1
            assert orban.ndir == 3
            assert orban.has_full_bz

            ngkpt, shifts = orban.ngkpt_and_shifts
            self.assert_equal(ngkpt, [4, 4, 4])
            assert shifts.shape == (1, 3)
            self.assert_equal(shifts.flatten(), [0, 0, 0])

            self.assert_almost_equal(orban.get_value("isotropic", spin=0, ikpt=0, band=0), -1.2001081154385034)
            self.assert_almost_equal(orban.get_value("isotropic", spin=0, ikpt=1, band=1), -1.8096595705860796)

            ref_omlamb = np.zeros((3, 3))
            self.assert_almost_equal(orban.get_omlamb(), ref_omlamb)

            for report_type in ['T', 'B', 'TB']:
                orban.report_eigvals(report_type=report_type)

            orb = orban.orb_files[0]
            assert orb.to_string(verbose=2)
            params = orb.params
            #assert params["orban_ntau"] == 6
            target_atom, nucdipmom = orb.target_atom_nucdipmom
            assert target_atom == 0
            self.assert_almost_equal(nucdipmom, [1, 0, 0])

            if self.has_matplotlib():
                orban.plot_fatbands(os.path.join(root, "bandso_DS1_GSR.nc"), show=False)
