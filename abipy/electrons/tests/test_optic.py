# coding: utf-8
"""Tests for optic module."""
from __future__ import division, print_function, unicode_literals, absolute_import

#import numpy as np
#import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy import abilab


class OpticTest(AbipyTest):

    def test_opticfile(self):
        """Testing OpticFile"""
        path = "/Users/gmatteo/git_repos/abinit_quick_prs/build_gcc/tests/Test_suite/tutorespfn_toptic_3-toptic_4/toptic_4_OPTIC.nc"
        with abilab.abiopen(path) as optic:
            repr(optic); str(optic)
            assert optic.to_string(verbose=2)
            assert optic.structure.formula == "Ga1 As1"
            assert optic.kptopt == 2
            self.assert_almost_equal(optic.broadening, 0.002)
            self.assert_almost_equal(optic.domega, 0.0003)
            self.assert_almost_equal(optic.maxomega, 0.3)
            self.assert_almost_equal(optic.scissor, 0.000)
            self.assert_almost_equal(optic.tolerance, 0.002)
            assert optic.ntemp == 1

            assert optic.reader.linopt_computed_components == ["xx"]
            assert not optic.reader.shg_computed_components
            assert not optic.reader.leo_computed_components
            assert not optic.reader.leo2_computed_components

            # Test plot methods
            if self.has_matplotlib():
                assert optic.plot_linear_optic(show=False)
                #assert optic.plot_shg(show=False)
                #assert optic.plot_leo(show=False)
                #assert optic.plot_leo2(show=False)

            if self.has_nbformat():
                optic.write_notebook(nbpath=self.get_tmpname(text=True))
