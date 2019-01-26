# coding: utf-8
"""Tests for ddk module."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.ddk import DdkFile
from abipy import abilab


#class DdkTest(AbipyTest):

    #def test_ddk_api(self):
    #    """Testing DDK API"""
    #    with abilab.abiopen(abidata.ref_file("gaas_444_dir1_DDK.nc")) as ddk:
    #        repr(ddk); str(ddk)
    #        assert ddk.to_string(verbose=2)

    #        assert ddk.structure.formula == "Ga1 As1"
    #        assert ddk.ebands.nsppol == 1 and ddk.ebands.nspden == 1 and ddk.ebands.nspinor == 1
    #        assert ddk.ebands.nband == 20
    #        #assert ddk.idir == 1 and ddk.ipert == len(ddk.structure) + 1
    #        assert ddk.kptopt == 2
    #        ksamp = ddk.kpoints.ksampling
    #        self.assert_equal(ksamp.kptrlatt_orig.ravel(), [4, 0, 0, 0, 4, 0, 0, 0, 4])
    #        self.assert_equal(ksamp.shifts_orig.ravel(),
    #                [0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5])
    #        assert ddk.kpoints.is_ibz
    #        ddk.kpoints.check_weights()
    #        assert np.all(ddk.hdr["qptn"] == 0)
    #        #assert ddk.xc == "LDA"

    #        assert ddk.params["nband"] == 20

    #        if self.has_nbformat():
    #            assert ddk.write_notebook(nbpath=self.get_tmpname(text=True))


#class DdkAnalyzerTest(AbipyTest):
#
#    def test_ddk_analyzer(self):
#        """Testing DDK analyzer."""
#
#        # TODO: Remove files
#        ddk_paths = abidata.ref_files(
#            "gaas_444_dir1_DDK.nc",
#            "gaas_444_dir2_DDK.nc",
#            "gaas_444_dir3_DDK.nc",
#        )
#
#        with DdksAnalyzer(ddk_paths) as dka:
#            repr(dka); str(dka)
#            assert dka.to_string(verbose=2)
#            assert dka.nsppol == 1 and dka.nspden == 1 and dka.nspinor == 1
