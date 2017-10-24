# coding: utf-8
"""Tests for ddk module."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.ddk import DdkFile, DdksAnalyzer
from abipy import abilab


class DdkTest(AbipyTest):

    def test_ddk_api(self):
        """Testing DDK API"""
        path = "/Users/gmatteo/git_repos/abinit_quick_prs/build_gcc/tests/Test_suite/tutorespfn_toptic_3-toptic_4/toptic_3o_DS4_DDK.nc"

        with abilab.abiopen(path) as ddk:
            repr(ddk); str(ddk)
            assert ddk.to_string(verbose=2)

            assert ddk.ebands.nsppol == 1 and ddk.ebands.nspden == 1 and ddk.ebands.nspinor == 1
            assert ddk.ebands.nband == 9
            assert ddk.idir == 1 and ddk.ipert == len(ddk.structure) + 1
            assert ddk.kptopt == 2
            assert ddk.kpoints.is_ibz
            ddk.kpoints.check_weights()
            assert np.all(ddk.hdr["qptn"] == 0)
            #assert ddk.xc == "LDA"

            #if self.has_matplotlib():
            #    assert f.plot_freq(gvec1=0, gvec2=None, waxis="real", cplx_mode="re-im", show=False)
            #    assert f.plot_freq(gvec1=[0, 0, 0], gvec2=[1, 0, 0], waxis="imag", cplx_mode="re-im", show=False)

            #if self.has_nbformat():
            #    assert ddk.write_notebook(nbpath=self.get_tmpname(text=True))


class DdkAnalyzerTest(AbipyTest):

    def test_ddk_analyzer(self):
        """Testing DDK analyzer."""

        ddk_paths = [
"/Users/gmatteo/git_repos/abinit_quick_prs/build_gcc/tests/Test_suite/tutorespfn_toptic_3-toptic_4/toptic_3o_DS4_DDK.nc",
"/Users/gmatteo/git_repos/abinit_quick_prs/build_gcc/tests/Test_suite/tutorespfn_toptic_3-toptic_4/toptic_3o_DS5_DDK.nc",
"/Users/gmatteo/git_repos/abinit_quick_prs/build_gcc/tests/Test_suite/tutorespfn_toptic_3-toptic_4/toptic_3o_DS6_DDK.nc"]

        with DdksAnalyzer(ddk_paths) as dka:
            repr(dka); str(dka)
            assert dka.to_string(verbose=2)
            assert dka.nsppol == 1 and dka.nspden == 1 and dka.nspinor == 1

            #dka.v_skb
            #r = dka.get_doses()
            #assert np.all(r.edos.mesh == e.vdos.mesh)

            #if self.has_matplotlib():
            #    assert dka.plot_vdos(show=False)

            #if self.has_nbformat():
            #    assert dka.write_notebook(nbpath=self.get_tmpname(text=True))
