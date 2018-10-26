"""Tests for Fold2Bloch module."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.fold2bloch import Fold2BlochNcfile


class Fold2BlochTest(AbipyTest):

    def test_fold2bloch_h6(self):
        """Test Fold2Bloch API"""
        with Fold2BlochNcfile(abidata.ref_file("h6_FOLD2BLOCH.nc")) as fb:
            repr(fb); str(fb)
            assert fb.to_string(verbose=1)
            assert fb.nsppol == 1 and fb.nspden == 1
            assert fb.nss == 1 and fb.uf_nkpt == 252
            assert fb.ebands is not None
            assert fb.structure.formula == "H6"
            assert fb.params["nspden"] == fb.nspden
            r = fb.get_spectral_functions()
            nw = len(r.mesh)
            assert r.sfw.shape == (fb.nss, fb.uf_nkpt, nw)
            assert r.int_sfw.shape == (fb.nss, fb.uf_nkpt, nw)

            if self.has_nbformat():
                assert fb.write_notebook(nbpath=self.get_tmpname(text=True))

            if self.has_matplotlib():
                kbounds = [0, 1/2, 0, 0, 0, 0, 0, 0, 1/2]
                klabels = ["Y", r"$\Gamma$", "X"]
                assert fb.plot_unfolded(kbounds, klabels, dist_tol=1e-12, verbose=1,
                                        colormap="afmhot", facecolor="black")
