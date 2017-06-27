"""Tests for Fold2Bloch module."""
from __future__ import division, print_function, unicode_literals, absolute_import

import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.electrons import Fold2Bloch


class Fold2BlochTest(AbipyTest):

    def test_api(self):
        """Test Fold2Bloch API"""
        filepath = ...
        with Fold2Bloch(filepath) as fb:
            repr(fb); str(fb)
            assert fb.to_string(verbose=1)
            assert fb.nsppol == 1 and fb.nspden == 1
            assert fb.nss == 1 and fb.uf_nkpt == 1
            assert fb.ebands is None
            assert fb.structure.formula == "H6"

            if self.has_nbformat():
                assert fb.write_notebook(nbpath=self.get_tmpname(text=True))

            if self.has_matplotlib():
                assert fb.plot_unfolded(klabels=None, ylims=None, dist_tol=1e-12, verbose=1,
                                        colormap="afmhot", facecolor="black")
