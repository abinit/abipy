"""Tests for wannier90 module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import os
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.wannier90 import Wannier90Input


class TestWinFile(AbipyTest):

    def test_win(self):
        """Testing win API"""
        filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "example01_gaas.wout")
        win = Wannier90Input.from_abinit_file(filepath)
        repr(win); str(win)
        assert win.to_string(verbose=2)
        assert win.structure.formula == "Ga1 As1"

        #if self.has_matplotlib():
        #    assert wout.plot(show=False)
        #    assert wout.plot_centers_spread(show=False)

        #if self.has_nbformat():
        #    assert wout.write_notebook(nbpath=self.get_tmpname(text=True))
