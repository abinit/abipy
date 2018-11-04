"""Tests for msqdos module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy import abilab


class MsqdTest(AbipyTest):

    def test_from_ddb(self):

        filepath = "~/git_repos/abipy/lifetimes/alphaSiO2/run.abo_PHDOS.nc"
        with abilab.abiopen(filepath) as phdos:
            msqd_dos = phdos.msqd_dos

        print(msqd_dos)
        #for fmt in ("cartesian", "cif", "ustar", "beta", "B"):
        for fmt in ("cartesian", "cif"):
            #print("Format:", fmt)
            df = msqd_dos.get_dataframe(temp=100, view="all", select_symbols="Si", fmt=fmt)
            abilab.print_dataframe(df)

        #msqd_dos.write_cif_file("foo.cif", temp=300)

        if self.has_matplotlib():
            assert msqd_dos.plot(show=False)
            assert msqd_dos.plot_tensor(show=False)
            assert msqd_dos.plot_uiso(show=False)
            #assert msqd_dos.plot(view="all")
            #assert msqd_dos.plot_tensor(view="all")
            #assert msqd_dos.plot_uiso(view="all")
            #assert msqd_dos.plot_uiso(view="all", what="vel")
