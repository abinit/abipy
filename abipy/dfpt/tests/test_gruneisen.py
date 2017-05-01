"""Tests for Grunesein module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy import abilab


class GrunsFileTest(AbipyTest):

    def test_gruns_ncfile(self):
        """Testsing GrunsFile."""
        filepath = "/Users/gmatteo/git_repos/abinit_eph/build_gcc/tests/grunesein/run.abo_GRUNS.nc"
        with abilab.abiopen(filepath) as ncfile:
            repr(ncfile); str(ncfile)
            #assert ncfile.structure.formula == "Si2"
            assert ncfile.iv0 == 2

            d = ncfile.doses
            assert d is ncfile.doses
            assert "wmesh" in d
            assert len(ncfile.phbands_qpath_vol) == 3
            assert d.qpoints.is_ibz

            df = ncfile.to_dataframe()
            assert "grun" in df and "freq" in df

            if self.has_matplotlib():
                ncfile.plot_doses(title="DOSes")
                ncfile.plot_doses(with_idos=False, xlims=None)

                # Arrow up for positive values, down for negative values.
                ncfile.plot_phbands_with_gruns(title="bands with gamma markers + DOSes")
                ncfile.plot_phbands_with_gruns(with_doses=None, gamma_fact=2, units="cm-1", match_bands=False)

                plotter = ncfile.get_plotter()
                plotter.combiboxplot()

            if self.has_nbformat():
                assert ncfile.write_notebook(nbpath=self.get_tmpname(text=True))
