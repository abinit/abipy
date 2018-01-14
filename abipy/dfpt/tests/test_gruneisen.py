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
        with abilab.abiopen(abidata.ref_file("mg2si_GRUNS.nc")) as ncfile:
            repr(ncfile); str(ncfile)
            assert ncfile.structure.formula == "Mg2 Si1"
            assert ncfile.iv0 == 1
            #assert len(ncfile.volumes) == 3
            assert not ncfile.params

            d = ncfile.doses
            assert d is ncfile.doses
            assert "wmesh" in d
            assert len(ncfile.phbands_qpath_vol) == 3
            assert d.qpoints.is_ibz

            df = ncfile.to_dataframe()
            assert "grun" in df and "freq" in df

            if self.has_matplotlib():
                assert ncfile.plot_doses(title="DOSes", show=False)
                assert ncfile.plot_doses(with_idos=False, xlims=None, show=False)

                # Arrow up for positive values, down for negative values.
                assert ncfile.plot_phbands_with_gruns(title="bands with gamma markers + DOSes", show=False)
                assert ncfile.plot_phbands_with_gruns(with_doses=None, gamma_fact=2, units="cm-1", match_bands=False)

                plotter = ncfile.get_plotter()
                assert plotter.combiboxplot(show=False)
                assert plotter.animate()

            if self.has_nbformat():
                assert ncfile.write_notebook(nbpath=self.get_tmpname(text=True))
