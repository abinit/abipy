"""Tests for electrons.arpes module"""
from __future__ import print_function, division, unicode_literals, absolute_import

#import os
#import numpy as np
import abipy.data as abidata

#from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.electrons.arpes import ArpesPlotter


class TestArpesPlotter(AbipyTest):

    def test_arpes_plotter_api(self):
        """Testing ArpesPlotter API."""
        path = abidata.ref_file("si_nscf_GSR.nc")
        plotter = ArpesPlotter.model_from_ebands(path)
        repr(plotter); str(plotter)
        assert plotter.to_string(verbose=2)

        if self.has_matplotlib():
            assert plotter.plot_ekmap_itemp(itemp=0, estep=0.05, show=False)
            assert plotter.plot_ekmap_temps(temp_inds=range(plotter.ntemp), show=False)
            assert plotter.plot_ak_vs_temp(show=False)
            assert plotter.plot_3dlines(itemp=0, estep=0.05, band_inds=[1, 2, 3], show=False)
            assert plotter.plot_surface(itemp=0, estep=0.05, show=False)

        if self.has_nbformat():
            assert plotter.write_notebook(nbpath=self.get_tmpname(text=True))
