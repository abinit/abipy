"""Tests for electrons.bse module"""
from __future__ import print_function, division

import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.bse import *


class TestMDF_Reader(AbipyTest):

    def test_MDF_reading(self):
        """Test MdfReader."""
        with MdfReader(abidata.ref_file("tbs_4o_DS2_MDF.nc")) as r:
            assert len(r.wmesh) == r.read_dimvalue("number_of_frequencies")
            assert len(r.qpoints) == r.read_dimvalue("number_of_qpoints")

            exc_mdf = r.read_exc_mdf()
            rpanlf_mdf = r.read_rpanlf_mdf()
            gwnlf_mdf = r.read_gwnlf_mdf()

            if self.has_matplotlib():
                exc_mdf.plot(show=False)

            # Test Plotter.
            plotter = MdfPlotter()
            plotter.add_mdf("EXC", exc_mdf)
            plotter.add_mdf("KS-RPA", rpanlf_mdf)
            plotter.add_mdf("GW-RPA", gwnlf_mdf)
            if self.has_matplotlib():
                plotter.plot(show=False)

    def test_mdf_api(self):
        """Test MdfFile API"""
        with MdfFile(abidata.ref_file("tbs_4o_DS2_MDF.nc")) as mdf_file:
            print(mdf_file)
            assert len(mdf_file.structure) == 2

            exc_tsr = mdf_file.get_tensor("exc")
            rpa_tsr = mdf_file.get_tensor("rpa")
            gw_tsr = mdf_file.get_tensor("gwrpa")

            rpa = mdf_file.get_mdf("rpa")
            str(rpa)
            rpa.to_string(with_info=True)
            assert rpa.num_qpoints == 6
            assert rpa.num_qpoints == len(rpa.qfrac_coords)
            assert mdf_file.qpoints == rpa.qpoints
            assert np.all(mdf_file.qfrac_coords == rpa.qfrac_coords)
            assert  mdf_file.params.get("nsppol") == 1

            if self.has_matplotlib():
                # Test plot_mdfs
                mdf_file.plot_mdfs(cplx_mode="Im", mdf_type="all", qpoint=None, show=False)
                mdf_file.plot_mdfs(cplx_mode="RE", mdf_type="all", qpoint=0, show=False)
                mdf_file.plot_mdfs(cplx_mode="re", mdf_type="all", qpoint=mdf_file.qpoints[0], show=False)

            if self.has_nbformat():
                mdf_file.write_notebook(nbpath=self.get_tmpname(text=True))


class MultipleMdfPlotterTest(AbipyTest):

    def test_multiplemdf_plotter(self):
        """Testing MultipleMdfPlotter."""
        mdf_paths = abidata.ref_files("si_444_MDF.nc", "si_666_MDF.nc", "si_888_MDF.nc")
        plotter = MultipleMdfPlotter()
        for f in mdf_paths:
            plotter.add_mdf_file(f, f)
        repr(plotter)
        str(plotter)
        assert plotter._can_use_basenames_as_labels()
        assert len(plotter._get_qpoints()) == 6

        if self.has_matplotlib():
            xlim, ylim = (2, 3), (1, None)
            plotter.plot(mdf_type="exc", qview="avg", xlim=xlim, ylim=ylim, show=False)
            plotter.plot(mdf_type="exc", qview="all", show=False)
            #plotter.plot_mdftypes(qview="avg", xlim=xlim, ylim=ylim, show=False)
            #plotter.plot_mdftypes(qview="all", xlim=xlim, ylim=ylim, show=False)

        if self.has_ipywidgets():
            assert plotter.ipw_select_plot() is not None
