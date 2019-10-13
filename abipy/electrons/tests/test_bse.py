"""Tests for electrons.bse module"""
import os
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.electrons.bse import *


class TestMDF_Reader(AbipyTest):

    def test_MDF_reading(self):
        """Testing MdfReader."""
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
            with self.assertRaises(ValueError):
                plotter.add_mdf("GW-RPA", gwnlf_mdf)
            repr(plotter); str(plotter)

            if self.has_matplotlib():
                plotter.plot(show=False)

            #if self.has_ipywidgets():
            #    assert plotter.ipw_plot() is not None

    def test_mdf_api(self):
        """Test MdfFile API"""
        with MdfFile(abidata.ref_file("tbs_4o_DS2_MDF.nc")) as mdf_file:
            repr(mdf_file); str(mdf_file)
            assert len(mdf_file.structure) == 2
            assert mdf_file.params["nsppol"] == 1

            exc_tsr = mdf_file.get_tensor("exc")
            rpa_tsr = mdf_file.get_tensor("rpa")
            gw_tsr = mdf_file.get_tensor("gwrpa")

            rpa = mdf_file.get_mdf("rpa")
            repr(rpa); str(rpa)
            rpa.to_string(with_info=True, verbose=2)
            assert rpa.num_qpoints == 6
            assert rpa.num_qpoints == len(rpa.qfrac_coords)
            assert mdf_file.qpoints == rpa.qpoints
            assert np.all(mdf_file.qfrac_coords == rpa.qfrac_coords)
            assert  mdf_file.params.get("nsppol") == 1

            tensor_exc = mdf_file.get_tensor("exc")
            tensor_exc.symmetrize(mdf_file.structure)

            if self.has_matplotlib():
                # Test plot_mdfs
                mdf_file.plot_mdfs(cplx_mode="Im", mdf_type="all", qpoint=None, show=False)
                mdf_file.plot_mdfs(cplx_mode="RE", mdf_type="all", qpoint=0, show=False)
                mdf_file.plot_mdfs(cplx_mode="re", mdf_type="all", qpoint=mdf_file.qpoints[0], show=False)
                tensor_exc.plot(title="Si macroscopic dielectric tensor (Reduced coord)")
                tred = tensor_exc.to_array(red_coords=True)
                assert len(tred) == 300
                tcart = tensor_exc.to_array(red_coords=False)
                assert len(tcart) == 300

            if self.has_nbformat():
                mdf_file.write_notebook(nbpath=self.get_tmpname(text=True))


class MultipleMdfPlotterTest(AbipyTest):

    def test_multiplemdf_plotter(self):
        """Testing MultipleMdfPlotter."""
        mdf_paths = abidata.ref_files("si_444_MDF.nc", "si_666_MDF.nc", "si_888_MDF.nc")
        plotter = MultipleMdfPlotter()
        for f in mdf_paths:
            plotter.add_mdf_file(f, f)
        repr(plotter); str(plotter)
        assert plotter._can_use_basenames_as_labels()
        assert len(plotter._get_qpoints()) == 6

        if self.has_matplotlib():
            xlims, ylims = (2, 3), (1, None)
            plotter.plot(mdf_type="exc", qview="avg", xlims=xlims, ylims=ylims, show=False)
            plotter.plot(mdf_type="exc", qview="all", show=False)
            #plotter.plot_mdftypes(qview="avg", xlims=xlims, ylims=ylims, show=False)
            #plotter.plot_mdftypes(qview="all", xlims=xlims, ylims=ylims, show=False)

        if self.has_ipywidgets():
            assert plotter.ipw_select_plot() is not None


class MdfRobotTest(AbipyTest):

    def test_mdf_robot(self):
        """Testing MDF robot."""
        robot = abilab.MdfRobot.from_dir(os.path.join(abidata.dirpath, "refs", "si_bse_kpoints"))
        assert len(robot) == 3

        robot = abilab.MdfRobot()
        robot.scan_dir(os.path.join(abidata.dirpath, "refs", "si_bse_kpoints"))
        assert len(robot) == 3
        repr(robot); str(robot)

        df = robot.get_dataframe(with_geo=True)
        assert df is not None
        df_params = robot.get_params_dataframe()
        assert "nsppol" in df_params

        plotter = robot.get_multimdf_plotter()
        if self.has_matplotlib():
            assert plotter.plot(show=False)
            assert robot.plot(show=False)
            #robot.plot_conv_mdf(self, hue, mdf_type="exc_mdf", **kwargs):

        if self.has_nbformat():
            robot.write_notebook(nbpath=self.get_tmpname(text=True))

        robot.close()
