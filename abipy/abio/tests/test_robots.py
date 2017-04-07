# coding: utf-8
"""Test for Robots"""
from __future__ import unicode_literals, division, print_function

import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.core.testing import AbipyTest
from abipy.abio.robots import *


class RobotTest(AbipyTest):

    def test_base_robot_class(self):
        """Testing base robot class"""
        # With context.
        with Robot() as robot:
            repr(robot)
            str(robot)
            assert len(robot) == 0 and not robot.exceptions

        assert Robot.class_for_ext("DDB") is DdbRobot
        with self.assertRaises(ValueError):
            Robot.for_ext("foo")

    def test_gsr_robot(self):
        """Testing GSR robot"""
        gsr_path = abidata.ref_file("si_scf_GSR.nc")
        robot = Robot.for_ext("GSR")()
        robot.add_file("gsr0", gsr_path)
        assert len(robot.ncfiles) == 1
        assert robot.EXT == "GSR"
        repr(robot)
        str(robot)

	# Cannot have same label
        with self.assertRaises(ValueError):
            robot.add_file("gsr0", gsr_path)

        assert len(robot) == 1 and not robot.exceptions

        robot.add_file("gsr1", abilab.abiopen(gsr_path))
        assert len(robot) == 2
        robot.show_files()

        ebands_plotter = robot.get_ebands_plotter()
        edos_plotter = robot.get_edos_plotter()

        if self.has_matplotlib():
            ebands_plotter.gridplot(show=False)
            edos_plotter.gridplot(show=False)

	# Get pandas dataframe.
        df = robot.get_dataframe()
        assert "energy" in df
        self.assert_equal(df["ecut"].values, 6.0)
        self.assert_almost_equal(df["energy"].values, -241.2364683)

        # FIXME
        #eos = robot.eos_fit()
        if self.has_nbformat():
            robot.write_notebook(nbpath=self.get_tmpname(text=True))

        robot.close()

    def test_sigres_robot(self):
        """Testing SIGRES robot."""
        filepaths = abidata.ref_files(
            "si_g0w0ppm_nband10_SIGRES.nc",
            "si_g0w0ppm_nband20_SIGRES.nc",
            "si_g0w0ppm_nband30_SIGRES.nc",
        )
        assert SigresRobot.class_handles_filename(filepaths[0])
        assert len(filepaths) == 3

        with SigresRobot.from_files(filepaths) as robot:
            assert robot.start is None
            start = robot.trim_paths(start=None)
            assert robot.start == start
            for p, _ in robot:
                assert p == os.path.relpath(p, start=start)

            assert robot.EXT == "SIGRES"
            repr(robot)
            str(robot)
            df_sk = robot.merge_dataframes_sk(spin=0, kpoint=[0, 0, 0])
            qpdata = robot.get_qpgaps_dataframe()
            if self.has_seaborn():
                robot.plot_conv_qpgap(x_vars="sigma_nband")

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))

            robot.pop_label(os.path.relpath(filepaths[0], start=start))
            assert len(robot) == 2
            robot.pop_label("foobar")

    def test_mdf_robot(self):
        """Testing MDF robot."""
        robot = MdfRobot.from_dir(os.path.join(abidata.dirpath, "refs", "si_bse_kpoints"))
        assert len(robot) == 3

        robot = MdfRobot()
        robot.scan_dir(os.path.join(abidata.dirpath, "refs", "si_bse_kpoints"))
        assert len(robot) == 3
        repr(robot)
        str(robot)

        df = robot.get_dataframe(with_geo=True)
        assert df is not None

        plotter = robot.get_multimdf_plotter()
        if self.has_matplotlib():
            plotter.plot()
            #robot.plot_conv_mdf(self, hue, mdf_type="exc_mdf", **kwargs):

        if self.has_nbformat():
            robot.write_notebook(nbpath=self.get_tmpname(text=True))

        robot.close()

    #def test_ddb_robot(self):
    #    """Testing DDB robots."""
    #    filepaths = abidata.ref_files(
    #    )
    #assert not DdbRobot.class_handles_filename(filepaths[0])
    #    assert len(filepaths) == 3
    #    repr(robot)
    #     str(robot)
    #     df = robot.get_dataframe_at_qpoint(self, qpoint=None, asr=2, chneut=1, dipdip=1, with_geo=True, **kwargs):

    #     phbands_plotter = robot.get_phbands_plotter()
    #     phdos_plotter = robot.get_phbands_plotter()
    #     if self.has_matplotlib():
    #       phbands_plotter.gridplot(show=False)
    #       phdos_plotter.gridplot(show=False)

    #     with DdbRobot.from_files(filepaths) as robot:
    #         assert robot.EXT == "DDB"
    #         robot.get_qpoints_union()
    #         df = robot.get_dataframe_at_qpoint(qpoint=None, asr=2, chneut=1, dipdip=1)
    #         if self.has_seaborn():
    #             robot.plot_conv_qpgap(x_vars="sigma_nband", show=False)

    #if self.has_nbformat():
    #    robot.write_notebook(nbpath=self.get_tmpname(text=True))

    #robot.close()
