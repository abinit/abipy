# coding: utf-8
"""Test for Robots"""
from __future__ import unicode_literals, division, print_function

import sys
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.core.testing import AbipyTest
from abipy.abio.robots import *


class RobotTest(AbipyTest):

    def test_base_robot_class(self):
        """Testing base robot class"""
        # With context.
        with Robot() as robot:
            print(robot)
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
        with self.assertRaises(ValueError):
            robot.add_file("gsr0", gsr_path)

        assert robot.EXT == "GSR"
        assert len(robot) == 1 and not robot.exceptions
        print(robot)
        print(repr(robot))
        robot.show_files()

        plotter = robot.get_ebands_plotter()
        if self.has_matplotlib():
            plotter.gridplot()

        df = robot.get_dataframe()
        assert df is not None
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
        assert len(filepaths) == 3

        with SigresRobot.from_files(filepaths) as robot:
            assert robot.EXT == "SIGRES"
            df_sk = robot.merge_dataframes_sk(spin=0, kpoint=[0, 0, 0])
            qpdata = robot.get_qpgaps_dataframe()
            if self.has_seaborn():
                robot.plot_conv_qpgap(x_vars="sigma_nband")

    def test_mdf_robot(self):
        """Testing MDF robot."""
        #robot = MdfRobot.from_files(abidata.ref_files("si_444_MDF.nc", "si_666_MDF.nc", "si_888_MDF.nc"))
        robot = MdfRobot.from_dir(os.path.join(abidata.dirpath, "refs", "si_bse_kpoints"))
        assert len(robot) == 3

        df = robot.get_dataframe(with_geo=True)
        assert df is not None

        plotter = robot.get_multimdf_plotter()
        if self.has_matplotlib():
            plotter.plot()
            #robot.plot_conv_mdf(self, hue, mdf_type="exc_mdf", **kwargs):

    #def test_ddb_robot(self):
    #    """Testing DDB robots."""
    #    filepaths = abidata.ref_files(
    #        "si_g0w0ppm_nband10_SIGRES.nc",
    #        "si_g0w0ppm_nband20_SIGRES.nc",
    #        "si_g0w0ppm_nband30_SIGRES.nc",
    #    )
    #    assert len(filepaths) == 3

    #    plotter = robot.get_phbands_plotter()
    #    if self.has_matplotlib():
    #       plotter.gridplot(show=False)

    #    with DdbRobot.from_files(filepaths) as robot:
    #        assert robot.EXT == "DDB"
    #        robot.get_qpoints_union()
    #        df = robot.get_dataframe_at_qpoint(qpoint=None, asr=2, chneut=1, dipdip=1)
    #        if self.has_seaborn():
    #            robot.plot_conv_qpgap(x_vars="sigma_nband", show=False)
