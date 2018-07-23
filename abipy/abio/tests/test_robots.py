# coding: utf-8
"""Test for Robots"""
from __future__ import unicode_literals, division, print_function

import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab

from abipy.core.testing import AbipyTest
from abipy.abio.robots import Robot


class RobotTest(AbipyTest):

    def test_base_robot_class(self):
        """Testing base robot class"""
        # This is for the abstract interface.
        class MyRobot(Robot):
            EXT = "FOOBAR.nc"
            def write_notebook(self, nbpath=None):
                raise NotImplementedError()
            def yield_figs(self, **kwargs):  # pragma: no cover
                raise NotImplementedError()

        # With context.
        with MyRobot() as robot:
            repr(robot); str(robot)
            assert robot.to_string(verbose=2)
            assert robot._repr_html_()
            assert len(robot) == 0 and not robot.exceptions
            robot.show_files()

        assert Robot.class_for_ext("DDB") is abilab.DdbRobot
        assert Robot.class_for_ext("OPTIC") is abilab.OpticRobot
        with self.assertRaises(ValueError):
            Robot.class_for_ext("foo")

        if self.has_nbformat():
            assert robot.get_baserobot_code_cells()
