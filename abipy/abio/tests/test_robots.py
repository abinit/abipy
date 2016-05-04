# coding: utf-8
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

        assert Robot.for_ext("DDB") is DdbRobot
        with self.assertRaises(ValueError): Robot.for_ext("foo")

    def test_gsr_robot(self):
        """Testing GSR robot"""
        gsr_path = abidata.ref_file("si_scf_GSR.nc")
        robot = Robot.for_ext("GSR")()
        robot.add_file("gsr0", gsr_path)
        with self.assertRaises(ValueError): robot.add_file("gsr0", gsr_path)

        assert len(robot) == 1 and not robot.exceptions
        print(robot)
        print(repr(robot))
        robot.show_files()
        robot.close()

        plotter = robot.get_ebands_plotter()
        # FIXME
        #eos = robot.eos_fit()
        #frame = robot.get_dataframe()

    #def test_sigres_robot(self):
    #def test_mdf_robot(self):
    #def test_ddb_robot(self):


if __name__ == '__main__':
    import unittest
    unittest.main()
