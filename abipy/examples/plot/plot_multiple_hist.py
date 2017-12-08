#!/usr/bin/env python
r"""
Multiple Structural relaxations
===============================

This example shows how to analyze the results of multiple 
structure relaxations with the HIST robot.
"""
from abipy import abilab 
import abipy.data as abidata

files = [
    abidata.ref_file("sic_relax_HIST.nc"),
    abidata.ref_file("sic_relax_HIST.nc"),
]

robot = abilab.HistRobot()
for i, f in enumerate(files):
    robot.add_file("file %s" % i, f)

# or alternatively use:
#
#   robot = HistRobot.from_files(files)

print(robot.get_dataframe())

for what in robot.what_list:
    robot.gridplot(what=what)
