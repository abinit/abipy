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

# Build HIST robot from paths, use user-selected labels.
robot = abilab.HistRobot.from_files(files, labels=["sic_relaxation", "same_hist"])

# Get dataframe with most important results.
print(robot.get_dataframe())

robot.gridplot()

# Use combiplot to compare the results
# (meaningless in this case as we are using the same HIST file).
robot.combiplot()
