#!/usr/bin/env python
r"""
Multiple Bethe-Salpeter calculations
====================================

This example shows how to analyze multiple MDF files.
"""
import abipy.data as abidata
from abipy import abilab

# Read data from multiple files.
mdf_paths = abidata.ref_files("si_444_MDF.nc", "si_666_MDF.nc", "si_888_MDF.nc")
robot = abilab.MdfRobot.from_files(mdf_paths)

# Build MultipleMdfPlotter
plotter = robot.get_multimdf_plotter()

# Plot the dielectric function with excitonic effects.
plotter.plot(mdf_type="exc", qview="avg",
             title="Real and Imaginary part (averaged over q-points)", tight_layout=True)

# Plot the dielectric function computed at the RPA level with KS energies.
# Show q-point dependence.
plotter.plot(mdf_type="rpa", qview="all",
             title="Real and Imaginary part for individual q-points", tight_layout=True)

robot.close()
