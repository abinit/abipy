#!/usr/bin/env python
r"""
Band structure plot
===================

This example shows how to analyze the wannier90 results using
the ABIWAN.nc netcdf file produced by Abinit when calling wannier90 in library mode.
Use `abiopen FILE.wout` for a command line interface and
the `--expose` option to generate matplotlib figures automatically.
"""
import os
import abipy.data as abidata

from abipy.abilab import abiopen

# Open the ABIWAN file
filepath = "/Users/gmatteo/git_repos/abinit_wannier/_build/tests/W90_REFS/tutoplugs_tw90_4/tw90_4o_DS3_ABIWAN.nc"
abiwan = abiopen(filepath)
print(abiwan)

# Plot the matrix elements of the KS Hamiltonian in real space in the Wannier Gauge.
abiwan.hwan.plot(title="Matrix elements in real space")

# Interpolate bands with Wannier using an automatically selected k-path
# with 3 points for the smallest segment.
ebands_wan = abiwan.interpolate_ebands(line_density=3)
ebands_wan.plot(title="Wannier-interpolated")

# To compare the interpolated bands with ab-initio results,
# pass a file with the ab-initio bands to get_plotter_from_ebands
# to build an ElectronBandsPlotter.
#plotter = abiwan.get_plotter_from_ebands("out_GSR.nc")
#plotter.combiplot()
