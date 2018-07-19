#!/usr/bin/env python
r"""
Band structure interpolation with Wannier functions
===================================================

This example shows how to analyze the wannier90 results using
the ABIWAN.nc netcdf file produced by Abinit when calling wannier90 in library mode.
Use `abiopen FILE.wout` for a command line interface and
the `--expose` option to generate matplotlib figures automatically.
"""
import os
import abipy.data as abidata

from abipy.abilab import abiopen

# Open the ABIWAN file
filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "tutoplugs_tw90_4", "tw90_4o_DS3_ABIWAN.nc")
abiwan = abiopen(filepath)
print(abiwan)

# Plot the matrix elements of the KS Hamiltonian in real space in the Wannier Gauge.
abiwan.hwan.plot(title="Matrix elements in real space")

# Interpolate bands with Wannier functions using an automatically selected k-path
# with 5 points for the smallest segment.
ebands_wan_kpath = abiwan.interpolate_ebands(line_density=5)
ebands_wan_kpath.plot(title="Wannier-interpolated")

# Interpolate bands in the IBZ defined by ngkpt
ebands_wan_kmesh = abiwan.interpolate_ebands(ngkpt=[8, 8, 8])
edos = ebands_wan_kmesh.get_edos()
ebands_wan_kpath.plot_with_edos(edos, title="Wannier-interpolated bands and DOS")

# To compare the interpolated bands with ab-initio results,
# pass a file with the ab-initio bands to the get_plotter_from_ebands method
# that will return an ElectronBandsPlotter object.

#plotter = abiwan.get_plotter_from_ebands("out_GSR.nc")
#plotter.combiplot()
