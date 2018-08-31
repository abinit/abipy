#!/usr/bin/env python
r"""
Wannier90 wout file
===================

This example shows how to plot the convergence
of the wannierization cycle using the .wout file produced by wannier90.
Use `abiopen FILE.wout` for a command line interface and
the `--expose` option to generate matplotlib figures automatically.
"""
import os
import abipy.data as abidata

from abipy.abilab import abiopen

# Open the wout file
filepath = os.path.join(abidata.dirpath, "refs", "wannier90", "example01_gaas.wout")
wout = abiopen(filepath)
print(wout)

# Plot the convergence of the Wannierise cycle.
wout.plot(title="Wannierise cycle")

# Plot the convergence of the Wannier centers and spread
# as function of iteration number
wout.plot_centers_spread(title="MLWF centeres and spread")
