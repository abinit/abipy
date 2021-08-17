#!/usr/bin/env python
r"""
ElectronDosPlotter
==================

This example shows how to plot several Electron Dos on a grid
using the results stored in the ni_666k_GSR.nc file.
"""

#%%
# We start by defining a list with the paths to the GSR.nc files
# In this case, for simplicity, we use the same file but we must
# use different labels when adding them to the plotter with the add_edos method.

from abipy import abilab
import abipy.data as abidata

edos_paths = [abidata.ref_file("ni_666k_GSR.nc")]

plotter = abilab.ElectronDosPlotter()
plotter.add_edos("Ni", edos_paths[0])
plotter.add_edos("Same Ni", edos_paths[0])

#%%
# To visualize the results on a grid with matplotlib, use:
plotter.gridplot(title="Electron DOS on grid")

#%%
# Plotly version:
plotter.gridplotly(title="Electron DOS on grid")

