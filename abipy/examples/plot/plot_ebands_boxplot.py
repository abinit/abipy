#!/usr/bin/env python
r"""
Electron bands box plot
=======================

This example shows how to use seaborn or plotly to draw a box plot
showing the distributions of eigenvalues with respect to the band index.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open the file with energies computed with a homogeneous sampling
# of the BZ and extract the band structure.
with abiopen(abidata.ref_file("si_scf_GSR.nc")) as gsr:
    ebands = gsr.ebands

#%%
# Use `swarm=True` to show the datapoints on top of the boxes
ebands.boxplot(swarm=True)

#%%
# For the plotly version, use:
ebands.boxplotly()
