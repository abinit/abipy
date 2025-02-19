#!/usr/bin/env python
r"""
Orbital magnetism
=================

This example shows how to plot ...
See ...
"""
from abipy.electrons.orbmag import OrbmagAnalyzer
#import abipy.data as abidata

#%%
# We start by defining a list with the paths to the PHDOS.nc files
# In this case, for simplicity, we use the same file but we must
# use different labels when adding them to the plotter with the add_phdos method.

import os
root = "/Users/giantomassi/git_repos/ABIPY_WORKS/JOE_ORBMAG"
filepaths = [os.path.join(root, s) for s in ["gso_DS12_ORBMAG.nc", "gso_DS22_ORBMAG.nc", "gso_DS32_ORBMAG.nc"]]

orban = OrbmagAnalyzer(filepaths)
print(orban)

choices = ['S','T','B','TB']
orban.report_eigvals(report_type="S")

#%%
# Plot the phonon frequencies. Note that the labels for the q-points
orban.plot_fatbands(os.path.join(root, "bandso_DS1_GSR.nc"))
