#!/usr/bin/env python
r"""
Projected phonon DOS
====================

This example shows how to plot several phonon DOS on a grid.
We phonon DOS file produced by anaddb:

    trf2_5.out_PHDOS.nc: phonon DOS compute with anaddb.

See also tutorial/lesson_rf2.html
"""

#%%
# We start by defining a list with the paths to the PHDOS.nc files
# In this case, for simplicity, we use the same file but we must
# use different labels when adding them to the plotter with the add_phdos method.

from abipy import abilab
import abipy.data as abidata

phdos_paths = [abidata.ref_file("trf2_5.out_PHDOS.nc")]

plotter = abilab.PhononDosPlotter()
plotter.add_phdos("AlAs", phdos_paths[0])
plotter.add_phdos("Same AlAs", phdos_paths[0])

#%%
# To visualize the results on a grid with matplotlib, use:
plotter.gridplot(units="cm-1", tight_layout=True, title="Phonon DOS in cm-1 on grid")

#%%
# For the plotly version:
plotter.gridplotly(units="cm-1", title="Phonon DOS in cm-1 on grid")
