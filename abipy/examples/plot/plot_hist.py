#!/usr/bin/env python
r"""
Structural relaxation
=====================

This example shows how to analyze the results of a
structure relaxation using the HIST.nc file.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

#%%
# Open the HIST file.
# (alternatively one can use the shell and `abiopen.py OUT_HIST.nc -nb`
# to open the file in jupyter notebook.

hist = abiopen(abidata.ref_file("sic_relax_HIST.nc"))

#%%
# The structure at the end of the structural relaxation.

print(hist.final_structure)

#%%
# Plot the evolution of the lattice parameters, forces, etotal, etc.

hist.plot(tight_layout=True)

#%%
# To plot the same results with plotly and the "dark" template, use:

hist.plotly(template="plotly_dark")
# sphinx_gallery_thumbnail_number = 2

#%%
# Plot the total energy at the different relaxation steps.

hist.plot_energies()

#%%
# For the plotly version use:
hist.plotly_energies()


hist.close()
