#!/usr/bin/env python
r"""
X-ray diffraction pattern
=========================

This example shows how to plot the X-ray diffraction pattern with pymatgen
"""
from abipy import abilab
import abipy.data as abidata

# Extract the structure (supports multiple formats e.g netcdf, Abinit input, cif, POSCAR)
# Also available via the `abistruct.py xrd FILE` command line interface.
structure = abilab.Structure.from_file(abidata.ref_file("si_scf_WFK.nc"))

#%%
# Plot the XRD spectrum.
structure.plot_xrd()
