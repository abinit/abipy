#!/usr/bin/env python
r"""
Brillouin zone
==============

This example shows how to display the Brillouin zone
with pymatgen and matplotlib.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open the WKF file.
wfk_file = abiopen(abidata.ref_file("si_scf_WFK.nc"))

# Extract the crystalline structure.
structure = wfk_file.structure

# Visualize the BZ.
structure.plot_bz()

# Close the wfk file
wfk_file.close()
