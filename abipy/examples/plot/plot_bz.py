#!/usr/bin/env python
#
# This example shows how to display the Brillouin zone 
# with pymatgen and matplotlib.
from abipy.abilab import abiopen
import abipy.data as data

# Open the WKF file.
wfk_file = abiopen(data.ref_file("si_scf_WFK-etsf.nc"))

# Extract the crystalline structure.
structure = wfk_file.structure

# Visualize the BZ.
structure.show_bz()
