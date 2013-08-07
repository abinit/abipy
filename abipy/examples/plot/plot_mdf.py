#!/usr/bin/env python
#
# This example shows how to plot the macroscopic 
# dielectric function (MDF) computed in the Bethe-Salpeter code
from abipy import *
import abipy.data as data

# Open the MDF file produced in the tutorial.
mdf_file = MDF_File(data.ref_file("tbs_4o_DS2_MDF.nc"))

# Plot the imaginary part of the macroscopic 
# dielectric function (EXC, RPA, GWRPA) betwee 2 and 5 eV.
title = "Si absorption spectrum: EXC vs RPA"
mdf_file.plot_mdfs(title=title, xlim=(2, 5))
