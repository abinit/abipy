#!/usr/bin/env python
#
# This example shows how to plot the macroscopic 
# dielectric function (MDF) computed in the Bethe-Salpeter code
from abipy import abiopen
import abipy.data as data

# Open the MDF file produced in the tutorial.
mdf_file = abiopen(data.ref_file("tbs_4o_DS2_MDF.nc"))

# Plot the imaginary part of the macroscopic 
# dielectric function (EXC, RPA, GWRPA) betwee 2 and 5 eV.
title = "Si absorption spectrum: EXC vs RPA"
mdf_file.plot_mdfs(title=title, xlim=(2, 5))

# Plot the 6 different components of the macroscopic dielectric tensor 
title = "Si macroscopic dielectric tensor"
tensor_exc = mdf_file.get_tensor("exc")
tensor_exc.plot(title=title)

