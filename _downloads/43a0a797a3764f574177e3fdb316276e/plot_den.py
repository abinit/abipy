#!/usr/bin/env python
r"""
Density
=======

This example shows how to analyze the electronic density
stored in the DEN.nc file.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open the DEN.nc file
ncfile = abiopen(abidata.ref_file("si_DEN.nc"))

# The DEN file has a `Density`, a `Structure` and an `ElectronBands` object
print(ncfile.structure)

# To plot the KS eigenvalues.
#ncfile.ebands.plot()

density = ncfile.density
print(density)

# To visualize the total charge wih vesta
#visu = density.visualize("vesta"); visu()

# To plot the density along the line connecting
# the first and the second in the structure:
density.plot_line(point1=0, point2=1)

# alternatively, one can define the line in terms of two points
# in fractional coordinates:
density.plot_line(point1=[0, 0, 0], point2=[2.25, 2.25, 2.25], num=300)

# To plot the density along the lines connect the firt atom in the structure
# and all the neighbors within a sphere of radius 3 Angstrom:
density.plot_line_neighbors(site_index=0, radius=3)
