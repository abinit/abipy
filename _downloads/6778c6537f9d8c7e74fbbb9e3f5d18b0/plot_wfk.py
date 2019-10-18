#!/usr/bin/env python
r"""
Wavefunction file
=================

This example shows how to analyze the wavefunctions
stored in the WFK.nc file.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open the DEN.nc file
ncfile = abiopen(abidata.ref_file("si_nscf_WFK.nc"))
print(ncfile)

# The DEN file has a `Structure` an `ElectronBands` object and wavefunctions.
#print(ncfile.structure)

# To plot the KS eigenvalues.
#ncfile.ebands.plot()

# Extract the wavefunction for the first spin, the first band and k=[0.5, 0, 0]
wave = ncfile.get_wave(spin=0, kpoint=[0.5, 0, 0], band=0)
print(wave)

# This is equivalent to
#wave = ncfile.get_wave(spin=0, kpoint=0, band=0)
# because [0.5, 0, 0] is the first k-point in the WFK file

# To visualize the total charge wih vesta
#visu = wave.visualize_ur2("vesta"); visu()

# To plot the wavefunction along the line connecting
# the first and the second in the structure:
#wave.plot_line(point1=0, point2=1)

#wave.plot_line(point1=0, point2=1, with_krphase=True)

# alternatively, one can define the line in terms of two points
# in fractional coordinates:
wave.plot_line(point1=[0, 0, 0], point2=[0, 4, 0], with_krphase=False, num=400)
wave.plot_line(point1=[0, 0, 0], point2=[0, 4, 0], with_krphase=True, num=400)

# To plot the wavefunction along the lines connect the firt atom in the structure
# and all the neighbors within a sphere of radius 3 Angstrom:
#wave.plot_line_neighbors(site_index=0, radius=3)

ncfile.close()
