#!/usr/bin/env python
r"""
MgB2 fatbands
=============

This example shows how to plot the L-projected fatbands of MgB2
using the FATBANDS.nc files produced by abinit with prtdos 3.
See also PhysRevLett.86.4656
"""
import abipy.abilab as abilab
import abipy.data as abidata

# Open the file (alternatively one can use the shell and `abiopen.py FILE -nb`
# to open the file in a jupyter notebook
# This file has been produced on a k-path so it's not suitable for DOS calculations.
fbnc_kpath = abilab.abiopen(abidata.ref_file("mgb2_kpath_FATBANDS.nc"))

# Print file info (dimensions, variables ...)
# Note that prtdos = 3, so LM decomposition is not available.
print(fbnc_kpath)


# Plot the k-points belonging to the path.
fbnc_kpath.ebands.kpoints.plot()

# NC files have contributions up to L=4 (g channel)
# but here we are intererested in s,p,d terms only so
# we use the optional argument lmax
lmax = 2

# Plot the electronic fatbands grouped by atomic type.
fbnc_kpath.plot_fatbands_typeview(lmax=lmax, tight_layout=True)

# Plot the electronic fatbands grouped by L.
fbnc_kpath.plot_fatbands_lview(lmax=lmax, tight_layout=True)

# Now we read another FATBANDS file produced on 18x18x18 k-mesh
fbnc_kmesh = abilab.abiopen(abidata.ref_file("mgb2_kmesh181818_FATBANDS.nc"))

print(fbnc_kmesh)
#fbnc_kmesh.ebands.kpoints.plot()

# Plot the L-PJDOS grouped by atomic type.
fbnc_kmesh.plot_pjdos_typeview(lmax=lmax, tight_layout=True)

# Plot the L-PJDOS grouped by L.
fbnc_kmesh.plot_pjdos_lview(lmax=lmax, tight_layout=True)

# Now we use the two netcdf files to produce plots with fatbands + PJDOSEs.
# The data for the DOS is taken from pjdosfile.
# sphinx_gallery_thumbnail_number = 6
fbnc_kpath.plot_fatbands_with_pjdos(pjdosfile=fbnc_kmesh, lmax=lmax,
                                    view="type", tight_layout=True)

# fatbands + PJDOS grouped by L
fbnc_kpath.plot_fatbands_with_pjdos(pjdosfile=fbnc_kmesh, lmax=lmax,
                                    view="lview", tight_layout=True)

fbnc_kpath.close()
fbnc_kmesh.close()
