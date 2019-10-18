#!/usr/bin/env python
r"""
Phonon Bands with LO-TO
=======================

This example shows how to plot the phonon band structure of AlAs
including the LO-TO splitting. See tutorial/lesson_rf2.html
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open PHBST file produced by anaddb and extract the phonon bands object.
# (alternatively one can use the shell and `abiopen.py OUT_PHBST.nc -nb`
# to open the file in a jupyter notebook.
with abiopen(abidata.ref_file("ZnSe_hex_886.out_PHBST.nc")) as ncfile:
    phbands = ncfile.phbands

# Phonon frequencies with non analytical contributions, if calculated, are saved
# in the anaddb.nc file produced by anaddb. The results should be fetched from there
# and added to the phonon bands.
phbands.read_non_anal_from_file(abidata.ref_file("ZnSe_hex_886.anaddb.nc"))

# Notice that all the directions starting from or arriving at gamma that are used
# in the path should explicitely calculated, even if the values are the same.

# Plot the phonon frequencies. Note that the labels for the q-points
# are found automatically by searching in an internal database.
phbands.plot(title="ZnSe with LO-TO splitting")
