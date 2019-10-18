#!/usr/bin/env python
r"""
Band structure plot
===================

This example shows how to plot a band structure
using the eigenvalues stored in the GSR file
produced at the end of the GS run.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Here we use one of the GSR files shipped with abipy.
# Replace filename with the path to your GSR file or your WFK file.
filename = abidata.ref_file("si_nscf_GSR.nc")

# Open the GSR file and extract the band structure.
# (alternatively one can use the shell and `abiopen.py OUT_GSR.nc -nb`
# to open the file in a jupyter notebook.
with abiopen(filename) as ncfile:
    ebands = ncfile.ebands

# Plot the band energies. Note that the labels for the k-points
# are found automatically in an internal database.
# Show fundamental and direct gaps.
#ebands.plot(with_gaps="fd", title="Silicon band structure")
ebands.plot(with_gaps=True, title="Silicon band structure")

# Plot the BZ and the k-point path.
ebands.kpoints.plot()
