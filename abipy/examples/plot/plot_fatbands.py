#!/usr/bin/env python
#
# This example shows how to plot the so-called fatbands for electrons
import numpy as np
import collections

from abipy import *
import abipy.data as data

# Open the WKF file.
wfk_file = WFK_File(data.ref_file("si_nscf_WFK-etsf.nc"))

# Extract the band structure. 
ebands = wfk_file.ebands

# Define the mapping reduced_coordinates -> name of the k-point.
klabels = {
    (0.5,  0.0,  0.0) : "L",
    (0.0,  0.0,  0.0) : "$\Gamma$",
    (0.0,  0.5,  0.5) : "X",
}

# Set the width 
ones = np.ones(ebands.shape)
widths = collections.OrderedDict([
    ("Si-3s", ones),
    ("Si-3p", 2*ones),
])

for key, value in widths.items():
    ebands.set_width(key, value)

# Plot the fatbands
ebands.plot_fatbands()
