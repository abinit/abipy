#!/usr/bin/env python
#
# This example shows how to compute the gaussian DOS from
# the eigenvalues stored in the WFK file.
import matplotlib.pyplot as plt

import abipy.data as data
from abipy.abilab import abiopen

# Open the wavefunction file computed with a homogeneous sampling of the BZ 
# and extract the band structure on the k-mesh.

gs_wfk = abiopen(data.ref_file("si_scf_WFK-etsf.nc"))

gs_ebands = gs_wfk.ebands

# Compute the DOS with the Gaussian method.
width = 0.1
step  = 0.01

edos = gs_ebands.get_edos(method="gaussian", step=step, width=width)

# Plot DOS and IDOS
edos.plot(title="Total IDOS and DOS of Silicon")
