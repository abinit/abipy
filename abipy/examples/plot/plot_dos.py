# This example shows how to compute the gaussian DOS from
# the eigenvalues stored in the WFK file.
from abipy import *
import matplotlib.pyplot as plt

# Open the wavefunction file computed with a homogeneous sampling of the BZ 
# and extract the band structure on the k-mesh.
gs_filename = get_ncfile("si_WFK-etsf.nc")

gs_wfk = WFK_File(gs_filename)

gs_bands = gs_wfk.get_bands()

# Compute the DOS with the Gaussian method.
width = 0.1
step  = 0.01

dos = gs_bands.get_dos(method="gaussian", step=step, width=width)

# Plot DOS and IDOS
dos.plot(title="Total IDOS and DOS of Silicon")
