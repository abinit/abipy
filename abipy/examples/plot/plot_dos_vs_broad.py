# This example shows how to compute and plot several
# gaussian DOS by varying the broadening parameters.
from abipy import *
import matplotlib.pyplot as plt

# Open the wavefunction file computed with a homogeneous sampling of the BZ 
# and extract the band structure on the k-mesh.
gs_filename = get_ncfile("si_WFK-etsf.nc")

gs_wfk = WFK_File(gs_filename)

gs_bands = gs_wfk.get_bands()

# Compute the DOS with the Gaussian method.
# Plot bands and DOS.
widths = [0.1, 0.2, 0.3]
step = 0.1

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

for width in widths:
   dos = gs_bands.get_dos(method="gaussian", step=step, width=width)
   dos.plot_ax(ax, label="$\sigma = %s$ [eV]" % width)

ax.grid(True)
ax.set_xlabel("Energy [eV]")
ax.set_ylabel("DOS")
ax.legend()

plt.show()
