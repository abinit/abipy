#!/usr/bin/env python
#
# This example shows how to generate an energy-dependent scissors operator
# by fitting the GW QPState corrections as function of the KS eigenvalues
# We then use the scissors operator to correct the KS band structure 
# computed on a high symmetry k-path. Finally, the LDA and the QPState band
# structure are plotted with matplotlib.
from abipy import *
import abipy.data as data

# Get the quasiparticle results from the SIGRES.nc database.
sigma_file = SIGRES_File(data.ref_file("tgw1_9o_DS4_SIGRES.nc"))
qplist_spin = sigma_file.qplist_spin

# Construct the scissors operator
domains = [[-10, 6.02], [6.3, 18]]
scissors = qplist_spin[0].build_scissors(domains, bounds=None)

# Read the KS band energies computed on the k-path
ks_bands = WFK_File(data.ref_file("si_nscf_WFK-etsf.nc")).ebands

# Read the KS band energies computed on the Monkhorst-Pack (MP) mesh
# and compute the DOS with the Gaussian method
ks_mpbands = WFK_File(data.ref_file("si_WFK-etsf.nc")).ebands
ks_dos = ks_mpbands.get_edos()

# Apply the scissors operator first on the KS band structure 
# along the k-path then on the energies computed with the MP mesh. 
qp_bands = ks_bands.apply_scissors(scissors)

qp_mpbands = ks_mpbands.apply_scissors(scissors)

# Compute the DOS with the modified QPState energies.
qp_dos = qp_mpbands.get_edos()

# Plot the LDA and the QPState band structure with matplotlib.
plotter = ElectronBandsPlotter()

plotter.add_ebands("LDA", ks_bands, dos=ks_dos)

plotter.add_ebands("LDA+scissors(e)", qp_bands, dos=qp_dos)

# Define the mapping reduced_coordinates -> name of the k-point.
klabels = {
    (0.5,  0.0,  0.0): "L",
    (0.0,  0.0,  0.0): "$\Gamma$",
    (0.0,  0.5,  0.5): "X",
}

plotter.plot(title="Silicon band structure", klabels=klabels)
