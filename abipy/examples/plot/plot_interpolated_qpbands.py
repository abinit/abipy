#!/usr/bin/env python
"""
This example shows how to interpolate the GW corrections and use the interpolated
values to correct the KS band structure computed on a high symmetry k-path and
the KS energies of a k-mesh. Finally, the KS and the GW results are plotted with matplotlib.
"""
import abipy.data as abidata
from abipy.abilab import abiopen, ElectronBandsPlotter

# Get quasiparticle results from the SIGRES.nc database.
sigres = abiopen(abidata.ref_file("si_g0w0ppm_nband30_SIGRES.nc"))

# Read the KS band energies computed on the k-path
ks_ebands_kpath = abiopen(abidata.ref_file("si_nscf_GSR.nc")).ebands

# Read the KS band energies computed on the Monkhorst-Pack (MP) mesh
# and compute the DOS with the Gaussian method
ks_ebands_kmesh = abiopen(abidata.ref_file("si_scf_GSR.nc")).ebands
ks_edos = ks_ebands_kmesh.get_edos()

# Interpolate QP corrections and apply them on top of the KS band structures.
# QP band energies are returned in r.qp_ebands_kpath and r.qp_ebands_kmesh.
r = sigres.interpolate(lpratio=5,
                       ks_ebands_kpath=ks_ebands_kpath,
                       ks_ebands_kmesh=ks_ebands_kmesh)
qp_edos = r.qp_ebands_kmesh.get_edos()

# Shortcut.
#r = sigres.interpolate(ks_ebands_kpath=abidata.ref_file("si_nscf_GSR.nc"),
#                       ks_ebands_kmesh=abidata.ref_file("si_scf_GSR.nc"))
#ks_edos = r.ks_ebands_kmesh.get_edos()
#qp_edos = r.qp_ebands_kmesh.get_edos()

# Plot the LDA and the QPState band structure with matplotlib.
plotter = ElectronBandsPlotter()

plotter.add_ebands("LDA", ks_ebands_kpath, dos=ks_edos)
plotter.add_ebands("GW (interpolated)", r.qp_ebands_kpath, dos=qp_edos)

# By default, the two band energies are shifted wrt to *their* fermi level.
# Use e=0 if you don't want to shift the eigenvalus
# so that it's possible to visualize the QP corrections.
plotter.combiplot(title="Silicon band structure")

plotter.gridplot(title="Silicon band structure")
