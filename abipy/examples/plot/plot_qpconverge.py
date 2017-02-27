#!/usr/bin/env python
"""
This example shows how to visualize the convergence of the
QP results stored in the SIGRES produced by the GW code (sigma run).
"""
from abipy.abilab import SigresPlotter
import abipy.data as abidata

# List of SIGRES files computed with different values of nband.
filenames = [
    "si_g0w0ppm_nband10_SIGRES.nc",
    "si_g0w0ppm_nband20_SIGRES.nc",
    "si_g0w0ppm_nband30_SIGRES.nc",
]

filepaths = [abidata.ref_file(fname) for fname in filenames]

# Create the plotter and add `filepaths` to the plotter.
plotter = SigresPlotter()
plotter.add_files(filepaths)

# Plot the convergence of the QP gaps.
#plotter.plot_qpgaps(title="QP gaps vs sigma_nband", hspan=0.05)

# Plot the convergence of the QP energies.
###plotter.plot_qpenes(title="QP energies vs sigma_nband", hspan=0.05)

plotter.plot_qps_vs_e0(tight_layout=True)
