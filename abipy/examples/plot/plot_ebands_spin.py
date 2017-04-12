#!/usr/bin/env python
"""
This example shows how to plot a band structure
using the eigenvalues stored in the GSR file produced by abinit at the end of the GS run.
"""
from abipy import abilab
import abipy.data as abidata

# Open the GSR file and extract the band structure.
# (alternatively one can use the shell and `abiopen.py OUT_GSR.nc -nb` to open the file in a jupyter notebook.
with abilab.abiopen(abidata.ref_file("ni_666k_GSR.nc")) as ncfile:
    ni_ebands_kmesh = ncfile.ebands

with abilab.abiopen(abidata.ref_file("ni_kpath_GSR.nc")) as ncfile:
    ni_ebands_kpath = ncfile.ebands

# Energy limits in eV for plots.
elims = [-10, 2]

ni_ebands_kpath.plot(ylims=elims, title="Ni band structure")

# Compute DOS with Gaussian method.
ni_edos = ni_ebands_kmesh.get_edos()

# TODO: Improve plots
ni_edos.plot(xlims=elims)
ni_edos.plot_dos_idos(xlims=elims)
ni_edos.plot_up_minus_down(xlims=elims)

ni_ebands_kpath.plot_with_edos(ni_edos, ylims=elims, title="Ni band structure + DOS")

#plotter = abilab.ElectronBandsPlotter()
#plotter.add_ebands("k-mesh", ni_ebands_kmesh)
#plotter.add_ebands("k-path", ni_ebands_kpath)
#plotter.combiboxplot(hue="spin")
