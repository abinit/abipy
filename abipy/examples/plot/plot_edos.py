#!/usr/bin/env python
r"""
Electron DOS
============

This example shows how to compute the DOS from
the eigenvalues stored in the WFK file with the gaussian method.
"""
import abipy.data as abidata
from abipy.abilab import abiopen

# Open the wavefunction file computed with a homogeneous sampling of the BZ
# and extract the band structure on the k-mesh.
with abiopen(abidata.ref_file("si_scf_WFK.nc")) as gs_wfk:
    gs_ebands = gs_wfk.ebands

# Compute the DOS with the Gaussian method (default)
edos = gs_ebands.get_edos(method="gaussian", step=0.01, width=0.1)

#%%
# To plot electron DOS and IDOS with matplotlib use:
edos.plot(title="Silicon DOS")

#%%
# For the plotly version use:
edos.plotly(title="Silicon DOS")

#%%
# To plot electron DOS and IDOS with matplotlib use:
edos.plot_dos_idos(title="DOS and Integrated DOS")

#%%
# For the plotly version use:
edos.plotly_dos_idos(title="DOS and Integrated DOS")
