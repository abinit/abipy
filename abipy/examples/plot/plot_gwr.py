#!/usr/bin/env python
r"""
Visualization of GWR results
============================

This example shows how to visualize the results produced by the GWR code.
"""
from abipy import abilab
import abipy.data as abidata

#%%
# Open the GWR.nc file
# Here we use one of the GSR files shipped with abipy.
# Replace filename with the path to your GSR file or your WFK file.

from abipy.electrons.gwr import GwrFile
gwr = GwrFile(abidata.ref_file("t01o_DS3_GWR.nc"))
print(gwr)

kpoint = 0
include_bands = "gap"

#%%
# Plot Sigma_nk(iw) along the imaginary axis for given k-point.
gwr.plot_sigma_imag_axis(kpoint=kpoint, include_bands=include_bands)

#%%
# Plot Sigma_nk(iw) along the real axis for given k-point.
gwr.plot_sigma_real_axis(kpoint=kpoint)

#%%
# Plot the spectral function A(w) along the real axis.
gwr.plot_spectral_functions()

#%%
# Plot qp-data versus KS energy
gwr.plot_qps_vs_e0()
