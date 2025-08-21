#!/usr/bin/env python
r"""
e-ph matrix along a path in the BZ
==================================

This example shows how to use the GPATH.nc file produced by ABINIT
to plot e-ph matrix elements along a path in the BZ.
"""
import os
import abipy.data as abidata

from abipy.eph.gpath import GpathFile

#%%
# Here we use one of the GPATH files shipped with abipy
# with the e-ph matrix g(k,q) with q along a path and k fixed.
# Replace the call to abidata.ref_file("") with the path to your GPATH.nc

gpath_q = GpathFile(abidata.ref_file("teph4zpr_9o_DS1_GPATH.nc"))
print(gpath_q)

#%%
# To plot the e-ph matrix elements as a function of q.
gpath_q.plot_g_qpath(band_range=None,
                     which_g="avg",
                     with_qexp=0,
                     scale=1,
                     gmax_mev=250,
                     ph_modes=None,
                     with_phbands=True,
                     with_ebands=False,
                     )

#%%
# For the meaning of the different arguments, see the docstring
print(gpath_q.plot_g_qpath.__doc__)

#%%
# Here we read another file
# with the e-ph matrix g(k,q) with k along a path and q fixed.

gpath_k = GpathFile(abidata.ref_file("teph4zpr_9o_DS2_GPATH.nc"))
print(gpath_k)

#%%
# To plot the e-ph matrix elements as a function of q.
gpath_k.plot_g_kpath(band_range=None,
                     which_g="avg",
                     scale=1,
                     gmax_mev=250,
                     ph_modes=None,
                     with_ebands=True,
                     )

#%%
# For the meaning of the different arguments, see the docstring
print(gpath_k.plot_g_kpath.__doc__)
