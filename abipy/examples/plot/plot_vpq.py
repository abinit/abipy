#!/usr/bin/env python
r"""
e-ph matrix along a path in the BZ
==================================

This example shows how to use the VPQ.nc file produced by ABINIT
to plot e-ph matrix elements along a path in the BZ.
"""
import os
import abipy.data as abidata

from abipy.eph.vpq import VpqFile

filepath = "/Users/giantomassi/git_repos/abinit/_build/tests/POLARON/varpeq6/out_VPQ.nc"

#%%
# Here we use one of the GPATH files shipped with abipy
# with the e-ph matrix g(k,q) with q along a path and k fixed.
# Replace the call to abidata.ref_file("") with the path to your GPATH.nc

varpeq = VpqFile(abidata.ref_file(filepath))
print(varpeq)

#%%
# To plot the e-ph matrix elements as a function of q.
for polaron in varpeq.polaron_spin:
    df = polaron.get_final_results_df(with_params=True)
    polaron.plot_scf_cycle()
    #polaron.plot_ank_with_ebands(ebands_kpath, ebands_kmesh=None)
    #polaron.plot_bqnu_with_ddb("in_DDB", with_phdos=True)
    #polaron.plot_bqnu_with_phbands(phbands_qpath)


#%%
# For the meaning of the different arguments, see the docstring
#print(varpeq.plot_g_qpath.__doc__)
