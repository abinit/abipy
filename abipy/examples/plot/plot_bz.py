#!/usr/bin/env python
r"""
Brillouin zone
==============

This example shows how to display the Brillouin zone
with matplotlib or plotly.
"""

#%%
# Open the WKF file and extract the crystalline structure.

from abipy.abilab import abiopen
import abipy.data as abidata

wfk_file = abiopen(abidata.ref_file("si_scf_WFK.nc"))
structure = wfk_file.structure

#%%
# To visualize the BZ with matplotlib, use:

structure.plot_bz()

#%%
# For the plotly version, use:

structure.plotly_bz()

#%%
# Remember to close the wfk file with:
wfk_file.close()
