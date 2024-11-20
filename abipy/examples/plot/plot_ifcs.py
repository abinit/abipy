#!/usr/bin/env python
r"""
Interatomic Force Constants
============================

This example shows how to plot the longitudinal part
of the Interatomic Force Constants in real-space
starting from a DDB file.

See e.g. <https://journals.aps.org/prb/pdf/10.1103/PhysRevB.50.13035>

.. important::

    Note that a **manager.yml** configuration file and an abinit installation are required
    to run this script as AbiPy needs to invoke anaddb to compute the IFCs from the DDB file.
"""

#%%
# Open the DDB file for AlAs

import os
import abipy.data as abidata

from abipy import abilab

filepath = os.path.join(abidata.dirpath, "refs", "alas_phonons", "trf2_3_DDB")
ddb = abilab.abiopen(filepath)

#%%
# Call anaddb to compute the Interatomic Force Constants
# Default args are: asr=2, chneut=1, dipdip=1.

ifc = ddb.anaget_ifc()

#%%
# Plot the total longitudinal IFCs in local coordinates,

ifc.plot_longitudinal_ifc(title="Total Longitudinal IFCs")

#%%
# Plot the short range longitudinal IFCs in local coordinates,

ifc.plot_longitudinal_ifc_short_range(title="Longitudinal SR-IFCs")

#%%
# Plot the Ewald part of the IFCs in local coordinates,

ifc.plot_longitudinal_ifc_ewald(title="Longitudinal LR-IFCs")
