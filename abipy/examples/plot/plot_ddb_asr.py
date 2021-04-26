#!/usr/bin/env python
r"""
Phonon bands with/without the ASR
=================================

This example shows how to plot a phonon band structure
with and without enforcing the acoustic sum rule (ASR).
Both matplotlib and plotly are supported.

.. important::

    Note that a **manager.yml** configuration file and an abinit installation are required
    to run this script as AbiPy needs to invoke anaddb to compute phonons from the DDB file.
"""

#%%
# Open the DDB file with:

from abipy import abilab
import abipy.data as abidata
filepath = abidata.ref_file("mp-1009129-9x9x10q_ebecs_DDB")
ddb = abilab.abiopen(filepath)

#%%
# The ``ddb.anacompare_asr`` method computes the phonon bands and the DOS by calling anaddb
# with different values of asr and returns a PhononBandsPlotter object:
# To make the computation faster, we use the **advanced** options dipdip -1.
# This option should produce results similar to dipdip 1 yet make sure to test
# the effect of this variable before using it in production.

plotter = ddb.anacompare_asr(asr_list=(0, 2), dipdip=-1)
print(plotter)

#%%
# To plot the bands on the same figure with matplotlib, use:

plotter.combiplot()

#%%
# For the plotly version, use:

plotter.combiplotly()

#%%
# To disable the DOS computation, set ``nqsmall` to 0:

plotter = ddb.anacompare_asr(asr_list=(0, 2), nqsmall=0, ndivsm=10, dipdip=-1)

#%%
# To plot the bands on different subplots with matplotlib, use:

plotter.gridplot()

#%%
# For the plotly version, use:
plotter.gridplotly()

#%%
# Finally, remember to close the file with:

ddb.close()
