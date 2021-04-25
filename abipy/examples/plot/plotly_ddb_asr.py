#!/usr/bin/env python
r"""
Phonon bands with/without the ASR (Plotly version)
==================================================

This example shows how to plot a phonon band structure
with and without enforcing the acoustic sum rule (ASR).

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
# The ddb.anacompare_asr method computes the phonon bands and DOS by calling anaddb
# with different values of asr and returns a PhononBandsPlotter object.
# To make the computation faster we used the **advanced** options dipdip -1.

plotter = ddb.anacompare_asr(asr_list=(0, 2), dipdip=-1)
plotter.combiplotly()

#%%
# To disable the DOS computation, set nqsmall to 0:

plotter = ddb.anacompare_asr(asr_list=(0, 2), nqsmall=0, ndivsm=10, dipdip=-1)
plotter.gridplotly()

#%%
# Remember to close the file with:

ddb.close()

# sphinx_gallery_thumbnail_path = '_static/plotly_logo.png'
