#!/usr/bin/env python
r"""
Phonon Bands with/without ASR
=============================

This example shows how to plot a phonon band structure
with and without imposing the acoustic sum rule at the Gamma point.
"""
from abipy import abilab
import abipy.data as abidata

# Open DDB file
filepath = abidata.ref_file("mp-1009129-9x9x10q_ebecs_DDB")
ddb = abilab.abiopen(filepath)

# This method computes the phonon bands and DOS by calling anaddb
# with different values of asr and returns a PhononBandsPlotter object.
plotter = ddb.anacompare_asr(asr_list=(0, 2))
plotter.combiplot()

# Set nqsmall to 0 to disable DOS computation.
plotter = ddb.anacompare_asr(asr_list=(0, 2), nqsmall=0, ndivsm=10)
plotter.gridplot()
