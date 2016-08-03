#!/usr/bin/env python
"""
This example shows how to plot the phonon band structure of AlAs.
See tutorial/lesson_rf2.html
"""

# FIXME: LO-TO splitting and phonon displacements instead of eigenvectors.
from abipy.abilab import abiopen 
import abipy.data as abidata

# Open the PHBST file produced by anaddb and get the phonon bands.
with abiopen(abidata.ref_file("trf2_5.out_PHBST.nc")) as ncfile:
    phbands = ncfile.phbands

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2)
with abiopen(abidata.ref_file("trf2_5.out_PHDOS.nc")) as ncfile:
    phdos = ncfile.phdos

# plot phonon bands and DOS.
phbands.plot_with_phdos(phdos, title="AlAs Phonon bands and DOS")
