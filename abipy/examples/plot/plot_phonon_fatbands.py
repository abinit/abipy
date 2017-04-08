#!/usr/bin/env python
"""
This example shows how to plot the phonon fatbands of AlAs.
See tutorial/lesson_rf2.html
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open the PHBST file produced by anaddb and get the phonon bands.
with abiopen(abidata.ref_file("trf2_5.out_PHBST.nc")) as ncfile:
    phbands = ncfile.phbands

# Mapping reduced coordinates -> labels
#qlabels = {
#    (0,0,0): "$\Gamma$",
#    (0.375, 0.375, 0.75): "K",
#    (0.5, 0.5, 1.0): "X",
#    (0.5, 0.5, 0.5): "L",
#    (0.5, 0.0, 0.5): "X",
#    (0.5, 0.25, 0.75): "W",
#}

# Plot the phonon band structure.
phbands.plot_fatbands(title="AlAs phonon fatbands without LO-TO splitting") #, qlabels=qlabels)

# Plot the phonon band structure + PJDOS
phdos_path = abidata.ref_file("trf2_5.out_PHDOS.nc")
phbands.plot_fatbands(units="Thz", phdos_file=phdos_path, 
                      title="AlAs phonon fatbands with PJDOS")
