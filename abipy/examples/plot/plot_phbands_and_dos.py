#!/usr/bin/env python
r"""
Phonon Band structures
======================

This example shows how to plot the phonon band structure of AlAs.
See tutorial/lesson_rf2.html
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open the PHBST file produced by anaddb and get the phonon bands.
with abiopen(abidata.ref_file("trf2_5.out_PHBST.nc")) as ncfile:
    phbands = ncfile.phbands

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2)
with abiopen(abidata.ref_file("trf2_5.out_PHDOS.nc")) as ncfile:
    phdos = ncfile.phdos

# plot phonon bands and DOS with LO-TO.
# `plot_phonons_lo_to.py` shows how to treat the LO-TO splitting.
phbands.plot(title="AlAs Phonon bands and DOS in eV")

# plot phonon bands with DOS.
phbands.plot_with_phdos(phdos, units="cm-1",
                        title="AlAs Phonon bands + DOS in cm-1")

# plot phonon DOS.
phdos.plot(units="cm-1", title="Phonon DOS and IDOS in cm-1")

# Plot the phonon band structure with different color for each line.
phbands.plot_colored_matched(units="cm-1",
                             title="AlAs with different color for each line.")
