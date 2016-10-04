#!/usr/bin/env python
"""
This example shows how to plot the projected phonon DOS of AlAs.
See tutorial/lesson_rf2.html
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2) and plot data
with abiopen(abidata.ref_file("trf2_5.out_PHDOS.nc")) as phdos_file:
    phdos_file.plot_pjdos_type(title="AlAs type-projected phonon DOS")
