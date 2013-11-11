#!/usr/bin/env python
#
# This example shows how to plot the projected phonon DOS of AlAs.
# See tutorial/lesson_rf2.html
from abipy.phonons import PHDOS_File
import abipy.data as abidata

# Read the Phonon DOS from the netcd file produced by anaddb (prtdos 2)
phdos_file = PHDOS_File(abidata.ref_file("trf2_5.out_PHDOS.nc"))

# Plot data.
phdos_file.plot_pjdos_type(title="AlAs type-projected phonon DOS")
