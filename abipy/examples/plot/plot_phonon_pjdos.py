#!/usr/bin/env python
"""
This example shows how to plot the projected phonon DOS of AlAs.
See tutorial/lesson_rf2.html
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Read the Phonon DOS from the netcd file produced by anaddb with prtdos 2 and plot the data
# (alternatively one can use the shell and `abiopen.py OUT_PHDOS.nc -nb` to open the file in a jupyter notebook.
with abiopen(abidata.ref_file("trf2_5.out_PHDOS.nc")) as phdos_file:
    phdos_file.plot_pjdos_type(title="AlAs type-projected phonon DOS")
