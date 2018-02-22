#!/usr/bin/env python
r"""
Projected phonon DOS
====================

This example shows how to plot the projected phonon DOS of AlAs.
See tutorial/lesson_rf2.html
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Read the Phonon DOS from the netcdf file produced by anaddb with prtdos 2
# (alternatively one can use the shell and `abiopen.py OUT_PHDOS.nc -nb`
# to open the file in a jupyter notebook.
with abiopen(abidata.ref_file("trf2_5.out_PHDOS.nc")) as phdos_file:

    # Plot PJDOS.
    phdos_file.plot_pjdos_type(units="cm-1", title="AlAs type-projected phonon DOS")

    # To have the projection along the cartesian directions (summed over atomic types)
    phdos_file.plot_pjdos_cartdirs_type(units="Thz", stacked=True,
            title="Type-projected ph-DOS decomposed along the three Cartesian directions.")

    # To plot the PJDOS for all the inequivalent sites.
    phdos_file.plot_pjdos_cartdirs_site(view="inequivalent", stacked=True)
