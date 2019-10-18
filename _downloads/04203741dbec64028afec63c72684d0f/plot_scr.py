#!/usr/bin/env python
r"""
Dielectric function with LFE
============================

This examples shows how to plot the macroscopic dielectric function
computed in the GW code (optdriver 3)
"""
import abipy.data as abidata
from abipy.abilab import abiopen

with abiopen(abidata.ref_file("sio2_SCR.nc")) as ncfile:
    # The SCR file contains a structure and electron bands in the IBZ.
    # We can thus use the ebands object to plot bands + DOS.
    print(ncfile)

    edos = ncfile.ebands.get_edos()
    ncfile.ebands.plot_with_edos(edos, title="KS energies used to compute the SCR file.")

    # sphinx_gallery_thumbnail_number = 2
    ncfile.plot_emacro(title="Macroscopic dielectric function of $SiO_2$ with local-field effects.")

    ncfile.plot_eelf(title="Electron Energy Loss Function of $SiO_2$")
