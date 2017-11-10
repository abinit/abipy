#!/usr/bin/env python
r"""
GW and spectral function
========================

This examples shows how to plot the GW self-energy
and the spectral function of Al at the gamma point (first band)
See lesson tgw2_4
"""
import abipy.data as abidata
from abipy.abilab import abiopen

sigma_file = abiopen(abidata.ref_file("al_g0w0_sigmaw_SIGRES.nc"))

sigmaw = sigma_file.get_sigmaw(spin=0, kpoint=(0,0,0), band=0)
sigmaw.plot()
