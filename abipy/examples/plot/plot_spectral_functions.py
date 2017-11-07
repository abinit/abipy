#!/usr/bin/env python
r"""
G0W0 spectral function
======================

This examples shows how to plot the G0W0 spectral functions A(w)
at the gamma point, for the first band). See lesson tgw2_4
"""
import abipy.data as abidata
from abipy.abilab import abiopen

# Open the file with the GW results
sigma_file = abiopen(abidata.ref_file("al_g0w0_sigmaw_SIGRES.nc"))

# Plot A(w) for the first spin, the gamma point, and bands in [0, 1, 2, 3]
sigma_file.plot_spectral_functions(spin=0, kpoint=(0, 0, 0), bands=range(0, 4))

sigma_file.close()
