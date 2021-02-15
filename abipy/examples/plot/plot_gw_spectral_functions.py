#!/usr/bin/env python
r"""
G0W0 spectral function
======================

This examples shows how to plot the G0W0 spectral functions A(w) at the gamma point.
See also lesson tgw2_4
"""
import abipy.data as abidata
from abipy.abilab import abiopen

# Open the file with the GW results
sigres = abiopen(abidata.ref_file("al_g0w0_sigmaw_SIGRES.nc"))

# Plot A(w) for the first spin, the gamma point, and all bands
sigres.plot_spectral_functions()

# Only bands in [0, 1, 2]
sigres.plot_spectral_functions(spin=0, kpoint=(0, 0, 0), include_bands=range(0, 3))

sigres.close()
