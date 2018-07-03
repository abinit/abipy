#!/usr/bin/env python
r"""
Lobster COHPCAR
===============

This example shows how to analyze the COHPCAR file
produced by Lobster code <http://schmeling.ac.rwth-aachen.de/cohp/>

Use `abiopen.py FILE` with --expose or --print for a command line interface
and --notebook to generate a jupyter notebook.
"""
import os
import abipy.data as abidata

from abipy.abilab import abiopen

dirpath = os.path.join(abidata.dirpath, "refs", "lobster_gaas")
filename = os.path.join(dirpath, "GaAs_COHPCAR.lobster.gz")

# Open the COHPCAR.lobster file (same API for COOPCAR.lobster)
cohp = abiopen(filename)
print(cohp)

# Plot COHP.
cohp.plot(title="GaAs COHP")

# Plot integrated COHP.
cohp.plot(what="i", title="GaAs integrated COHP")

# Plot total overlap for all sites listed in `from_site_index`
cohp.plot_site_pairs_total(from_site_index=[0], title="COHP total overlap for site index 0")

# Plot partial crystal orbital projections for all sites listed in `from_site_index`
cohp.plot_site_pairs_partial(from_site_index=[0], title="COHP with orbital projections from site index 0")
