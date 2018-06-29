#!/usr/bin/env python
r"""
Lobster COHPCAR
===============

This example shows how to analyze the COHPCAR file
produced by Lobster code <http://schmeling.ac.rwth-aachen.de/cohp/index.php?menuID=1>
"""
from abipy.abilab import abiopen
import abipy.data as abidata

#filename = abidata.ref_file("si_nscf_GSR.nc")
import os
dirpath = "/Users/gmatteo/git_repos/abipy/abipy/test_files"
filename = os.path.join(dirpath, "GaAs_COHPCAR.lobster.gz")

# Open the COHPCAR.lobster file (same API for COOPCAR.lobster)
cohp = abiopen(filename)
print(cohp)

# Plot COHP.
cohp.plot(title="GaAs COHP")

# Plot integrate COHP.
cohp.plot(what="i", title="GaAs integrated COHP")

# Plot COHP total overlap for all sites listed in `from_site_index`
cohp.plot_site_pairs_total(from_site_index=[0], title="COHP total overlap for site index 0")

# Plot partial crystal orbital projections (DOS or IDOS) for all sites listed in `from_site_index`
cohp.plot_site_pairs_partial(from_site_index=[0], title="COHP with orbital projections from site index 0")
