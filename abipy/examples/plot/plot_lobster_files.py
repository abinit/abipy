#!/usr/bin/env python
r"""
Lobster COHP/COP/DOS
====================

This example shows how to analyze the output files
produced by Lobster code <http://schmeling.ac.rwth-aachen.de/cohp/index.php?menuID=1>
"""
from abipy.abilab import LobsterAnalyzer
import abipy.data as abidata

#filename = abidata.ref_file("si_nscf_GSR.nc")

import os
dirpath = "/Users/gmatteo/git_repos/abipy/abipy/test_files"
#filename = os.path.join(dirpath, "GaAs_COHPCAR.lobster.gz")

# Open the all the lobster files produced in directory dirpath
# with the (optional) prefix GaAs_
lobana = LobsterAnalyzer.from_dir(dirpath, prefix="GaAs_")

print(lobana)

# Plot COOP + COHP + DOS.
lobana.plot(title="COOP + COHP + DOS")

# Plot COHP for all sites in from_site_index and Lobster DOS.
lobana.plot_coxp_with_dos(from_site_index=[0, 1])

# Plot orbital projections.
lobana.plot_coxp_with_dos(from_site_index=[0], with_orbitals=True)

#lobana.plot_with_ebands(ebands="out_GSR.nc")
