#!/usr/bin/env python
r"""
Lobster COHP/COP/DOS
====================

This example shows how to analyze the output files
produced by Lobster code <<http://schmeling.ac.rwth-aachen.de/cohp/>

Use `abiview.py lobster DIRPATH` for a command line interface.
"""
import os
import abipy.data as abidata

from abipy.abilab import LobsterAnalyzer

dirpath = os.path.join(abidata.dirpath, "refs", "lobster_gaas")

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
