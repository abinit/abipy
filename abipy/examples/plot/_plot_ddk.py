#!/usr/bin/env python
r"""
DDK.nc
======

This example shows how to analyze DDK.nc files 
containing the velocity matrix elements.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Here we use the DDK.nc files shipped with abipy.
# Replace ddk_paths with the paths to your files.
ddk_paths = abidata.ref_files(
    "gaas_444_dir1_DDK.nc",
    "gaas_444_dir2_DDK.nc",
    "gaas_444_dir3_DDK.nc",
)

from abipy.electrons.ddk import DdksAnalyzer

with DdksAnalyzer(ddk_paths) as van:
    print(van)
    #ddk_doses = van.get_doses(method="gaussian", step=0.1, width=0.2)
    #ddk_doses.plot()
    #ddk_doses.plot_with_ebands(ebands_kpath)
