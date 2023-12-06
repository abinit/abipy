#!/usr/bin/env python
r"""
Plot the LRUJ results
=====================

This example shows how to parse the output file produced by lruj and plot the results


See also https://docs.abinit.org/tutorial/lruj/#43-execution-of-the-lruj-post-processinng-utility
"""
from abipy.electrons.lruj import LrujAnalyzer, LrujResults
import abipy.data as abidata
import os

# Initialize LrujResults from the main output file of lruj
outfile = abidata.ref_file("lruj_data/lruj.out")
lr = LrujResults.from_file(outfile)

#%%
# Plot the fits.

lr.plot(degrees="all", insetdegree=4, ptcolor0='blue', ptitle="Hello World", fontsize=9)

#filepaths = [
#    "tlruj_2.o_DS1_LRUJ.nc",
#    "tlruj_2.o_DS2_LRUJ.nc",
#    "tlruj_2.o_DS3_LRUJ.nc",
#    "tlruj_2.o_DS4_LRUJ.nc",
#]
#
#root = "tutorial_tlruj_1-tlruj_2-tlruj_3/"
#filepaths = [os.path.join(root, p) for p in filepaths]
#
#lruj = LrujAnalyzer(verbose=1)
#lruj.add_ncpaths("foo", filepaths)
#lruj.add_ncpaths("bar", filepaths)

#print(lruj)
#lruj.plot()
