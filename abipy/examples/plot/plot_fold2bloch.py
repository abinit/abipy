#!/usr/bin/env python
r"""
Unfolding with fold2bloch
=========================

This example shows how to plot the results produced by fold2bloch.
<http://www.abinit.org/doc/helpfiles/for-v8.0/tutorial/lesson_fold2Bloch.html>
"""
from __future__ import division, print_function

from abipy import abilab
import abipy.data as abidata
import numpy as np

with abilab.abiopen(abidata.ref_file("h6_FOLD2BLOCH.nc")) as ncfile:
    print(ncfile)
    # Plot folded bands
    ncfile.ebands.plot(title="Folded bands")

    # Plot unfolded bands along the path defined by kbounds.
    kbounds = [0, 1/2, 0, 0, 0, 0, 0, 0, 1/2]
    klabels = ["Y", r"$\Gamma$", "X"]
    # sphinx_gallery_thumbnail_number = 2
    ncfile.plot_unfolded(kbounds, klabels, title="Unfolded bands")
