#!/usr/bin/env python
r"""
Multiple phonon bands
=========================

This example shows how to plot several phonon band structures on a grid.

We use two files produced by anaddb:

    trf2_5.out_PHBST.nc: phonon frequencies on a q-path in the BZ (used to plot the band dispersion)
    trf2_5.out_PHDOS.nc: phonon DOS compute with anaddb.
"""
from abipy import abilab
import abipy.data as abidata

import os
paths = [
    #"mgb2_444k_0.01tsmear_DDB",
    #"mgb2_444k_0.02tsmear_DDB",
    #"mgb2_444k_0.04tsmear_DDB",
    "mgb2_888k_0.01tsmear_DDB",
    #"mgb2_888k_0.02tsmear_DDB",
    "mgb2_888k_0.04tsmear_DDB",
    "mgb2_121212k_0.01tsmear_DDB",
    #"mgb2_121212k_0.02tsmear_DDB",
    "mgb2_121212k_0.04tsmear_DDB",
]

paths = [os.path.join(abidata.dirpath, "refs", "mgb2_phonons_nkpt_tsmear", f) for f in paths]

robot = abilab.DdbRobot()
for i, path in enumerate(paths):
    robot.add_file(path, path)

robot.remap_labels(lambda ddb: "nkpt: %s, tsmear: %.2f" % (ddb.header["nkpt"], ddb.header["tsmear"]))

r = robot.anaget_phonon_plotters(nqsmall=2)

r.phbands_plotter.gridplot_with_hue("tsmear")

r.phbands_plotter.gridplot_with_hue("tsmear", with_dos=True)


#r.phbands_plotter.gridplot_with_hue("nkpt")
