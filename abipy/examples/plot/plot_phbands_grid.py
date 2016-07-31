#!/usr/bin/env python
"""
This example shows how to plot several phonon band structures on a grid.
We use two files produced by anaddb:
    trf2_5.out_PHBST.nc: phonon frequencies on a q-path in the BZ (used to plot the band dispersion)
    trf2_5.out_PHDOS.nc: phonon DOS compute with anaddb.

See also tutorial/lesson_rf2.html
"""
from abipy.abilab import phbands_gridplot
from abipy.data import ref_file

# To plot a grid with two band structures:
phb_objects = 2 * [ref_file("trf2_5.out_PHBST.nc")]
phbands_gridplot(phb_objects, titles=["AlAs", "Same AlAs"])

# To plot a grid with band structures + DOS, use the optional argument `phdos_objects`
# The first subplot will get the band dispersion from phb_objects[0] and the dos from phdos_objects[0]
# phdos_kwargs is an optional dictionary with the arguments that will be passed to `get_phdos` to compute the DOS.
phb_objects = 3 * [ref_file("trf2_5.out_PHBST.nc")]
phdos_objects = 3 * [ref_file("trf2_5.out_PHDOS.nc")]

phbands_gridplot(phb_objects, phdos_objects=phdos_objects, phdos_kwargs=None,
                 titles=["AlAs phbands + DOS", "Same data", "Same data"])
