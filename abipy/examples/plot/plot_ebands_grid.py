#!/usr/bin/env python
# This example shows how to plot several band structures on a grid
# We use two GSR files:
#       si_scf_GSR.n: energies on a homogeneous sampling of the BZ (can be used to compute DOS)
#       si_nscf_GSR.nc: energies on a k-path in the BZ (use to plot band dispersions)

from abipy.abilab import ebands_gridplot
from abipy.data import ref_file

# To plot a grid with two band structures:
eb_objects = [ref_file("si_scf_GSR.nc"), ref_file("si_nscf_GSR.nc")]
ebands_gridplot(eb_objects, titles=["si_scf_GSR.nc", "si_nscf_GSR.nc"])

# Use the optional argument `edos_objects` to plot a grid with band structures + DOS.
# The first subplot will get the band dispersion from eb_objects[0] and the dos from edos_objects[0]
# edos_kwargs is an optional dictionary with the arguments that will be passed to `get_dos` to compute the DOS.
eb_objects = 2 * [ref_file("si_nscf_GSR.nc")]
edos_objects = 2 * [ref_file("si_scf_GSR.nc")]

ebands_gridplot(eb_objects, edos_objects=edos_objects, edos_kwargs=None,
                titles=["Si ebands + DOS", "Same data"])
