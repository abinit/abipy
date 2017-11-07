#!/usr/bin/env python
r"""
LM-fatbands
===========

This example shows how to plot the LM-projected fatbands.
The FATBANDS file must have benn produced with prtdos 3 and prtdosm 1.
"""
import abipy.abilab as abilab
import abipy.data as abidata

fbnc_kpath = abilab.abiopen(abidata.ref_file("ni_kpath_FATBANDS.nc"))
print(fbnc_kpath)

# NC files have contributions up to L = 4 (g channel)
# but here we are intererested in s,p,d terms only so we use the optional argument lmax
lmax = 2

# we are not interested in a small energy window around the Fermi level.
elims = [-1.5, +1]

# and a subset of bands (remember that in python we start to count from 0)
blist = list(range(4, 10))

# Plot fatbands with LM character up to lmax.
# The grid contains (lmax + 1) columns, each column has (2l + 1) subplots
# corresponding to the LM character for M in [-l, -l-1, ... 0, 1, +l].
fbnc_kpath.plot_fatbands_mview(iatom=0, fact=1.5, lmax=lmax,
                               ylims=elims, blist=list(range(4, 10)),
                               title="LM fatbands for atom index 0")

fbnc_kpath.close()
