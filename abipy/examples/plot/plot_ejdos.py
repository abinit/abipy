#!/usr/bin/env python
"""
This example shows how plot the different contributions 
to the electronic joint density of states of Silicon
"""
import abipy.data as abidata
from abipy.abilab import abiopen

# Extract the bands computed with the SCF cycle on a Monkhorst-Pack mesh.
with abiopen(abidata.ref_file("si_scf_WFK.nc")) as wfk_file:
    ebands = wfk_file.ebands

# Select the valence and conduction bands to include in the JDOS
# Here we include valence bands from 0 to 3 and the first conduction band (4).
vrange = range(0,4)
crange = range(4,5)

# Plot data.
ebands.plot_ejdosvc(vrange, crange)

ebands.plot_ejdosvc(vrange, crange, cumulative=False)
