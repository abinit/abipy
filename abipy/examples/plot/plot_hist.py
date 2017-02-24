#!/usr/bin/env python
"""
This example shows how to analyze the results of a structure relaxation
using the HIST.nc file.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open the HIST file.
hist = abiopen(abidata.ref_file("sic_relax_HIST.nc"))

# Structure at the end of the structural relaxation.
print(hist.final_structure)

hist.plot(tight_layout=True)

hist.plot_energies(tight_layout=True)

hist.close()
