#!/usr/bin/env python
from __future__ import division

import abipy.data as data
from abipy.abilab import Structure

# One can export the crystalline structure reported in the Netcdf file
# produced by abinit with a few lines of code. 
# See also abipy/scripts/abistruct.py for a handy command line interface.

# Read the structure from a netcdf file
filepath = data.ref_file("si_nscf_GSR.nc")
structure = abilab.Structure.from_file(filepath)

# Call convert to get the string representation in the new format.
for format in ["cif", "POSCAR", "cssr", "json"]:
    print((" Abinit --> %s " % format).center(80, "*"))
    s = structure.convert(format=options.format)
    print(s)
