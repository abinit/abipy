#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import abipy.abilab as abilab 
import abipy.data as abidata

inp = abilab.AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)
structure = inp.set_structure_from_file(abidata.cif_file("si.cif"))

inp.set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0], dtset=1)

inp.set_kpath(ndivsm=5, dtset=2)
#inp.set_kpath(ndivsm=5, kptbounds=[[0,0,0], [0.5, 0.0, 0.0]], dtset=2)

# Global variables
inp.ecut = 10

# Dataset 1
inp.tolvrs1 = 1e-9

# Dataset 2
inp.tolwfr2 = 1e-18
inp.getden2 = -1

print(inp)
