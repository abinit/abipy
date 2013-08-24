from __future__ import division

import numpy as np

import abipy.data as data
from abipy.core import Structure
from abipy.htc.input import AbiInput

inp = AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=1)

structure = Structure.from_file(data.cif_file("si.cif"))

# Supercell
new_structure = structure.copy()
new_structure.make_supercell(scaling_matrix=[1,2,3])

inp.set_structure(new_structure)
print(inp)

# From 94% to 106% of the initial volume.
v0 = structure.volume
volumes = v0 * np.arange(94, 108, 2) / 100.

inp = AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=len(volumes))

for idt, vol in enumerate(volumes):
    new_structure = structure.copy()
    new_structure.scale_lattice(vol)
    #print("vo", v0, new_structure.volume)
    if idt == 0:
        inp.set_structure(new_structure, dtset=idt)
    inp.set_structure(new_structure, dtset=idt+1)

print(inp)
