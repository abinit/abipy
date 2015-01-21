#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import numpy as np

import abipy.data as abidata
from abipy.abilab import Structure, StructureModifier, AbiInput

structure = Structure.from_file(data.cif_file("si.cif"))

modifier = StructureModifier(structure)

# From 98% to 102% of the initial volume.
new_structures = modifier.scale_lattice(vol_ratios=np.arange(98, 104, 2) / 100.)

inp = AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(new_structures))

for (idt, new_structure) in enumerate(new_structures):
    inp.set_structure(new_structure, dtset=idt+1)

print(inp)

# Supercell
inp = AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=1)

new_structure = modifier.make_supercell(scaling_matrix=[1,2,3])
inp.set_structure(new_structure)

print(inp)

structure = Structure.from_abivars({
    "acell" : [10, 10, 10],
    "rprim" : np.eye(3),
    "natom" : 2,
    "ntypat": 1,
    "typat" : [1, 1],
    "znucl" : 14,
    "xred"  : [[0,0,0], [1/2., 0, 0]],
    }
    )

# We need a new modifier since we have a new structure.
modifier = StructureModifier(structure)

frac_displ = np.zeros((len(structure), 3))
frac_displ[1,0] = 1.0
etas = [0.01, 0.02, 0.03]

new_structures = modifier.displace(frac_displ, etas)

# Frozen Phonon at q=0
inp = AbiInput(pseudos=abidata.pseudos("14si.pspnc"), ndtset=len(new_structures))

for (idt, new_structure) in enumerate(new_structures):
    inp.set_structure(new_structure, dtset=idt+1)

print(inp)
