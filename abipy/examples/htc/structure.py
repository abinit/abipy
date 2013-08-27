from __future__ import division

import numpy as np

import abipy.data as data
from abipy.core import Structure
from abipy.htc.input import AbiInput

structure = Structure.from_file(data.cif_file("si.cif"))

# From 98% to 102% of the initial volume.
v0 = structure.volume
volumes = v0 * np.arange(98, 104, 2) / 100.

inp = AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=len(volumes))

#scaler = LatticeScaler(structure, vol_ratios)
#new_structures = scaler.build_new_structures()
#
#for i, new_structure in enumerate(new_structures):
#    inp.set_structure(new_structure, dtset=i+1)


for idt, vol in enumerate(volumes):
    new_structure = structure.copy()
    new_structure.scale_lattice(vol)
    inp.set_structure(new_structure, dtset=idt+1)

#print(inp)

# Supercell
inp = AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=1)

new_structure = structure.copy()
new_structure.make_supercell(scaling_matrix=[1,2,3])

inp.set_structure(new_structure)
#print(inp)

# Frozen Phonon at q=0
inp = AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=3)

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

for dtset in range(1,inp.ndtset+1):
    new_structure = structure.copy()
    new_structure.translate_sites(indices=[0], vector=[0.01*dtset,0,0], frac_coords=True)
    inp.set_structure(new_structure, dtset=dtset)

print(inp)
