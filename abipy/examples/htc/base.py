from __future__ import division

import abipy.data as data
from abipy.htc.input import AbiInput

# Create an ABINIT input file with 1 dataset. 
# Pseudos are located in the pseudo_dir directory.
inp = AbiInput(pseudos="14si.pspnc", pseudo_dir=data.pseudo_dir, ndtset=1)

# One can set the value of the variables directly with the syntax.
inp.ecut = 10.
inp.tolwfr = 1e-8

# It's possible to use strings but use them only for special cases such as:
inp.istwfk = '*1'       

# One can create a dictionary mapping keywords to values 
unit_cell = {
    "acell": 3*[10.217],       
    'rprim': [ [.0, .5, .5],
               [.5, .0, .5],
               [.5, .5, .0]],
    'ntypat': 1,
    'znucl': [14,],
    'natom': 2,
    'typat': [1, 1],
    'xred': [ [.0, .0, .0],
              [.25,.25,.25] ]
}

# and set the variables in the input file with the call:
inp.set_variables(**unit_cell)

# Alternatively, it's possible to create a dictionary on the fly with the syntax.
inp.set_variables(kptopt=1, 
                  ngkpt=[2, 2, 2], 
                  nshiftk=1, 
                  shiftk=[0.0, 0.0, 0.0]
                  )

# To print the input to stdout use:
print(inp)

# To write the input to a file, use:
# inp.write("run.abi")

# A slightly more complicated example: input file with two datasets
inp = AbiInput(pseudos="14si.pspnc", pseudo_dir=data.pseudo_dir, ndtset=2)

# Global variable common to all datasets.
inp.tolwfr = 1e-8

# To specify values for the different datasets, one can use the syntax
inp.ecut1 = 10
inp.ecut2 = 20

# or by passing the index of the dataset to set_variables via the dtset argument.
inp.set_variables(ngkpt=[2,2,2], tsmear=0.004, dtset=1)
inp.set_variables(kptopt=[4,4,4], tsmear=0.008, dtset=2)

print(inp)
