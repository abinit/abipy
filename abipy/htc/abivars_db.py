"""Database with the names of the input variables used in Abinit and in other main programs."""
from __future__ import print_function, division

import os
import json

with open(os.path.join(os.path.dirname(__file__), "abinit_vars.json")) as fh:
    ABI_VARNAMES = json.load(fh)

# Add internal variables i.e. those variables that are not declared in dtset%
ABI_VARNAMES += [
    "acell",
    "xred",
    "rprim",
    "kptbounds",
    "ndivsm",
    "qpt",
]

# Unit names.
ABI_UNITS = [
    'au',
    'Angstr',
    'Angstrom',
    'Angstroms',
    'Bohr',
    'Bohrs',
    'eV',
    'Ha',
    'Hartree',
    'Hartrees',
    'K',
    'Ry',
    'Rydberg',
    'Rydbergs',
    'T',
    'Tesla',
]

# Operators.
ABI_OPS = ['sqrt', 'end', '*', '/']
