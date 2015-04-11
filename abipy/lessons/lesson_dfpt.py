#!/usr/bin/env python
"""DFPT lesson: phonon band structure of AlAs with Born effective charges."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import numpy as np
import abipy.abilab as abilab
import abipy.data as abidata  


def make_scf_input(ecut=2, ngkpt=(4, 4, 4)):
    """
    This function constructs an `AbinitInput` for performing a 
    GS-SCF calculation in crystalline AlAs.

    Args:
        ecut: cutoff energy in Ha.
        ngkpt: 3 integers specifying the k-mesh for the electrons.

    Return:
        `AbinitInput` object 
    """
    # Initialize the AlAs structure from an internal database. Use the pseudos shipped with AbiPy.
    gs_inp = abilab.AbinitInput(structure=abidata.structure_from_ucell("AlAs"),
                                pseudos=abidata.pseudos("13al.981214.fhi", "33as.pspnc"))
    
    # Set the value of the Abinit variables needed for GS runs.
    gs_inp.set_vars(
        nband=4,             
        ecut=ecut,         
        ngkpt=ngkpt,
        nshiftk=4,
        shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
                0.0, 0.5, 0.0,
                0.5, 0.0, 0.0,
                0.5, 0.5, 0.5],
        ixc=1,
        nstep=25,
        diemac=9.0,
        tolvrs=1.0e-10,
    )
    
    gs_inp.set_mnemonics(True)
    return gs_inp

