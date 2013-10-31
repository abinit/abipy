#!/usr/bin/env python
"""
Band structure of silicon in a distorted geometry (frozen phonon at q=0)
"""
from __future__ import division, print_function

import sys
import os
import numpy as np
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import Tester, enable_logging


def make_scf_nscf_inputs(structure):
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=2)
    structure = inp.set_structure(structure)

    # Global variables
    global_vars = dict(ecut=6,
                       nband=8,
                       timopt=-1,
                       paral_kgb=1,
                       #nstep=4, # This is not enough to converge. Used to test the automatic restart.
                       nstep=10,
                    )

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    inp[1].set_variables(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    inp[2].set_kpath(ndivsm=6, kptbounds=kptbounds)
    inp[2].set_variables(tolwfr=1e-12)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()

    return scf_input, nscf_input


def bands_flow(workdir):
    manager = abilab.TaskManager.from_user_config()

    # build the structures
    base_structure = abilab.Structure.from_file(data.cif_file("si.cif"))
    modifier = abilab.StructureModifier(base_structure)

    etas = [-0.1, 0, +0.1]
    ph_displ = np.reshape(np.zeros(3*len(base_structure)), (-1,3))
    ph_displ[0,:] = [+1, 0, 0]
    ph_displ[1,:] = [-1, 0, 0]

    displaced_structures = modifier.displace(ph_displ, etas, frac_coords=False)

    flow = abilab.AbinitFlow(workdir, manager)

    for structure in displaced_structures:
        # Create the workflow for the band structure calculation.
        scf_input, nscf_input = make_scf_nscf_inputs(structure)
                                                                   
        work = abilab.BandStructureWorkflow(scf_input, nscf_input)
        flow.register_work(work)

    return flow.allocate()


@enable_logging
def main():
    tester = Tester()
    flow = bands_flow(tester.workdir)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
