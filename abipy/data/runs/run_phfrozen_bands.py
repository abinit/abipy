#!/usr/bin/env python
"""
Band structure of silicon in a distorted geometry (frozen phonon at q=0)
"""
from __future__ import division, print_function, unicode_literals

import sys
import os
import numpy as np
import abipy.data as data  
import abipy.abilab as abilab


def make_scf_nscf_inputs(structure, paral_kgb=1):
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=2)
    structure = inp.set_structure(structure)

    # Global variables
    global_vars = dict(ecut=6,
                       nband=8,
                       timopt=-1,
                       paral_kgb=0,
                       #nstep=4, # This is not enough to converge. Used to test the automatic restart.
                       nstep=10,
                    )

    inp.set_vars(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    inp[1].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    inp[2].set_kpath(ndivsm=6, kptbounds=kptbounds)
    inp[2].set_vars(tolwfr=1e-12)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()

    return scf_input, nscf_input


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir: 
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 
                                                                                                                         
    # build the structures
    base_structure = abilab.Structure.from_file(data.cif_file("si.cif"))
    modifier = abilab.StructureModifier(base_structure)

    etas = [-0.1, 0, +0.1]
    ph_displ = np.reshape(np.zeros(3*len(base_structure)), (-1,3))
    ph_displ[0,:] = [+1, 0, 0]
    ph_displ[1,:] = [-1, 0, 0]

    displaced_structures = modifier.displace(ph_displ, etas, frac_coords=False)

    flow = abilab.Flow(workdir, manager=options.manager)

    for structure in displaced_structures:
        # Create the work for the band structure calculation.
        scf_input, nscf_input = make_scf_nscf_inputs(structure)
                                                                   
        work = abilab.BandStructureWork(scf_input, nscf_input)
        flow.register_work(work)

    return flow.allocate()


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
