#!/usr/bin/env python
"""Phonon band structure of AlAs."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import numpy as np
import abipy.abilab as abilab
import abipy.data as abidata  


def scf_ph_inputs(paral_kgb=0):
    """
    This function constructs the input files for the phonon calculation: 
    GS input + the input files for the phonon calculation.
    """
    # Crystalline AlAs: computation of the second derivative of the total energy
    structure = abidata.structure_from_ucell("AlAs")
    pseudos = abidata.pseudos("13al.981214.fhi", "33as.pspnc")

    # List of q-points for the phonon calculation.
    qpoints = [
             0.00000000E+00,  0.00000000E+00,  0.00000000E+00, 
    ]

    qpoints = np.reshape(qpoints, (-1,3))

    # Global variables used both for the GS and the DFPT run.
    global_vars = dict(nband=4,             
                       ecut=3.0,         
                       ngkpt=[4, 4, 4],
                       shiftk=[0, 0, 0],
                       tolvrs=1.0e-8,
                       paral_kgb=paral_kgb,
                    )

    multi = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=1+len(qpoints))

    multi.set_structure(structure)
    multi.set_vars(global_vars)

    for i, qpt in enumerate(qpoints):
        # Response-function calculation for phonons.
        multi[i+1].set_vars(
            nstep=20,
            rfphon=1,                     # Will consider phonon-type perturbation
            nqpt=1,                       # One one wavevector is to be considered
            qpt=qpt,                      # The q-point
            rfatpol=[1, len(structure)],  # Only the first atom is displaced
            rfdir=[1, 1, 1],              # Along the first reduced coordinate axis
            #kptopt   2                   # Automatic generation of k points, taking
            )


    # Split input into gs_inp and ph_inputs
    return multi.split_datasets()


def build_flow(options):
    """
    Create an `AbinitFlow` for phonon calculations:

        1) One workflow for the GS run.

        2) nqpt workflows for phonon calculations. Each workflow contains 
           nirred tasks where nirred is the number of irreducible phonon perturbations
           for that particular q-point.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    all_inps = scf_ph_inputs()
    scf_input, ph_inputs = all_inps[0], all_inps[1:]

    flow = abilab.Flow(workdir, manager=options.manager)
    from pymatgen.io.abinit.works import build_oneshot_phononwork
    work = build_oneshot_phononwork(scf_input, ph_inputs)
    flow.register_work(work)

    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
