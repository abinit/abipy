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
             2.50000000E-01,  0.00000000E+00,  0.00000000E+00,
             5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
             2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
             5.00000000E-01,  2.50000000E-01,  0.00000000E+00,
            -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
             5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
            -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
            ]

    qpoints = np.reshape(qpoints, (-1,3))

    # Global variables used both for the GS and the DFPT run.
    global_vars = dict(nband=4,             
                       ecut=3.0,         
                       #ecut=12.0,
                       ngkpt=[4, 4, 4],
                       nshiftk=4,
                       shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
                               0.0, 0.5, 0.0,
                               0.5, 0.0, 0.0,
                               0.5, 0.5, 0.5],
                       #shiftk=[0, 0, 0],
                       paral_kgb=paral_kgb,
                       ixc=1,
                       nstep=25,
                       diemac=9.0,
                    )

    gs_inp = abilab.AbiInput(pseudos=pseudos)
    gs_inp.set_structure(structure)
    gs_inp.set_variables(**global_vars)
    gs_inp.set_variables(tolvrs=1.0e-18)

    # Get the qpoints in the IBZ. Note that here we use a q-mesh with ngkpt=(4,4,4) and shiftk=(0,0,0)
    # i.e. the same parameters used for the k-mesh in gs_inp.
    qpoints = gs_inp.get_ibz(ngkpt=(4,4,4), shiftk=(0,0,0), kptopt=1).qpoints
    #print("get_ibz", qpoints)

    ph_inputs = abilab.AbiInput(pseudos=pseudos, ndtset=len(qpoints))

    for ph_inp, qpt in zip(ph_inputs, qpoints):
        # Response-function calculation for phonons.
        ph_inp.set_structure(structure)
        ph_inp.set_variables(**global_vars)
        ph_inp.set_variables(
            rfphon=1,        # Will consider phonon-type perturbation
            nqpt=1,          # One wavevector is to be considered
            qpt=qpt,         # This wavevector is q=0 (Gamma)
            tolwfr=1.0e-20,
            kptopt=3,
            )

            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            #kptopt   2      # Automatic generation of k points, taking

    # Split input into gs_inp and ph_inputs
    all_inps = [gs_inp] 
    all_inps.extend(ph_inputs.split_datasets())
    return all_inps


def build_flow(options):
    """
    Create a `Flow` for phonon calculations:

        1) One workflow for the GS run.

        2) nqpt works for phonon calculations. Each work contains 
           nirred tasks where nirred is the number of irreducible phonon perturbations
           for that particular q-point.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)

    all_inps = scf_ph_inputs()
    scf_input, ph_inputs = all_inps[0], all_inps[1:]

    return abilab.phonon_flow(workdir, scf_input, ph_inputs, manager=manager)


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
