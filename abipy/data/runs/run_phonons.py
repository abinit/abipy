#!/usr/bin/env python
"""This script shows how to compute the phonon band structure of AlAs."""
from __future__ import division, print_function

import sys
import os
import abipy.abilab as abilab
import abipy.data as data  

from abipy.data.runs import enable_logging
from abipy.data.runs.qptdm_workflow import *


def scf_ph_inputs():
    """
    This function constructs the input files for the phonon calculation: 
    GS input + the input files for the phonon calculation.
    """
    # Crystalline AlAs: computation of the second derivative of the total energy
    structure = data.structure_from_ucell("AlAs")
    pseudos = data.pseudos("13al.981214.fhi", "33as.pspnc")

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
                       ngkpt=[4, 4, 4],
                       shiftk=[0, 0, 0],
                       tolvrs=1.0e-8,
                    )

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=1+len(qpoints))

    inp.set_structure(structure)
    inp.set_variables(**global_vars)

    for i, qpt in enumerate(qpoints):
        # Response-function calculation for phonons.
        inp[i+2].set_variables(
            #nstep=2,
            nstep=20,
            rfphon=1,        # Will consider phonon-type perturbation
            nqpt=1,          # One wavevector is to be considered
            qpt=qpt,         # This wavevector is q=0 (Gamma)
            )
            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            #kptopt   2      # Automatic generation of k points, taking

    # Split input into gs_inp and ph_inputs
    return inp.split_datasets()


def ph_flow():
    """
    Create an `AbinitFlow` for phonon calculations:

        1) One workflow for the GS run.

        2) nqpt workflows for phonon calculations. Each workflow contains 
           nirred tasks where nirred is the number of irreducible phonon perturbations
           for that particular q-point.
    """
    workdir = "PHONONS"
    manager = abilab.TaskManager.from_user_config()

    all_inps = scf_ph_inputs()
    scf_input, ph_inputs = all_inps[0], all_inps[1:]

    flow = abilab.phonon_flow(workdir, manager, scf_input, ph_inputs)
    return flow


@enable_logging
def main():
    """Build the flow for Phonon calculations and save the object in cpickle format."""
    flow = ph_flow()
    for task in flow.iflat_tasks():
        print(task, task.manager)

    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
