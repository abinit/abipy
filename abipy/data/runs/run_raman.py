#!/usr/bin/env python
from __future__ import division, print_function

import os
import sys 
import numpy as np
import abipy.abilab as abilab
import abipy.data as data  

from abipy.data.runs import Tester, decorate_main

def build_flow():
    pseudos = data.pseudos("14si.pspnc")

    # Get the unperturbed structure.
    base_structure = data.structure_from_ucell("si")

    etas = [-.001, 0, +.001]
    ph_displ = np.reshape(np.zeros(3*len(base_structure)), (-1,3))
    ph_displ[0,:] = [+1, 0, 0]
    ph_displ[1,:] = [-1, 0, 0]

    # Build new structures by displacing atoms according to the phonon displacement
    # ph_displ (in cartesian coordinates). The Displacement is normalized so that 
    # the maximum atomic diplacement is 1 Angstrom and then multiplied by eta.
    modifier = abilab.StructureModifier(base_structure)

    displaced_structures = modifier.displace(ph_displ, etas, frac_coords=False)

    # Initialize flow. Each workflow in the flow defines a complete BSE calculation for given eta.
    workdir = os.path.join(os.path.dirname(__file__), base_structure.formula + "_RAMAN")
    manager = abilab.TaskManager.from_file("taskmanager.yaml")

    flow = abilab.AbinitFlow(workdir, manager)

    # Generate the different shifts to average
    ndiv = 2
    shift1D = np.arange(1,2*ndiv+1,2)/(2*ndiv)
    allshifts = [[x,y,z] for x in shift1D for y in shift1D for z in shift1D]

    for structure, eta in zip(displaced_structures, etas):
        for shift in allshifts:
            flow.register_work(raman_workflow(structure, pseudos, shift))

    return flow.allocate()


def raman_workflow(structure, pseudos, shiftk):
    # Generate 3 different input files for computing optical properties with BSE.

    # Global variables
    global_vars = dict(
        ecut=12,
        istwfk="*1",
        chksymbreak=0,
    )

    # GS run
    scf_inp = abilab.AbiInput(pseudos=pseudos)

    scf_inp.set_structure(structure)
    scf_inp.set_variables(**global_vars)
    scf_inp.set_kmesh(ngkpt=[2,2,2], shiftk=shiftk)
    scf_inp.tolvrs = 1e-6

    # NSCF run
    nscf_inp = abilab.AbiInput(pseudos=pseudos)

    nscf_inp.set_structure(structure)
    nscf_inp.set_variables(**global_vars)
    nscf_inp.set_kmesh(ngkpt=[2,2,2], shiftk=shiftk)

    nscf_inp.set_variables(tolwfr=1e-12,
                           nband=12,
                           nbdbuf=4,
                           iscf=-2,
                           )

    # BSE run with Model dielectric function and Haydock (only resonant + W + v)
    # Note that SCR file is not needed here
    bse_inp = abilab.AbiInput(pseudos=pseudos)

    bse_inp.set_structure(structure)
    bse_inp.set_variables(**global_vars)
    bse_inp.set_kmesh(ngkpt=[2,2,2], shiftk=shiftk)

    bse_inp.set_variables(
        optdriver=99,
        ecutwfn=global_vars["ecut"],
        ecuteps=3,
        inclvkb=2,
        bs_algorithm=2,       # Haydock
        bs_haydock_niter=60,  # No. of iterations for Haydock
        bs_exchange_term=1,
        bs_coulomb_term=21,   # Use model W and full W_GG.
        mdf_epsinf=12.0,
        bs_calctype=1,        # Use KS energies and orbitals to construct L0
        soenergy="0.8 eV",
        bs_coupling=0,
        bs_loband=2,
        nband=8,
        bs_freq_mesh="0 10 0.1 eV",
        bs_hayd_term=0,      # No terminator
    )

    # Build the workflow representing a BSE run with model dielectric function.
    return abilab.BSEMDF_Workflow(scf_inp, nscf_inp, bse_inp)


@decorate_main
def main():
    # Build the list of workflows.
    flow = build_flow()
    flow.build_and_pickle_dump()
    return 0


if __name__ == "__main__":
    sys.exit(main())
