#!/usr/bin/env python
"""LDA+U band structure of NiO for several values of U-J."""
from __future__ import division, print_function

import sys
import os
import numpy as np
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import AbipyTest, MixinTest


class LdausFlowTest(AbipyTest, MixinTest):
    """
    Unit test for the flow defined in this module.  
    Users who just want to learn how to use this flow can ignore this section.
    """
    def setUp(self):
        super(LdausFlowTest, self).setUp()
        self.init_dirs()
        self.flow = build_flow()


def make_scf_nscf_dos_inputs(structure, pseudos, luj_params):
    # Input file taken from tldau_2.in
    inp = abilab.AbiInput(pseudos=pseudos, ndtset=3)
    inp.set_structure(structure)

    # Global variables
    global_vars = dict(
        # 
        ecut=15,
        pawecutdg=30,
        nband=40,
        occopt=7,
        tsmear=0.015,
        nstep=50,
        paral_kgb=0,
        #
        # Spin
        nsppol=1,
        nspden=2,
        nspinor=1,
        spinat=[0,  0,  1,
                0,  0, -1,
                0,  0,  0,
                0,  0,  0],
        # Kpoint Grid
        # The k point grid is not symmetric, but the calculations being 
        # for the ground-state, this is not a problem.
    )

    inp.set_variables(**global_vars)
    inp.set_variables(**luj_params.to_abivars())

    # GS run.
    inp[1].set_variables(
        iscf=17,
        toldfe=1.0e-8,
        ngkpt=[2, 2, 2],
        chksymbreak=0, 
    )

    # Band structure run.
    inp[2].set_kpath(ndivsm=6)
    inp[2].set_variables(tolwfr=1e-10)

    # Dos calculation.
    inp[3].set_variables(
        iscf=-3,   # NSCF calculation
        ngkpt=structure.calc_ngkpt(nksmall=8),      
        shiftk=[0.0, 0.0, 0.0],
        nshiftk=1,
        towlfr=1.e-8,
        #pawprtdos=1,
    )

    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input, dos_input = inp.split_datasets()

    return scf_input, nscf_input, dos_input


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed)
    workdir = os.path.basename(os.path.abspath(__file__).replace(".py", "")) if not options.workdir else options.workdir

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else options.manager

    flow = abilab.AbinitFlow(workdir, manager)

    # Create the workflow for the band structure calculation.
    structure = data.structure_from_ucell("NiO")
    pseudos = data.pseudos("28ni.paw", "8o.2.paw")

    # The code below set up the parameters for the LDA+U calculation in NiO.
    #usepawu   1
    #lpawu   2 -1
    #upawu  8.0 0.0 eV
    #jpawu  0.8 0.0 eV
    usepawu = 1
    u_values = [5.0, 8.0]

    for u in u_values:
        # Apply U-J on Ni only.
        luj_params = abilab.LdauParams(usepawu, structure)
        luj_params.luj_for_symbol("Ni", l=2, u=u, j=0.1*u, unit="eV")

        scf_input, nscf_input, dos_input = make_scf_nscf_dos_inputs(structure, pseudos, luj_params)
                                                                       
        work = abilab.BandStructureWorkflow(scf_input, nscf_input, dos_inputs=dos_input)
        flow.register_work(work)

    return flow.allocate()


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
