#!/usr/bin/env python
"""
This script shows how to perform a structural relaxation in two steps:

    1) Relaxation of atomic positions with unit cell parameters fixed.

    2) Full relaxation (atoms + cell) with the initial configuration read from step 1)
"""
from __future__ import division, print_function

import sys
import os
import abipy.data as abidata  
import abipy.abilab as abilab


def make_ion_ioncell_inputs():
    cif_file = abidata.cif_file("si.cif")
    structure = abilab.Structure.from_file(cif_file)

    # Perturb the structure (random perturbation of 0.1 Angstrom)
    structure.perturb(distance=0.1)

    pseudos = abidata.pseudos("14si.pspnc")

    global_vars = dict(
        ecut=4,  
        ngkpt=[4,4,4], 
        shiftk=[0,0,0],
        nshiftk=1,
        chksymbreak=0,
        paral_kgb=0,
    )

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    inp.set_structure(structure)

    # Global variables
    inp.set_variables(**global_vars)

    # Dataset 1 (Atom Relaxation)
    inp[1].set_variables(
        optcell=0,
        ionmov=2,
        tolrff=0.02,
        tolmxf=5.0e-5,
        ntime=50,
        #ntime=5, To test the restart
        dilatmx=1.1, # FIXME: abinit crashes if I don't use this
    )

    # Dataset 2 (Atom + Cell Relaxation)
    inp[2].set_variables(
        optcell=1,
        ionmov=2,
        ecutsm=0.5,
        dilatmx=1.1,
        tolrff=0.02,
        tolmxf=5.0e-5,
        strfact=100,
        ntime=50,
        #ntime=5, To test the restart
        )

    ion_inp, ioncell_inp = inp.split_datasets()
    return ion_inp, ioncell_inp


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)

    # Create the flow
    flow = abilab.AbinitFlow(workdir, manager)

    # Create a relaxation workflow and add it to the flow.
    ion_inp, ioncell_inp = make_ion_ioncell_inputs()

    work = abilab.RelaxWorkflow(ion_inp, ioncell_inp)
    flow.register_work(work)

    #bands_work = abilab.BandStructureWorkflow(scf_input, nscf_input)
    return flow.allocate()


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())

