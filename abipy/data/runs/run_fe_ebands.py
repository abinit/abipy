#!/usr/bin/env python
"""
Calculation of the band structure of Fe with and without magnetization.
See tutorial/Input/tspin_1.in
"""
from __future__ import division, print_function, unicode_literals

import os
import sys
import abipy.data as data  
import abipy.abilab as abilab


def make_scf_nscf_inputs(nsppol):
    """Generate two input files for the GS and the NSCF run for given nsppol"""
    inp = abilab.AbiInput(pseudos=data.pseudos("26fe.pspnc"), ndtset=2)

    # Fe normal bcc structure for test of a ferromagnetic calculation
    structure = data.structure_from_ucell("Fe-fm")
    inp.set_structure(structure)

    # Global variables
    global_vars = dict(nsppol=nsppol,
                       ecut=18,
                       nband=8,
                       occopt=3,
                       tsmear=0.01,
                       paral_kgb=0,
                    )
    if nsppol == 2:
        global_vars.update(spinat=[0.0, 0.0, 4.0])

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[4,4,4], shiftk=[0.5,0.5,0.5])
    inp[1].set_variables(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    inp[2].set_kpath(ndivsm=4)
    inp[2].set_variables(tolwfr=1e-8)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()

    return scf_input, nscf_input


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)

    # Create the Flow.
    flow = abilab.AbinitFlow(workdir=workdir, manager=manager)

    # Create the task defining the calculation and run and register it in the flow
    for nsppol in [1,2]:
        scf_input, nscf_input = make_scf_nscf_inputs(nsppol)
        work = abilab.BandStructureWorkflow(scf_input, nscf_input)
        flow.register_work(work)

    return flow.allocate()


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
