#!/usr/bin/env python
"""
This script shows how to compute the band structure of Fe with and without magnetization.
See tutorial/Input/tspin_1.in
"""
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import Tester, enable_logging

def make_scf_nscf_inputs(nsppol):
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
                    )
    if nsppol == 2:
        global_vars.update(spinat=[0.0, 0.0, 4.0])

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[4,4,4], shiftk=[0.5,0.5,0.5])
    inp[1].set_variables(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    inp[2].set_kpath(ndivsm=6)
    inp[2].set_variables(tolwfr=1e-10)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()

    return scf_input, nscf_input

def build_bands_flow(workdir):
    manager = abilab.TaskManager.from_user_config()

    flow = abilab.AbinitFlow(workdir, manager)

    # Create the task defining the calculation and run.
    for nsppol in [1,2]:
        scf_input, nscf_input = make_scf_nscf_inputs(nsppol)
        work = abilab.BandStructureWorkflow(scf_input, nscf_input)
        flow.register_work(work)

    return flow.allocate()


@enable_logging
def main():
    tester = Tester()
    flow = build_bands_flow(tester.workdir)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    import sys
    sys.exit(main())
