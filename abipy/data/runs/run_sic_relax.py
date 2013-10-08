#!/usr/bin/env python
from __future__ import division, print_function

import sys
import numpy as np
import os

import abipy.abilab as abilab
import abipy.data as data

from abipy.data.runs import Tester, decorate_main

def build_flow():
    global_vars = dict(
        istwfk="*1",
        chksymbreak=0,
        #accesswff=3
    )
    all_ecuts = np.arange(20,28,4)
    all_ngkpts = [3 * [nk] for nk in np.arange(2,6,2)]
    shiftk = [[0.5,0.5,0.5],[0.5,0.0,0.0],[0.0,0.5,0.0],[0.0,0.0,0.5]]

    pseudos = data.pseudos("14si.pspnc","6c.pspnc")
    structure = data.structure_from_ucell("sic")
    print(structure)

    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1)
    workdir = "relax_SiC"

    flow = abilab.AbinitFlow(workdir, manager)

    for ecut in all_ecuts:
        for ngkpt in all_ngkpts:
    
            gs_inp = abilab.AbiInput(pseudos=pseudos)
            gs_inp.set_structure(structure)
            gs_inp.set_variables(**global_vars)

            gs_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
            gs_inp.tolvrs=1e-16
            gs_inp.ecut = ecut
            
            nscf_inp = abilab.AbiInput(pseudos=pseudos)
            nscf_inp.set_structure(structure)
            nscf_inp.set_variables(**global_vars)

            nscf_inp.set_kpath(ndivsm=20)
            nscf_inp.tolwfr = 1e-22
            nscf_inp.ecut = ecut

            relax_inp = abilab.AbiInput(pseudos=pseudos)
            relax_inp.set_structure(structure)
            relax_inp.set_variables(**global_vars)

            relax_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
            relax_inp.toldff = 1e-6
            relax_inp.tolmxf = 1e-5
            relax_inp.strfact = 100
            relax_inp.ecutsm = 0.5
            relax_inp.dilatmx = 1.15
            relax_inp.ntime = 100
            relax_inp.ecut = ecut

            #relax_inp[1].set_variables(
            #    ionmov = 2,
            #    optcell = 0,
            #)

            #relax_inp[2].set_variables(
            #    ionmov = 2,
            #    optcell = 1,
            #)
  
            relax_inp.set_variables(
                ionmov = 2,
                optcell = 1,
            )

            # Initialize the workflow.
            work = abilab.Workflow()
            gs_task = work.register(gs_inp)  

            nscf_task = work.register(nscf_inp, deps={gs_task: "DEN"})

            work.register(relax_inp, deps={gs_task: "WFK"})

            flow.register_work(work)  

    return flow.allocate()

@decorate_main
def main():
    flow = build_workflow()
    flow.build_and_pickle_dump()
    return 0


if __name__=="__main__":
    sys.exit(main())

