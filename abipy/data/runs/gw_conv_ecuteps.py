#!/usr/bin/env python
from __future__ import division, print_function

import sys
import numpy as np
import os

import abipy.abilab as abilab
import abipy.data as abidata

global_vars = dict(
    istwfk="*1",
    paral_kgb=0,
    irdwfk=1,
    gwpara=2
    #accesswff=3
)

ecut = 5
ngkpt = [4,4,4]
shiftk = [[0.5,0.5,0.5],[0.5,0.0,0.0],[0.0,0.5,0.0],[0.0,0.0,0.5]]

structure = abidata.structure_from_ucell("SiC")
print(structure)
pseudos = abidata.pseudos("14si.pspnc","6c.pspnc")

manager = abilab.TaskManager.from_user_config()


def build_bands_flow():
    # Get the WFK file.
    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)

    inp.set_structure(structure)
    inp.set_variables(**global_vars)
    inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)

    inp[1].set_variables(
		ecut=ecut,
        nband=10,
        tolvrs=1.e-8,
    )

    inp[2].set_variables(
		ecut=ecut,
        nband=100,
        tolwfr=1.e-8,
        iscf=-2
    )

    # Get the SCF and the NSCF input.
    scf_input, nscf_input = inp.split_datasets()

    workdir = "GS"

    #flow = abilab.AbinitFlow(manager=manager, workdir="GS")

    # Instantiate the TaskManager from `taskmanager.yml`.
    #manager = abilab.TaskManager.from_user_config()
    # Build the flow.
    flow = abilab.bandstructure_flow(workdir, manager, scf_input, nscf_input)
    #flow.register_task(inp)
    return flow.allocate()


def build_workflows():
    
    all_ecuteps = np.concatenate([np.arange(2,8,2)])
    print(all_ecuteps)

    # Initialize the workflow.

    workdir="gw_SiC_convecuteps_2"

    flow = abilab.AbinitFlow(manager=manager,workdir=workdir)

    #wfk_file = "/home/naps/ygillet/SiC/gw_SiC_wfk/work_0/task_0/outdata/out_WFK"
    wfk_file = os.path.abspath("GS/work_0/task_1/outdata/out_WFK")

    for ecuteps in all_ecuteps:
    
        scr_inp = abilab.AbiInput(pseudos=pseudos)
        scr_inp.set_structure(structure)
        scr_inp.set_variables(**global_vars)
        scr_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
        scr_inp.set_variables(
        	optdriver=3,
		    ecut=ecut,
        	ecutwfn=ecut,
        	nband=25,
        	symchi=1,
        	inclvkb=0,
        	ecuteps=ecuteps,
        	ppmfrq="16.7 eV",
        )
        
        sig_inp = abilab.AbiInput(pseudos=pseudos)
        sig_inp.set_structure(structure)
        sig_inp.set_variables(**global_vars)
        sig_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
        sig_inp.set_variables(
        	optdriver=4,
        	nband=25,
        	ecutwfn=ecut,
            ecutsigx=(4*ecut),
		    ecut=ecut,
        	symchi=1,
        	inclvkb=0,
        	ecuteps=ecuteps,
        	ppmfrq="16.7 eV",
        )

        work = abilab.Workflow()

        scr_task = work.register(scr_inp)
        sig_task = work.register(sig_inp, deps={scr_task:"SCR"})

        flow.register_work(work)

        flow.allocate()
        work.build()
        scr_task.inlink_file(wfk_file)
        sig_task.inlink_file(wfk_file)
           
    return flow.allocate()

def main():
    #flow = build_bands_flow()
    flow = build_workflows()
    return flow.build_and_pickle_dump()


if __name__=="__main__":
    sys.exit(main())

