#!/usr/bin/env python
from __future__ import division, print_function

import sys
import numpy as np
import os

import abipy.abilab as abilab
import abipy.data as data

# One can create a dictionary mapping keywords to values 
unit_cell = dict(
    acell=3*[8.19],
    rprim=[[.0, .5, .5],
               [.5, .0, .5],
               [.5, .5, .0]],
    ntypat=2,
    znucl=[6,14],
    natom=2,
    typat=[1, 2],
    xred=[ [.0, .0, .0],
           [.25,.25,.25] ]
)

global_vars = dict(
    istwfk="*1",
    paral_kgb=0,
    irdwfk=1,
    gwpara=2
    #accesswff=3
)

def build_workflows():
    
    all_ecuteps = np.concatenate([np.arange(20,32,2)])

    ecut = 30
    ngkpt = [4,4,4]

    structure = abilab.Structure.from_abivars(unit_cell)
    pseudos = data.pseudos("14si.pspnc","6c.pspnc")

    print(structure)

    shiftk = [[0.5,0.5,0.5],[0.5,0.0,0.0],[0.0,0.5,0.0],[0.0,0.0,0.5]]

    # Initialize the workflow.
    manager = abilab.TaskManager.from_user_config()
    workdir="gw_SiC_convecuteps_2"

    flow = abilab.AbinitFlow(manager=manager,workdir=workdir)

    wfk_file = "/home/naps/ygillet/SiC/gw_SiC_wfk/work_0/task_0/outdata/out_WFK"

    for ecuteps in all_ecuteps:
    
        scr_inp = abilab.AbiInput(pseudos=pseudos)
        scr_inp.set_structure(structure)
        scr_inp.set_variables(**global_vars)
        scr_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
        scr_inp.set_variables(
        	optdriver=3,
        	nband=500,
        	ecutwfn=ecut,
		ecut=ecut,
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
        	nband=500,
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
    works = build_workflows()
    works.build_and_pickle_dump()

    return 0


if __name__=="__main__":
    sys.exit(main())

