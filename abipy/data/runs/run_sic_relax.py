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
    ecut=12,
    istwfk="*1",
    chksymbreak=0,
    #accesswff=3
)

def build_workflows():
    
    all_ecuts = np.arange(20,28,4)
    all_ngkpts = [3 * [nk] for nk in np.arange(2,6,2)]

    structure = abilab.Structure.from_abivars(unit_cell)
    pseudos = data.pseudos("14si.pspnc","6c.pspnc")

    print(structure)

    shiftk = [[0.5,0.5,0.5],[0.5,0.0,0.0],[0.0,0.5,0.0],[0.0,0.0,0.5]]

    works = []    
    for ecut in all_ecuts:
        for ngkpt in all_ngkpts:
    
            gs_inp = abilab.AbiInput(pseudos=pseudos)
            gs_inp.set_structure(structure)
            gs_inp.set_variables(**global_vars)

            gs_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
            gs_inp.tolvrs=1e-16

            
            nscf_inp = abilab.AbiInput(pseudos=pseudos)
            nscf_inp.set_structure(structure)
            nscf_inp.set_variables(**global_vars)

            nscf_inp.set_kpath(ndivsm=20)
            nscf_inp.tolwfr = 1e-22

            #relax_inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
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
            manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1)

            workdir = "relax_SiC/ecut" + str(ecut) + "_ngkpt" + "x".join(str(nk)for nk in ngkpt)
            assert workdir not in [work.workdir for work in works]
            work = abilab.Workflow(workdir, manager)

            gs_link = work.register(gs_inp)  
            nscf_link = work.register(nscf_inp, links=gs_link.produces_exts("_DEN"))
            relax_link = work.register(relax_inp, links=gs_link.produces_exts("_WFK"))

            works.append(work)    

    return works

def main():
    works = build_workflows()
    for work in works:
        work.build_and_pickle_dump()

    return 0


if __name__=="__main__":
    sys.exit(main())

