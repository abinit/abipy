#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.abilab as abilab
import abipy.data as data  

from abipy.data.runs import decorate_main

from abipy.data.runs.qptdm_workflow import *

def make_base_input():
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))
    inp.set_structure(data.structure_from_ucell("si"))

    # Global variables
    global_vars = dict(ecut=12,
                       nsppol=1,
                       nband=20,
                       timopt=-1,
                       paral_kgb=0,
                       chksymbreak=0,
                       prtwf=0
                       #accesswff=3,
                    )

    inp.set_variables(**global_vars)

    # Simple GS run.
    inp.set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0])
    #inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    #inp.set_kmesh(ngkpt=[12,12,12], shiftk=[0.1,0.2,0.3])
    inp.tolvrs = 1e-8

    return inp


def all_inputs():
    structure = abilab.Structure.from_file(data.cif_file("si.cif"))
    pseudos = data.pseudos("14si.pspnc")

    ecut = 6

    global_vars = dict(
        ecut=ecut,
        timopt=-1,
        istwfk = "*1",
    )

    gs = abilab.AbiInput(pseudos=pseudos)
    gs.set_structure(structure)

    # This grid is the most economical, but does not contain the Gamma point.
    gs_kmesh = dict(
        ngkpt=[2,2,2],
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap. 
    gw_kmesh = dict(
        ngkpt=[2,2,2],
        shiftk=[0.0, 0.0, 0.0,  
                0.0, 0.5, 0.5,  
                0.5, 0.0, 0.5,  
                0.5, 0.5, 0.0]
    )

    # Dataset 1 (GS run)
    gs.set_variables(**global_vars)
    gs.set_kmesh(**gs_kmesh)
    gs.set_variables(tolvrs=1e-6,
                     nband=4,
                    )

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    nscf = abilab.AbiInput(pseudos=pseudos)
    nscf.set_structure(structure)
    nscf.set_variables(**global_vars)
    nscf.set_kmesh(**gw_kmesh)

    nscf.set_variables(iscf=-2,
                       getden=1,
                       tolwfr=1e-12,
                       nband=35,
                       nbdbuf=5,
                       )

    # Dataset3: Calculation of the screening.
    scr = abilab.AbiInput(pseudos=pseudos)
    scr.set_structure(structure)
    scr.set_kmesh(**gw_kmesh)
    scr.set_variables(**global_vars)

    scr.set_variables(
        optdriver=3,   
        getkss=2,      
        nband=25,    
        ecutwfn=ecut,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
        ppmfrq="16.7 eV",
    )

    return gs, nscf, scr

@decorate_main
def main():
    gs, nscf, scr_input = all_inputs()

    policy = dict(autoparal=0, max_ncpus=2)
    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=policy)

    workdir = "QPTDM"

    # This is to produce the out_WFK file
    #wfk_work = abilab.Workflow(workdir, manager)
    #gs_link = wfk_work.register(gs)
    #nscf_link = wfk_work.register(nscf, links=gs_link.produces_exts("DEN"))
    #wfk_work.start()
    #return 

    wfk_file = os.path.join(os.getcwd(), "out_WFK")
    qptdm_work = build_qptdm_workflow(workdir, manager, scr_input, wfk_file)
    #qptdm_work.connect(nscf_link.produces_exts("WFK"))
    #qptdm_work.connect(nscf_link, exts="WFK")

    qptdm_work.build_and_pickle_dump()

    return 

    works = AbinitWorks(workdir, manager)

    # One can register a workflow object.
    wfk_link = works.register(wfk_work)

    # Register a function that will be executed to build another workflow
    #works.register(qptdm_work, links=wfk_link.produces_exts("WFK"))
    #works.build_and_pickle_dump()

if __name__ == "__main__":
    import sys
    sys.exit(main())
