#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.abilab as abilab
import abipy.data as data  

from abipy.data.runs import decorate_main

from abipy.data.runs.qptdm_workflow import *

def all_inputs():
    structure = abilab.Structure.from_file(data.cif_file("si.cif"))
    pseudos = data.pseudos("14si.pspnc")

    ecut = ecutwfn = 6

    global_vars = dict(
        ecut=ecut,
        timopt=-1,
        istwfk = "*1",
    )

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=4)
    inp.set_structure(structure)
    inp.set_variables(**global_vars)

    gs, nscf, scr, sigma = inp[1:]

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
    gs.set_kmesh(**gs_kmesh)
    gs.set_variables(tolvrs=1e-6,
                     nband=4,
                    )

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    nscf.set_kmesh(**gw_kmesh)

    nscf.set_variables(iscf=-2,
                       tolwfr=1e-12,
                       nband=35,
                       nbdbuf=5,
                       )

    # Dataset3: Calculation of the screening.
    scr.set_kmesh(**gw_kmesh)

    scr.set_variables(
        optdriver=3,   
        nband=25,    
        ecutwfn=ecutwfn,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
        ppmfrq="16.7 eV",
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    sigma.set_kmesh(**gw_kmesh)

    sigma.set_variables(
            optdriver=4,
            nband=25,      
            ecutwfn=ecutwfn,
            ecuteps=4.0,
            ecutsigx=6.0,
            symsigma=1,
        )

    kptgw = [
         -2.50000000E-01, -2.50000000E-01,  0.00000000E+00,
         #-2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
         # 5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
         #-2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
         # 5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
         # 0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
      ]

    bdgw = [1,3]

    sigma.set_kptgw(kptgw, bdgw)

    #return gs, nscf, scr, sigma
    return inp.split_datasets()


def gw_works():
    workdir = "WORKS"

    gs, nscf, scr_input, sigma_input = all_inputs()
                                                                        
    policy = dict(autoparal=0, max_ncpus=2)
    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=policy)

    flow = g0w0_flow_with_qptdm(workdir, manager, gs, nscf, scr_input, sigma_input)
    flow.build_and_pickle_dump()

    return 0

    
def qptdm_work():
    workdir = "QPTDM"

    gs, nscf, scr_input, sigma_input = all_inputs()
                                                                        
    policy = dict(autoparal=0, max_ncpus=2)
    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=policy)

    # This is to produce the out_WFK file
    #wfk_work = Workflow(workdir, manager)
    #gs_link = wfk_work.register(gs)
    #nscf_link = wfk_work.register(nscf, deps=gs_link.produces_exts("DEN"))
    #wfk_work.start()
    #return 

    wfk_file = os.path.join(os.getcwd(), "out_WFK")
    qptdm_work = qptdm_workflow(wfk_file, scr_input, workdir, manager)

    qptdm_work.build_and_pickle_dump()
    return 0

@decorate_main
def main():
    # QPTDM
    #qptdm_work()

    # GW Works
    gw_works()


if __name__ == "__main__":
    import sys
    sys.exit(main())
