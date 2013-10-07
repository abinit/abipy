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
    print("pseudos",inp.pseudos)
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
            #symsigma=1,
            gwcalctyp=20,
        )

    kptgw = [
         -2.50000000E-01, -2.50000000E-01,  0.00000000E+00,
         -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
          5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
         -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
          5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
          0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
      ]

    bdgw = [1,8]

    sigma.set_kptgw(kptgw, bdgw)

    return inp.split_datasets()


def scf_ph_inputs():
    # Crystalline AlAs : computation of the second derivative of the total energy
    structure = data.structure_from_ucell("alas")
    pseudos = data.pseudos("13al.981214.fhi", "33as.pspnc")

    # List of q-points for the phonon calculation.
    qpoints = [
             0.00000000E+00,  0.00000000E+00,  0.00000000E+00, 
             2.50000000E-01,  0.00000000E+00,  0.00000000E+00,
             5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
             2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
             5.00000000E-01,  2.50000000E-01,  0.00000000E+00,
            -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
             5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
            -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
            ]
    qpoints = np.reshape(qpoints, (-1,3))

    # Global variables used both for the GS and the DFPT run.
    global_vars = dict(nband=4,             
                       ecut=3.0,         
                       ngkpt=[4, 4, 4],
                       shiftk=[0, 0, 0],
                       tolvrs=1.0e-8,
                    )

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=1+len(qpoints))

    inp.set_structure(structure)
    inp.set_variables(**global_vars)

    for i, qpt in enumerate(qpoints):
        # Response-function calculation for phonons.
        inp[i+2].set_variables(
            rfphon=1,        # Will consider phonon-type perturbation
            nqpt=1,          # One wavevector is to be considered
            qpt=qpt,         # This wavevector is q=0 (Gamma)
            )
            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            #kptopt   2     # Automatic generation of k points, taking

    # return gs_inp, ph_inputs
    return inp.split_datasets()


def ph_flow():
    workdir = "PHONONS"

    all_inps  = scf_ph_inputs()

    scf_input, ph_inputs = all_inps[0], all_inps[1:]
                                                                        
    #manager = abilab.TaskManager.from_file("taskmanager.yaml")
    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=dict(autoparal=0, max_ncpus=1))
    #manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=dict(autoparal=1, max_ncpus=2))

    flow = phonon_flow(workdir, manager, scf_input, ph_inputs)

    flow.build_and_pickle_dump()
    return 0


def gw_flow():
    workdir = "WORKS"

    gs, nscf, scr_input, sigma_input = all_inputs()
                                                                        
    #manager = abilab.TaskManager.from_file("taskmanager.yaml")
    #manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=dict(autoparal=0, max_ncpus=1))
    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=dict(autoparal=1, max_ncpus=2))

    flow = g0w0_flow_with_qptdm(workdir, manager, gs, nscf, scr_input, sigma_input)

    from pymatgen.io.abinitio.tasks import GW0_Task
    flow[2][0].__class__ = GW0_Task
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

    #wfk_file = os.path.join(os.getcwd(), "out_WFK")
    #qptdm_work = qptdm_workflow(wfk_file, scr_input, workdir, manager)

    #qptdm_work.build_and_pickle_dump()
    #return 0

    flow = g0w0_flow_with_qptdm(workdir, manager, gs, nscf, scr_input, sigma_input)
    flow.build_and_pickle_dump()

    return 0

@decorate_main
def main():
    # QPTDM
    qptdm_work()

    # GW Works
    #gw_flow()

    # Phonon Works
    #ph_flow()


if __name__ == "__main__":
    import sys
    sys.exit(main())
