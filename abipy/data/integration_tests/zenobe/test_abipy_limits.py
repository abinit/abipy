#!/usr/bin/env python
"""
This script shows how to perform a RAMAN calculation with excitonic effects 
included with the BSE formalism.
"""
from __future__ import division, print_function

import sys 
import os
import numpy as np
from abipy import abilab
import abipy.data as data  
from pymatgen.io.abinit.tasks import TaskPolicy
from pymatgen.io.abinit.abiobjects import KSampling
from abipy.dfpt import PhononBands, PhdosReader, PhdosFile

unit_cell = dict(
            acell=3*[10.20],
            natom=2,
           ntypat=1,
            typat=[1,1],
            znucl=[14],
             xred=[[0,0,0],
                   [1/4,1/4,1/4]],
	    rprim=[[0,0.5,0.5],
                   [0.5,0,0.5],
                   [0.5,0.5,0]],
)

global_vars = dict(
    paral_kgb=0,
    ecut=8,
    nstep=1,
    diemac=12,
    ixc=7,
    chksymbreak=0,
    chkprim=0,
#    nsym=1,
)

def raman_flow():

    # Get the unperturbed structure.
    base_structure = abilab.Structure.from_abivars(unit_cell)

    pseudos=["14si.pspnc"]

    workdir = os.path.join(os.path.dirname(__file__), "test_abipy_new")

    manager = abilab.TaskManager.from_user_config()
    #manager = abilab.TaskManager.from_file("bfo_manager.yml")

    policy = TaskPolicy(autoparal=0)
    gs_manager = manager.deepcopy()

    # Initialize flow. Each workflow in the flow defines a complete BSE calculation for given eta.
    flow = abilab.Flow(workdir, manager)

    # There will be kppa/natom kpoints in the unit cell !
    kppa = 3456  # ngkpt = [12,12,12] for the primitive cell
    kppa_gs = 1728 # ngkpt = [8,8,8] for the primitive cell

    etas = [-1,0,1] # Anyway, it rescales everything at the end :-)
    eta=0.01

    scale_matrix = [[-1,0,1],[-1,1,0],[-1,-1,0]]

    ph_tot = np.array([[0.01,0.011,0.021],[-0.01,0.041,-0.02]])
    modifier = abilab.StructureModifier(base_structure)
           
    displaced_structure = modifier.frozen_phonon([0.5,0.5,0.5],ph_tot,do_real=True,frac_coords=False,scale_matrix=scale_matrix)

    structure = displaced_structure

    ksampgs = KSampling.automatic_density(structure,kppa_gs,chksymbreak=0,shifts=[0,0,0])

    gs_inp = gs_input(structure,pseudos, ksampgs)
    wflow = abilab.Work()
    gs_t = wflow.register_scf_task(gs_inp)
    gs_t.set_manager(gs_manager)
    flow.register_work(wflow,workdir="gs_task")

    ksamp = KSampling.automatic_density(structure,kppa,chksymbreak=0,shifts=[1/4,1/4,1/4])
    flow.register_work(raman_workflow(structure, pseudos, gs_t, ksamp),workdir="bse_task")

    return flow.allocate()

def gs_input(structure, pseudos, ksamp):

    inp = abilab.AbinitInput(structure,pseudos=pseudos)

    inp.set_vars(**global_vars)

    inp.set_vars(
	tolvrs=1e-16,
        nband=2*len(structure), # 2 bands / atoms
        nstep=500,
        kptopt=1,
        iscf=5,
    )

    vars_ksamp = ksamp.to_abivars()
    vars_ksamp.pop("#comment",None)

    inp.set_vars(**vars_ksamp)

    return inp        

def raman_workflow(structure, pseudos, scf_t, ksamp):
    # Generate 3 different input files for computing optical properties with BSE.

    inp = abilab.MultiDataset(structure,pseudos=pseudos, ndtset=2)

    inp.set_vars(**global_vars)

    vars_ksamp = ksamp.to_abivars()
    vars_ksamp.pop("#comment",None)

    inp.set_vars(**vars_ksamp)

    # NSCF run
    inp[0].set_vars(
        iscf=-2,
       nband=10*len(structure),
       nbdbuf=2*len(structure),
       nstep=500,
       tolwfr=1.e-22,
    )

    inp[1].set_vars(
      gwmem=10,
      gwpara=2,
      optdriver=99,
      nband=5*len(structure), # 10 bands for 2 atoms
      bs_loband=1*len(structure), # start at 2 for 2 atoms
      bs_algorithm=2, # Haydock
      bs_calctype=1, # KS wavefunctions
      bs_coulomb_term=21, # Use model dielectric function
      mdf_epsinf=12.0,
      bs_exchange_term=1, # Exchange term included
      bs_freq_mesh="0 10 0.01 eV",
      bs_hayd_term=0,
      bs_coupling=0,
      bs_haydock_tol="-0.01 0",
      bs_haydock_niter=1000,
      inclvkb=2,
      ecuteps=4,
      soenergy="0.8 eV",
      ecutwfn=global_vars["ecut"],
    )
    
    nscf_inp, bse_inp = inp.split_datasets()

    workflow = abilab.Work()
    nscf_t = workflow.register(nscf_inp, deps={scf_t: "DEN"}, task_class=abilab.NscfTask)
    bse_t = workflow.register_bse_task(bse_inp, deps={nscf_t: "WFK"})

    return workflow

def main():
    # Define the flow, build files and dirs 
    # and save the object in cpickle format.
    flow = raman_flow()
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
