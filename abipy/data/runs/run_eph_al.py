#!/usr/bin/env python
"""
Electron-phonon calculations.
"""
from __future__ import division, print_function, unicode_literals

import os
import sys
import numpy as np
import abipy.data as abidata  
import abipy.abilab as abilab


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Preparatory run for E-PH calculations.
    # The sequence of datasets makes the ground states and
    # all of the independent perturbations of the single Al atom 
    # for the irreducible qpoints in a 4x4x4 grid.
    # Note that the q-point grid must be a sub-grid of the k-point grid (here 8x8x8)
    pseudos = abidata.pseudos("Al.oncvpsp")

    structure = abilab.Structure.from_abivars(
        acell=3*[7.5],
        rprim=[0.0, 0.5, 0.5, 
               0.5, 0.0, 0.5,
               0.5, 0.5, 0.0],
        typat=1,
        xred=[0.0, 0.0, 0.0],
        ntypat=1,
        znucl=13,
    )

    gs_inp = abilab.AbinitInput(structure, pseudos)

    gs_inp.set_vars(
        prtpot=1,
        istwfk="*1",
        ecut=6.0,
        nband=5,
        occopt=7,    # include metallic occupation function with a small smearing
        tsmear=0.04,
        tolvrs=1e7,
    )

    # The kpoint grid is minimalistic to keep the calculation manageable.
    # Use a centered grid for the kpoints
    gs_inp.set_kmesh(
        ngkpt=[8, 8, 8], 
        kptopt=3,
        shiftk=[0.0, 0.0, 0.0],
    )

    # Phonon calculation with 4x4x4
    qpoints = np.reshape([
         0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 
         2.50000000e-01,  0.00000000e+00,  0.00000000e+00,
       #  5.00000000e-01,  0.00000000e+00,  0.00000000e+00,
       #  2.50000000e-01,  2.50000000e-01,  0.00000000e+00,
       #  5.00000000e-01,  2.50000000e-01,  0.00000000e+00,
       # -2.50000000e-01,  2.50000000e-01,  0.00000000e+00,
       #  5.00000000e-01,  5.00000000e-01,  0.00000000e+00,
       # -2.50000000e-01,  5.00000000e-01,  2.50000000e-01,
        ], (-1,3))

    flow = abilab.Flow(workdir, manager=options.manager)
    work0 = flow.register_task(gs_inp, task_class=abilab.ScfTask)

    ph_work = abilab.PhononWork.from_scf_task(work0[0], qpoints)
    flow.register_work(ph_work)

    eph_inp = gs_inp.deepcopy()
    eph_inp.set_vars(
        optdriver=7,
        #getwfk   20      # Read GS wavefunctions from DS20_WFK
        #getddb   20      # Read DDB files from DS20_DDB
        ddb_ngqpt=[4, 4, 4],  # q-mesh used to produce the DDB file (must be consisten with DDB data)
        #eph_intmeth=2,       # Tetra
        #eph_fsewin="0.8 eV",   # Energy window around Ef
        #eph_mustar="0.12",     # mustar parameter
        # q-path for phonons and phonon linewidths.
        #ph_ndivsm 20
        #ph_nqpath 3
        #ph_qpath 
        #  0   0   0
        #  0.5 0   0
        #  0.5 0.5 0
        # phonon DOS obtained via Fourier interpolation
        #ph_intmeth 2          # Tetra for phonon DOS and A2F
        #ph_smear   0.001 eV
        #ph_wstep   0.0001 eV
        #ph_ngqpt   16 16 16   # q-mesh for Fourier interpolatation of IFC and a2F(w)
        #ph_nqshift 1
        #ph_qshift 
        #  0 0 0 
    )

    eph_work = abilab.Work()
    eph_task = eph_work.register_eph_task(eph_inp, deps={work0[0]: "WFK",
        ph_work: ["DDB", "DVDB"]})

    # ph_work: "DVDB"
    flow.register_work(eph_work)
    flow.allocate()

    # EPH does not support autoparal
    eph_task.with_fixed_mpi_omp(1,1)
                                                               
    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    import sys
    sys.exit(main())

"""
    # DATASET 1 : make ground state wavefunctions and density
    #multi[0].set_vars(
    #    tolwfr=1.0e-14
    #    rfphon=0        # for DS1 do _not_ do perturbation
    #    nqpt1=0         # for DS1 do _not_ do perturbation
    #    prtwf1=1        # need GS wavefunctions for further runs
    #    getwfk1=0       # Change default
    #    kptopt1=1       # Use symmetries to produce GS-WFK

    # general data for all phonon calculations:
    multi[1:9].set_vars(
        nqpt=1,            # 1 qpoint at a time
        getwfk=1,          # all other datasets get WFK from DS1
        prtwf=0,
        rfatpol=[1, 1],    # all atoms are perturbed
        rfdir=[1, 1, 1],   # all directions of perturbation
        rfphon=1,
        tolvrs=1e-7,
    )
"""
