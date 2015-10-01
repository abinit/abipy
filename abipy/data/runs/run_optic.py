#!/usr/bin/env python
"""Optical spectra with Optic."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import abipy.data as data  
import abipy.abilab as abilab


def build_flow(options, paral_kgb=0):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_") 

    multi = abilab.MultiDataset(structure=data.structure_from_ucell("GaAs"),
                                pseudos=data.pseudos("31ga.pspnc", "33as.pspnc"), ndtset=5)

    # Global variables
    kmesh = dict(ngkpt=[4, 4, 4], 
                 nshiftk=4,
                 shiftk=[[0.5, 0.5, 0.5],
                         [0.5, 0.0, 0.0],
                         [0.0, 0.5, 0.0],
                         [0.0, 0.0, 0.5]]
                )

    global_vars = dict(ecut=2, paral_kgb=paral_kgb)
    global_vars.update(kmesh)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_vars(
        tolvrs=1e-6,
        nband=4,
    )

    # NSCF run with large number of bands, and points in the the full BZ
    multi[1].set_vars(
        iscf=-2,
        nband=20, 
        nstep=25,
        kptopt=1,
        tolwfr=1.e-9,
        #kptopt=3,
    )

    # Fourth dataset : ddk response function along axis 1
    # Fifth dataset : ddk response function along axis 2
    # Sixth dataset : ddk response function along axis 3
    for dir in range(3):
        rfdir = 3 * [0]
        rfdir[dir] = 1

        multi[2+dir].set_vars(
            iscf=-3,
            nband=20,  
            nstep=1,
            nline=0,  
            prtwf=3,
            kptopt=3,
            nqpt=1, 
            qpt=[0.0, 0.0, 0.0],
            rfdir=rfdir,
            rfelfd=2,
            tolwfr=1.e-9,
        )

    scf_inp, nscf_inp, ddk1, ddk2, ddk3 = multi.split_datasets()

    # Initialize the flow.
    flow = abilab.Flow(workdir, manager=options.manager)

    bands_work = abilab.BandStructureWork(scf_inp, nscf_inp)
    flow.register_work(bands_work)

    ddk_work = abilab.Work()
    for inp in [ddk1, ddk2, ddk3]:
        ddk_work.register_ddk_task(inp, deps={bands_work.nscf_task: "WFK"})

    flow.register_work(ddk_work)

    # TODO
    # Check is the order of the 1WF files is relevant. Can we use DDK files ordered 
    # in an arbitrary way or do we have to pass (x,y,z)?

    """
    0.002         ! Value of the smearing factor, in Hartree
    0.0003  0.3   ! Difference between frequency values (in Hartree), and maximum frequency ( 1 Ha is about 27.211 eV)
    0.000         ! Scissor shift if needed, in Hartree
    0.002         ! Tolerance on closeness of singularities (in Hartree)
    1             ! Number of components of linear optic tensor to be computed
    11            ! Linear coefficients to be computed (x=1, y=2, z=3)
    2             ! Number of components of nonlinear optic tensor to be computed
    123 222       ! Non-linear coefficients to be computed
    """

    # Optic does not support MPI with ncpus > 1.
    optic_input = abilab.OpticInput(
        broadening=0.002,
        domega=0.0003,
        maxomega=0.3,
        scissor=0.000,
        tolerance=0.002,
        num_lin_comp=1,
        lin_comp=11,
        num_nonlin_comp=2,
        nonlin_comp=(123, 222),
    )

    optic_task = abilab.OpticTask(optic_input, nscf_node=bands_work.nscf_task, ddk_nodes=ddk_work, 
                                  manager=flow.manager.to_shell_manager(mpi_procs=1))
    flow.register_task(optic_task)

    return flow


def optic_flow_from_files():
    # Optic does not support MPI with ncpus > 1.
    manager = abilab.TaskManager.from_user_config()
    manager.set_mpi_procs(1)

    flow = abilab.Flow(workdir="OPTIC_FROM_FILE", manager=manager)
    
    ddk_nodes = [
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_0/outdata/out_1WF",
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_1/outdata/out_1WF",
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_2/outdata/out_1WF",
    ]
    nscf_node = "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_0/task_1/outdata/out_WFK"

    optic_task = abilab.OpticTask(optic_input, nscf_node=nscf_node, ddk_nodes=ddk_nodes)
    flow.register_task(optic_task)
                                                                                                                          
    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    #flow = optic_flow_from_files(options)
    #print("optic manager after allocate", flow[2][0].manager)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
