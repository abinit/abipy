#!/usr/bin/env python
"""Optical spectra with Optic."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import abipy.data as data  
import abipy.abilab as abilab


optic_input = """\
0.002         ! Value of the smearing factor, in Hartree
0.0003  0.3   ! Difference between frequency values (in Hartree), and maximum frequency ( 1 Ha is about 27.211 eV)
0.000         ! Scissor shift if needed, in Hartree
0.002         ! Tolerance on closeness of singularities (in Hartree)
1             ! Number of components of linear optic tensor to be computed
11            ! Linear coefficients to be computed (x=1, y=2, z=3)
2             ! Number of components of nonlinear optic tensor to be computed
123 222       ! Non-linear coefficients to be computed
"""

def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)

    structure = data.structure_from_ucell("GaAs")

    inp = abilab.AbiInput(pseudos=data.pseudos("31ga.pspnc", "33as.pspnc"), ndtset=5)

    inp.set_structure(structure)

    # Global variables
    kmesh = dict(ngkpt=[4, 4, 4], 
                 nshiftk=4,
                 shiftk=[[0.5, 0.5, 0.5],
                         [0.5, 0.0, 0.0],
                         [0.0, 0.5, 0.0],
                         [0.0, 0.0, 0.5]]
                )

    global_vars = dict(ecut=2,
                       paral_kgb=0,
                      )

    global_vars.update(kmesh)

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_variables(
        tolvrs=1e-6,
        nband=4,
    )

    # NSCF run with large number of bands, and points in the the full BZ
    inp[2].set_variables(
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

        inp[3+dir].set_variables(
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

    scf_inp, nscf_inp, ddk1, ddk2, ddk3 = inp.split_datasets()

    # Initialize the flow.
    flow = abilab.AbinitFlow(workdir, manager)

    bands_work = abilab.BandStructureWorkflow(scf_inp, nscf_inp)
    flow.register_work(bands_work)

    ddk_work = abilab.Workflow()
    for inp in [ddk1, ddk2, ddk3]:
        ddk_work.register_ddk_task(inp, deps={bands_work.nscf_task: "WFK"})

    flow.register_work(ddk_work)

    # TODO
    # Check is the order of the 1WF files is relevant. Can we use DDK files ordered 
    # in an arbitrary way or do we have to pass (x,y,z)?

    # Optic does not support MPI with ncpus > 1.
    shell_manager = manager.to_shell_manager(mpi_procs=1)

    optic_task = abilab.OpticTask(optic_input, nscf_node=bands_work.nscf_task, ddk_nodes=ddk_work, 
                                  manager=shell_manager)
    flow.register_task(optic_task)

    return flow.allocate()


def optic_flow_from_files():
    # Optic does not support MPI with ncpus > 1.
    manager = abilab.TaskManager.from_user_config()
    manager.set_mpi_procs(1)

    flow = abilab.AbinitFlow(workdir="OPTIC_FROM_FILE", manager=manager)
    
    ddk_nodes = [
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_0/outdata/out_1WF",
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_1/outdata/out_1WF",
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_2/outdata/out_1WF",
    ]
    nscf_node = "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_0/task_1/outdata/out_WFK"

    optic_task = abilab.OpticTask(optic_input, nscf_node=nscf_node, ddk_nodes=ddk_nodes)
    flow.register_task(optic_task)
                                                                                                                          
    return flow.allocate()


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    #flow = optic_flow_from_files(options)
    #print("optic manager after allocate", flow[2][0].manager)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
