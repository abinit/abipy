#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import Tester, decorate_main

def optics_flow():
    structure = data.structure_from_ucell("gaas")

    inp = abilab.AbiInput(pseudos=data.pseudos("31ga.pspnc", "33as.pspnc"), ndtset=5)

    inp.set_structure(structure)

    # Global variables
    kmesh = dict(ngkpt=[4, 4,  4], 
                 nshiftk=4,
                 shiftk=[[0.5, 0.5, 0.5],
                         [0.5, 0.0, 0.0],
                         [0.0, 0.5, 0.0],
                         [0.0, 0.0, 0.5]]
                )

    global_vars = dict(ecut=2,
                       #nband=8,
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
      #getwfk=2,  
      #getden=1,
    )

    #Fourth dataset : ddk response function along axis 1
    #Fifth dataset : ddk response function along axis 2
    #Sixth dataset : ddk response function along axis 3
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
         #getwfk=3,
         tolwfr=1.e-9,
        )

    print(inp)

    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=2)

    scf_inp, nscf_inp, ddk1, ddk2, ddk3 = inp.split_datasets()

    # Initialize the flow.
    flow = abilab.AbinitFlow(workdir="OPTICS", manager=manager)

    bands_work = abilab.BandStructureWorkflow(scf_inp, nscf_inp)
    flow.register_work(bands_work)

    ddk_work = abilab.Workflow()
    for inp in [ddk1, ddk2, ddk3]:
        ddk_work.register(inp, deps={bands_work.nscf_task: "WFK"}, task_class=abilab.DDK_Task)

    flow.register_work(ddk_work)

    # TODO
    # Check is the order of the 1WF files is relevant. Can we use DDK files ordered 
    # in an arbitrary way or do we have to pass (x,y,z)?

    optics_input = """\
/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTICS/work_1/task_0/outdata/out_1WF
/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTICS/work_1/task_1/outdata/out_1WF
/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTICS/work_1/task_2/outdata/out_1WF
/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTICS/work_0/task_1/outdata/out_WFK
0.002         ! Value of the smearing factor, in Hartree
0.0003  0.3   ! Difference between frequency values (in Hartree), and maximum frequency ( 1 Ha is about 27.211 eV)
0.000         ! Scissor shift if needed, in Hartree
0.002         ! Tolerance on closeness of singularities (in Hartree)
1             ! Number of components of linear optic tensor to be computed
11            ! Linear coefficients to be computed (x=1, y=2, z=3)
2             ! Number of components of nonlinear optic tensor to be computed
123 222       ! Non-linear coefficients to be computed
"""
    deps = {task: "1WF" for task in ddk_work}
    deps.update({bands_work.nscf_task: "WFK"})
    #deps = {task: "1WF" for task in ddk_work}.update(

    shell_manager = manager.to_shell_manager()

    optics_task = flow.register_task(optics_input, deps=deps, manager=shell_manager, task_class=abilab.OpticsTask)
    #optics_task.add_deps(deps)
    #flow.register_task(optics_task, deps=deps)

    flow.allocate()
    #print(optics_task)
    #print(optics_task.__class__)
    print("input",optics_task.make_input())
    return flow


@decorate_main
def main():
    flow = optics_flow()
    flow.build_and_pickle_dump()
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
