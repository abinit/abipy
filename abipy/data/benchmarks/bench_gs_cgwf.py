#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab


def make_base_input():
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))
    inp.set_structure(data.structure_from_ucell("si"))

    # Global variables
    global_vars = dict(ecut=18,
                       nsppol=2,
                       nband=20,
                       timopt=-1,
                       paral_kgb=0,
                       #accesswff=3,
                    )

    inp.set_variables(**global_vars)

    # Simple GS run.
    inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    inp.tolvrs = 1e-8

    return inp


def build_abinit_benchmark(workdir, base_input, max_ncpus, manager):
    """
    Initialize the benchmark

    Args:
        workdir:
            Working directory. 
        base_input:
            Initial `AbinitInput`. It will be modified by adding the 
            parameters reported by autoparal (if any).
        max_ncpus:
            Maximum number of CPUs that can be used to run each task of the workflow.
            This parameter will be used to get the autoparal configurations.
        manager:
            `TaskManager` object.
    """

    # Build a temporary workflow just to run ABINIT to get the autoparal configurations.
    simple_manager = manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=1, max_ncpus=max_ncpus))
    print(simple_manager)

    w = abilab.Workflow(workdir=os.path.join(workdir, "autoparal_run"), manager=simple_manager)
    w.register(base_input)
    w.build()

    # Get the configurations suggested by autoparal.
    auto_confs = simple_manager.autoparal(w[0])
    #w.rmtree()

    if auto_confs is None:
        raise ValueError("autoparal returned None")

    # Select the configurations with reasonable efficiency.
    #for constraint in constraints:
    #    auto_confs.filter(constraint)

    #if not auto_confs:
    #    raise ValueError("Empty set of autoparal configurations.")

    print(auto_confs)

    work = abilab.Workflow(workdir=workdir, manager=simple_manager)

    # Build the final workflow using the autoparal configurations selected above.
    for i, conf in enumerate(auto_confs):
        # Set the number of MPI nodes.
        new_manager = manager.deepcopy()
        print("new_manager:", new_manager.policy)

        new_manager.set_mpi_ncpus(conf.mpi_ncpus)

        #new_manager.set_omp_ncpus(conf.omp_ncpus)
        # Change the number of OpenMP threads.
        #if optimal.omp_ncpus > 1:
        #    self.qadapter.set_omp_ncpus(optimal.omp_ncpus)
        #else:
        #    self.qadapter.disable_omp()
        #print(new_manager)

        new_input = base_input.deepcopy()
        new_input.set_variables(**conf.vars)
        #print(new_input)

        work.register(new_input, manager=new_manager)

    return work


def main():
    base_input = make_base_input()

    max_ncpus = 2
    policy=dict(autoparal=0, max_ncpus=max_ncpus)

    manager = abilab.TaskManager(qtype="slurm",
       qparams=dict(
           ntasks=2,
           #partition="hmem",
           time="0:20:00",
       ),
       #setup="SetEnv intel13_intel",
       modules = ["intel/compilerpro/13.0.1.117", "fftw3/intel/3.3"],
       shell_env=dict(
         PATH=("/home/naps/ygillet/NAPS/src/abinit-7.4.3-public/tmp_intel13/src/98_main/:" +
               "/home/naps/ygillet/NAPS/intel13/bin:$PATH"),
         LD_LIBRARY_PATH="/home/naps/ygillet/NAPS/intel13/lib:$LD_LIBRARY_PATH",
       ),
       mpi_runner="mpirun",
       policy=policy
    )

    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=2, policy=policy)

    benchmark = build_abinit_benchmark("hello", base_input, max_ncpus=max_ncpus, manager=manager)
    print(benchmark)

    benchmark.build_and_pickle_dump()

if __name__ == "__main__":
    import sys
    sys.exit(main())
