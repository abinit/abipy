#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab


def make_base_input():
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))
    inp.set_structure(data.structure_from_ucell("si"))

    # Global variables
    global_vars = dict(ecut=12,
                       nsppol=1,
                       nband=20,
                       timopt=-1,
                       paral_kgb=0,
                       #accesswff=3,
                    )

    inp.set_variables(**global_vars)

    # Simple GS run.
    #inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    inp.set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0])
    inp.tolvrs = 1e-8

    return inp


def build_abinit_benchmark(workdir, base_input, manager):
    """
    Build an `AbinitWorkflow` used for benchmarking ABINIT.

    Args:
        workdir:
            Working directory. 
        base_input:
            Initial `AbinitInput` that will be modified by adding the 
            parameters reported by autoparal (if any).
        manager:
            `TaskManager` object.

    Return
        `AbinitWorflow` object.
    """
    max_ncpus = manager.policy.max_ncpus

    # Build a temporary workflow just to run ABINIT to get the autoparal configurations.
    simple_manager = manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=1, max_ncpus=max_ncpus))

    w = abilab.Workflow(workdir=os.path.join(workdir, "autoparal_run"), manager=simple_manager)
    w.register(base_input)
    w.build()

    # Get the configurations suggested by autoparal.
    auto_confs, optimal = simple_manager.autoparal(w[0])

    if auto_confs is None:
        raise ValueError("autoparal returned None")

    w.rmtree()

    # Select the configurations with reasonable efficiency.
    #for constraint in constraints:
    #    auto_confs.filter(constraint)
    auto_confs = [c for c in auto_confs if c.efficiency > 0.9]

    if not auto_confs:
        raise ValueError("Received empty set of autoparal configurations.")
    print(auto_confs)

    work = abilab.Workflow(workdir=workdir, manager=simple_manager)

    # Build the final workflow using the autoparal configurations selected above.
    for i, conf in enumerate(auto_confs):
        new_manager = manager.deepcopy()

        # Change the number of MPI nodes
        new_manager.set_mpi_ncpus(conf.mpi_ncpus)

        # Change the number of OpenMP threads.
        #new_manager.set_omp_ncpus(conf.omp_ncpus)

        #print("new_manager.policy", new_manager.policy)
        #print("new_manager.policy", new_manager.policy)

        # Change the ABINIT input.
        new_input = base_input.deepcopy()
        new_input.set_variables(**conf.vars)

        work.register(new_input, manager=new_manager)

    return work


def main():
    base_input = make_base_input()

    policy=dict(autoparal=0, max_ncpus=2)

    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=policy)
    #manager = abilab.TaskManager.from_file("taskmanager.yaml")

    benchmark = build_abinit_benchmark("hello", base_input, manager=manager)
    print(benchmark)

    benchmark.build_and_pickle_dump()

if __name__ == "__main__":
    import sys
    sys.exit(main())
