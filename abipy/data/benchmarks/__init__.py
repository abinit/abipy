from __future__ import division, print_function

import os
import abipy.abilab as abilab


def build_abinit_benchmark(workdir, template, manager):
    """
    Build an `AbinitWorkflow` used for benchmarking ABINIT.

    Args:
        workdir:
            Working directory. 
        template:
            Initial `AbinitInput` that will be modified by adding the 
            parameters reported by autoparal (if any).
        manager:
            `TaskManager` object.

    Return
        `AbinitWorflow` object.
    """
    max_ncpus = manager.policy.max_ncpus

    # Build a temporary workflow with a shell manager just to run ABINIT to get the autoparal configurations.
    shell_manager = manager.to_shell_manager(mpi_ncpus=1, policy=dict(autoparal=1, max_ncpus=max_ncpus))

    w = abilab.Workflow(workdir=os.path.join(workdir, "autoparal_run"), manager=shell_manager)
    w.register(template)
    w.build()

    # Get the configurations suggested by autoparal.
    auto_confs, optimal = shell_manager.autoparal(w[0])

    if auto_confs is None:
        raise ValueError("autoparal returned None")

    #w.rmtree()

    # Select the configurations with reasonable efficiency.
    #for constraint in constraints:
    #    auto_confs.filter(constraint)
    auto_confs = [c for c in auto_confs if c.efficiency > 0.9]

    if not auto_confs:
        raise ValueError("Received empty set of autoparal configurations.")
    print(auto_confs)

    work = abilab.Workflow(workdir=workdir, manager=manager)

    # Build the final workflow using the autoparal configurations selected above.
    # Each task will have its own manager and an optimized input.
    for i, conf in enumerate(auto_confs):
        new_manager = manager.deepcopy()

        # Change the number of MPI nodes
        new_manager.set_mpi_ncpus(conf.mpi_ncpus)

        # Change the number of OpenMP threads.
        #new_manager.set_omp_ncpus(conf.omp_ncpus)

        #print("new_manager.policy", new_manager.policy)

        # Add the ABINIT input for the parallel execution.
        new_input = template.deepcopy()
        new_input.set_variables(**conf.vars)

        work.register(new_input, manager=new_manager)

    return work
