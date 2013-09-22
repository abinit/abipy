#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab


def make_input():
    structure = data.structure_from_ucell("si")

    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))
    inp.set_structure(structure)

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


class AbinitBenchmark(object):

    def __init__(self, workdir, base_input, max_ncpus, manager):
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
        self.workdir = os.path.abspath(workdir)
        self.base_input = base_input
        self.manager = manager
        self.max_ncpus = max_ncpus

        # Build a temporary workflow just to run ABINIT to get the autoparal configurations.
        simple_manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=dict(autoparal=1, max_ncpus=max_ncpus))
        #print(simple_manager)

        w = abilab.Workflow(workdir=os.path.join(self.workdir, "autoparal_run"), manager=simple_manager)
        w.register(self.base_input)

        # Get the configurations suggested by autoparal.
        self.confs = simple_manager.autoparal(w[0])
        w.rmtree()

        if self.confs is None:
            raise ValueError("autoparal returned None")
        print(self.confs)

        # Select the configurations with reasonable efficiency.
        #self.confs = ["ciao", "bello"] #confs

        # Build the final workflow using the autoparal configurations selected above.
        self.work = work = abilab.Workflow(self.workdir, simple_manager)

        for i, conf in enumerate(self.confs):
            # Set the number of MPI nodes.
            new_manager = self.manager.deepcopy()
            print("new_manager",new_manager.policy)

            new_manager.set_mpi_ncpus(conf.mpi_ncpus)

            #new_manager.set_omp_ncpus(conf.omp_ncpus)
            # Change the number of OpenMP threads.
            #if optimal.omp_ncpus > 1:
            #    self.qadapter.set_omp_ncpus(optimal.omp_ncpus)
            #else:
            #    self.qadapter.disable_omp()
            #print(new_manager)

            new_input = self.base_input.deepcopy()
            new_input.set_variables(**conf.vars)
            #print(new_input)

            work.register(new_input, manager=new_manager)

    def build_and_pickle_dump(self):
       self.work.build_and_pickle_dump()

    #def __str__(self):


def main():
    base_input = make_input()

    max_ncpus = 2
    policy=dict(autoparal=0, max_ncpus=max_ncpus)

    manager = abilab.TaskManager(qtype="slurm",
       qparams=dict(
           ntasks=2,
           #partition="hmem",
           time="0:20:00",
           #account='nobody@nowhere.org',
           #ntasks_per_node=None,
           #cpus_per_task=None,
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

    benchmark = AbinitBenchmark("hello", base_input, max_ncpus=max_ncpus, manager=manager)
    print(benchmark)

    benchmark.build_and_pickle_dump()

if __name__ == "__main__":
    import sys
    sys.exit(main())
