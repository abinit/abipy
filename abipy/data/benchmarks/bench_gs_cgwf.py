#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.abilab as abilab
import abipy.data as data  

from abipy.data.benchmarks import build_abinit_benchmark

def make_base_input():
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))
    inp.set_structure(data.structure_from_ucell("si"))

    # Global variables
    global_vars = dict(ecut=12,
                       nsppol=1,
                       nband=20,
                       timopt=-1,
                       paral_kgb=1,
                       chksymbreak=0,
                       prtwf=0
                       #accesswff=3,
                    )

    inp.set_variables(**global_vars)

    # Simple GS run.
    #inp.set_kmesh(ngkpt=[4,4,4], shiftk=[0,0,0])
    inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    #inp.set_kmesh(ngkpt=[12,12,12], shiftk=[0.1,0.2,0.3])
    inp.tolvrs = 1e-8

    return inp

def main():
    max_ncpus = 10
    base_input = make_base_input()

    policy = dict(autoparal=0, max_ncpus=10)

    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1, policy=policy)
    #manager = abilab.TaskManager.from_user_config()

    benchmark = build_abinit_benchmark("gs_cgwf", base_input, manager=manager)

    return benchmark.build_and_pickle_dump()


if __name__ == "__main__":
    import sys
    sys.exit(main())
