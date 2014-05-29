#!/usr/bin/env python
from __future__ import division, print_function

import abipy.abilab as abilab
import abipy.data as data


def make_input():
    """Build a template input file for GS calculations with paral_kgb"""
    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))
    inp.set_structure(data.structure_from_ucell("Si"))

    # GS run with paral_kgb
    inp.set_kmesh(ngkpt=[1,1,1], shiftk=[0,0,0])

    # Global variables
    global_vars = dict(ecut=20,
                       nsppol=1,
                       nband=20,
                       paral_kgb=1,
                       npkpt=1,
                       npband=1,
                       npfft=1,
                       #
                       istwfk="*1",
                       timopt=-1,
                       chksymbreak=0,
                       prtwf=0,
                       prtden=0,
                       tolvrs=1e-8,
                       nstep=20,
                       )
    inp.set_variables(**global_vars)
    return inp


def build_flow():
    inp = make_input()

    manager = abilab.TaskManager.from_user_config()
    manager.set_autoparal(0)
    flow = abilab.AbinitFlow(workdir="paralkgb_benchmark", manager=manager)

    ncpu_list = [1, 2] #, 4, 8]

    work = abilab.Workflow()
    for ncpus in ncpu_list:
        #manager.set_max_ncpus(ncpus)
        manager.set_mpi_ncpus(ncpus)
        work.register(inp, manager=manager)
    flow.register_work(work)

    return flow.allocate()


def main():
    flow = build_flow()
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    import sys
    sys.exit(main())
