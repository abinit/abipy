#!/usr/bin/env python
from __future__ import division, print_function

import abipy.abilab as abilab
import abipy.data as data


def make_template():
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
                       fftalg=112,
                       #
                       istwfk="*1",
                       timopt=-1,
                       chksymbreak=0,
                       prtwf=0,
                       prtden=0,
                       tolvrs=1e-8,
                       nstep=10,
                       )
    inp.set_variables(**global_vars)
    return inp


def build_flow():
    template = make_template()

    flow = abilab.AbinitFlow(workdir="gs_pfft")
    manager = abilab.TaskManager.from_user_config()

    mpi_list = [2]
    fftalg_list = [312, 402, 401]
    ecut_list = list(range(400, 410, 10)) 

    for fftalg in fftalg_list: 
        work = abilab.Workflow()
        for npfft in mpi_list:
            manager.set_mpi_procs(npfft)
            for inp in abilab.input_gen(template, fftalg=fftalg, npfft=npfft, ecut=ecut_list):
                manager.set_autoparal(0)
                #manager.set_mpi_procs(mpi_procs)
                work.register(inp, manager=manager)
        flow.register_work(work)

    return flow.allocate()


def main():
    flow = build_flow()
    flow.build_and_pickle_dump()
    flow.make_scheduler().start()
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
