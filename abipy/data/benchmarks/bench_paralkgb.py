#!/usr/bin/env python
from __future__ import division, print_function

import abipy.abilab as abilab
import abipy.data as abidata


def make_input(paral_kgb=1, paw=False):
    """Build a template input file for GS calculations with paral_kgb"""
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    inp = abilab.AbiInput(pseudos=pseudos)

    inp.set_structure(abidata.structure_from_ucell("Si"))
    inp.set_kmesh(ngkpt=[1,1,1], shiftk=[0,0,0])

    # Global variables
    global_vars = dict(ecut=20,
                       nsppol=1,
                       nband=20,
                       paral_kgb=paral_kgb,
                       npkpt=1,
                       npband=1,
                       npfft=1,
                       #
                       istwfk="*1",
                       timopt=-1,
                       chksymbreak=0,
                       prtwf=0,
                       prtden=0,
                       tolvrs=1e-10,
                       nstep=50,
                       )
    inp.set_variables(**global_vars)
    return inp


def build_flow():
    inp = make_input(paral_kgb=1, paw=False)

    flow = abilab.AbinitFlow(workdir="paralkgb_benchmark")
    work = abilab.Workflow()
    #import yaml
    #d = yaml.load("taskmanager.yml")
    manager = abilab.TaskManager.from_user_config()

    mpi_list = [1, 2] #, 4, 8]
    for mpi_procs in mpi_list:
        manager.set_autoparal(0)
        manager.set_mpi_procs(mpi_procs)
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
