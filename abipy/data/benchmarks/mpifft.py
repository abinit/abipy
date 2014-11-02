#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import abipy.abilab as abilab
import abipy.data as abidata

from abipy.data.benchmarks import bench_main


def make_template(paral_kgb=1, paw=False):
    """Build a template input file for GS calculations with paral_kgb"""
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    inp = abilab.AbiInput(pseudos=pseudos, ndtset=1)
    inp.set_structure(abidata.structure_from_ucell("Si"))

    # GS run with paral_kgb
    inp.set_kmesh(ngkpt=[1,1,1], shiftk=[0,0,0])

    # Global variables
    ecut = 20
    global_vars = dict(
        ecut=ecut,
        pawecutdg=ecut*4,
        nsppol=1,
        nband=20,
        paral_kgb=paral_kgb,
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


def build_flow(options):
    template = make_template()

    flow = abilab.AbinitFlow(workdir="bench_mpifft")
    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)

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


@bench_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    flow.make_scheduler().start()
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
