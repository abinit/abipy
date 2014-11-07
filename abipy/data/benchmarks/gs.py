#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import abipy.abilab as abilab
import abipy.data as abidata

from abipy.data.benchmarks import bench_main


def make_input(paral_kgb=1, paw=False):
    """Build a template input file for GS calculations with paral_kgb"""
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    inp = abilab.AbiInput(pseudos=pseudos)

    inp.set_structure(abidata.structure_from_ucell("Si"))
    inp.set_kmesh(ngkpt=[2,2,2], shiftk=[0,0,0])

    # Global variables
    ecut = 10
    inp.set_variables(
        ecut=ecut,
        pawecutdg=ecut*4,
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

    return inp


def build_flow(options):
    inp = make_input(paral_kgb=options.paral_kgb, paw=options.paw)

    flow = abilab.AbinitFlow(workdir="bench_gs")
    work = abilab.Workflow()
    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)
    manager.set_autoparal(0)

    for mpi_procs in options.mpi_range:
	for qad in manager.qads:
		qad.min_cores = 1
		qad.max_cores = mpi_procs
        manager.set_mpi_procs(mpi_procs)
        work.register(inp, manager=manager)

    flow.register_work(work)
    return flow.allocate()


@bench_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    flow.rapidfire()
    #flow.make_scheduler().start()
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
