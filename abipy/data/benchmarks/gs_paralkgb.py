#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import abipy.abilab as abilab
import abipy.data as abidata

from abipy.data.benchmarks import bench_main, BenchmarkFlow


def make_input(paw=False):
    """
    Build and return an input file for GS calculations with paral_kgb=1
    """
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    structure = abidata.structure_from_ucell("Si")

    inp = abilab.AbinitInput(structure, pseudos)
    inp.set_kmesh(ngkpt=[1,1,1], shiftk=[0,0,0])

    # Global variables
    ecut = 20
    inp.set_vars(
        ecut=ecut,
        pawecutdg=ecut*4,
        nsppol=1,
        nband=20,
        paral_kgb=1,
        #npkpt=1,
        #npband=1,
        #npfft=1,
        #fftalg=112,
        #
        #istwfk="*1",
        timopt=-1,
        chksymbreak=0,
        prtwf=0,
        prtden=0,
        tolvrs=1e-8,
        nstep=10,
    )

    return inp


def build_flow(options):
    template = make_input()

    # Get the list of possible parallel configurations from abinit autoparal.
    max_ncpus = 10
    pconfs = template.abiget_autoparal_pconfs(max_ncpus, autoparal=1)
    print(pconfs)

    flow = BenchmarkFlow(workdir="bench_paralkgb")
    work = abilab.Work()

    for conf in pconfs:
        if conf.efficiency < 0.7: continue
        inp = template.deepcopy()
        inp.set_vars(conf.vars)

        manager = options.manager.deepcopy()
        manager.policy.autoparal = 0
        #manager.set_autoparal(0)
        manager.set_mpi_procs(conf.mpi_ncpus)

        work.register(inp, manager=manager)

    flow.register_work(work)
    return flow.allocate()


@bench_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
