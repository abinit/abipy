#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import abipy.abilab as abilab
import abipy.data as abidata

from abipy.benchmarks import bench_main, BenchmarkFlow


def make_input(paw=False):
    """Build a template input file for GS calculations with k-point parallelism """
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    structure = abidata.structure_from_ucell("Si")

    inp = abilab.AbinitInput(structure, pseudos)
    inp.set_kmesh(ngkpt=[1,1,1], shiftk=[0,0,0])

    # Global variables
    ecut = 10
    inp.set_vars(
        ecut=ecut,
        pawecutdg=ecut*4,
        nsppol=1,
        nband=20,
        paral_kgb=0,
        #istwfk="*1",
        #fftalg=312,
        timopt=-1,
        chksymbreak=0,
        prtwf=0,
        prtden=0,
        tolvrs=1e-10,
        nstep=50,
    )

    return inp


def build_flow(options):
    inp = make_input(paw=options.paw)
    nkpt = len(inp.abiget_ibz().points)

    flow = BenchmarkFlow(workdir="bench_gs_pureomp")
    work = abilab.Work()

    omp_range = options.omp_range
    print("Using omp-range:", omp_range)
    if omp_range is None:
        raise RuntimeError("--omp-range must be specified for this benchmark")

    for omp_threads in omp_range:
        manager = options.manager.deepcopy()
        manager.policy.autoparal = 0
        manager.set_mpi_procs(1)
        manager.set_omp_threads(omp_threads)
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
