#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import abipy.abilab as abilab
import abipy.data as abidata

from abipy.benchmarks import bench_main, BenchmarkFlow


def make_input():
    """Build a template input file for GS calculations with k-point parallelism """
    pseudos = abidata.pseudos("14si.pspnc")
    structure = abidata.structure_from_ucell("Si")

    inp = abilab.AbinitInput(structure, pseudos)
    inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])

    # Global variables
    ecut = 40
    inp.set_vars(
        ecut=ecut,
        nsppol=1,
        nband=40,
        paral_kgb=0,
        #istwfk="*1",
        timopt=-1,
        chksymbreak=0,
        prtwf=0,
        prtden=0,
        tolvrs=1e-10,
        nstep=50,
    )

    return inp


def build_flow(options):
    inp = make_input()

    mpi_range = options.mpi_range
    if mpi_range is None:
        nkpt = len(inp.abiget_ibz().points)
    	mpi_range = range(1, nkpt*inp["nsppol"] + 1) 
    	print("Using mpi_range:", mpi_range, " = nkpt * nsppol")
    else:
    	print("Using mpi_range from cmd line:", mpi_range)

    flow = BenchmarkFlow(workdir="bench_gs_useylm")

    omp_threads = 1
    for useylm in [0, 1]:
        work = abilab.Work()
        for mpi_procs in mpi_range:
            if not options.accept_mpi_omp(mpi_procs, omp_threads): continue
            manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
            work.register_scf_task(inp.new_with_vars(useylm=useylm), manager=manager)
        flow.register_work(work)

    return flow.allocate()


@bench_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
