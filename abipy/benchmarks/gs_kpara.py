#!/usr/bin/env python
"""Benchmark for k-point parallelism."""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_input(paw=False):
    """Build a template for GS calculations with k-point parallelism"""
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    structure = abidata.structure_from_ucell("Si")

    inp = abilab.AbinitInput(structure, pseudos)
    inp.set_kmesh(ngkpt=[12,12,12], shiftk=[0,0,0])

    # Global variables
    ecut = 40
    inp.set_vars(
        ecut=ecut,
        pawecutdg=ecut*4 if paw else None,
        nsppol=2,
        nband=40,
        paral_kgb=0,
        #istwfk="*1",
        timopt=-1,
        chksymbreak=0,
        prtwf=0,
        prtden=0,
        tolvrs=1e-10,
    )

    return inp


def build_flow(options):
    inp = make_input(paw=options.paw)

    mpi_list = options.mpi_list
    if mpi_list is None:
        nkpt = len(inp.abiget_ibz().points)
        nks = nkpt * inp["nsppol"]
        mpi_list = [p for p in range(1, nks+1) if nks % p == 0]
        #print("Using mpi_list:", mpi_list)
    else:
        print("Using mpi_list from cmd line:", mpi_list)

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)
    work = flowtk.Work()

    for mpi_procs, omp_threads in product(mpi_list, options.omp_list):
        if not options.accept_mpi_omp(mpi_procs, omp_threads): continue
        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        work.register_scf_task(inp, manager=manager)

    flow.register_work(work)
    return flow.allocate()


@bench_main
def main(options):
    if options.info:
        # print doc string and exit.
        print(__doc__)
        return

    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
