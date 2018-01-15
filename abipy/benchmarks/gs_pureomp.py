#!/usr/bin/env python
"""Analyze the scalability of the OpenMP sections in the GS part. 1 k-point, cg method."""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from abipy.benchmarks import bench_main, BenchmarkFlow


def make_input(paw=False):
    """Build a template input file for GS calculations with k-point parallelism """
    pseudos = abidata.pseudos("14si.pspnc", "8o.pspnc") if not paw else \
              abidata.pseudos("Si.GGA_PBE-JTH-paw.xml", "o.paw")
    structure = abidata.structure_from_ucell("SiO2-alpha")

    inp = abilab.AbinitInput(structure, pseudos)
    inp.set_kmesh(ngkpt=[1,1,1], shiftk=[0,0,0])

    # Global variables
    ecut = 24
    inp.set_vars(
        ecut=ecut,
        pawecutdg=ecut*2 if paw else None,
        nsppol=1,
        nband=28,
        paral_kgb=0,
        #istwfk="*1",
        #fftalg=312,
        timopt=-1,
        chksymbreak=0,
        prtwf=0,
        prtden=0,
        tolvrs=1e-10,
    )

    return inp


def build_flow(options):
    inp = make_input(paw=options.paw)
    nkpt = len(inp.abiget_ibz().points)

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)
    work = flowtk.Work()

    omp_list = options.omp_list
    if omp_list is None: omp_list = [1, 2, 4, 6]
    print("Using omp_list:", omp_list)

    mpi_procs = 1
    for omp_threads in omp_list:
        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        work.register(inp, manager=manager)

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
