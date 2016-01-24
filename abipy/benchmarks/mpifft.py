#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import abipy.abilab as abilab
import abipy.data as abidata

from abipy.benchmarks import bench_main


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
        npkpt=1,
        npband=1,
        npfft=1,
        fftalg=112,
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
    fftalg_list = [312, 402, 401]
    ecut_list = list(range(200, 610, 100)) 
    ecut_list = [400,]

    print("Using mpi_range:", options.mpi_range)
    if options.mpi_range is None:
	raise RuntimeError("This benchmark requires --mpi-range")

    template = make_input()
    flow = abilab.Flow(workdir="bench_mpifft")

    omp_threads = 1
    for fftalg in fftalg_list: 
        work = abilab.Work()
        for npfft in options.mpi_range:
            if not options.accept_mpi_omp(npfft, omp_threads): continue
            manager = options.manager.new_with_fixed_mpi_omp(npfft, omp_threads)
            for inp in abilab.input_gen(template, fftalg=fftalg, npfft=npfft, ecut=ecut_list):
                work.register_scf_task(inp, manager=manager)
        flow.register_work(work)

    return flow.allocate()


@bench_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":

    sys.exit(main())
