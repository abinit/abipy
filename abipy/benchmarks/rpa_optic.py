#!/usr/bin/env python
"""Benchmark for Optic calculations."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.data as data
import abipy.abilab as abilab
import abipy.flowtk as flowtk

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_base_flow(options):
    multi = abilab.MultiDataset(structure=data.structure_from_ucell("GaAs"),
                                pseudos=data.pseudos("31ga.pspnc", "33as.pspnc"), ndtset=5)

    # Global variables
    kmesh = dict(ngkpt=[4, 4, 4],
                 nshiftk=4,
                 shiftk=[[0.5, 0.5, 0.5],
                         [0.5, 0.0, 0.0],
                         [0.0, 0.5, 0.0],
                         [0.0, 0.0, 0.5]]
    )

    paral_kgb = 1
    global_vars = dict(ecut=2, paral_kgb=paral_kgb)
    global_vars.update(kmesh)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_vars(
        tolvrs=1e-6,
        nband=4,
    )

    # NSCF run with large number of bands, and points in the the full BZ
    multi[1].set_vars(
        iscf=-2,
        nband=20,
        nstep=25,
        kptopt=1,
        tolwfr=1.e-9,
        #kptopt=3,
    )

    # Fourth dataset: ddk response function along axis 1
    # Fifth dataset: ddk response function along axis 2
    # Sixth dataset: ddk response function along axis 3
    for dir in range(3):
        rfdir = 3 * [0]
        rfdir[dir] = 1

        multi[2+dir].set_vars(
            iscf=-3,
            nband=20,
            nstep=1,
            nline=0,
            prtwf=3,
            kptopt=3,
            nqpt=1,
            qpt=[0.0, 0.0, 0.0],
            rfdir=rfdir,
            rfelfd=2,
            tolwfr=1.e-9,
        )

    scf_inp, nscf_inp, ddk1, ddk2, ddk3 = multi.split_datasets()

    # Initialize the flow.
    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.manager)

    bands_work = flowtk.BandStructureWork(scf_inp, nscf_inp)
    flow.register_work(bands_work)
    flow.exclude_from_benchmark(bands_work)

    ddk_work = flowtk.Work()
    for inp in [ddk1, ddk2, ddk3]:
        ddk_work.register_ddk_task(inp, deps={bands_work.nscf_task: "WFK"})

    flow.register_work(ddk_work)
    flow.exclude_from_benchmark(ddk_work)

    return flow


def build_flow(options):
    flow = make_base_flow(options)

    optic_input = abilab.OpticInput(
        broadening=0.002,
        domega=0.0003,
        maxomega=0.3,
        scissor=0.000,
        tolerance=0.002,
        num_lin_comp=1,
        lin_comp=11,
        num_nonlin_comp=2,
        nonlin_comp=(123, 222),
    )

    mpi_list = options.mpi_list
    if mpi_list is None:
        mpi_list = [1,2,4,8]
        print("Using mpi_list:", mpi_list)
    else:
        print("Using mpi_list from cmd line:", mpi_list)

    work = flowtk.Work()
    for mpi_procs, omp_threads in product(mpi_list, options.omp_list):
        if not options.accept_mpi_omp(mpi_procs, omp_threads): continue
        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        optic_task = flowtk.OpticTask(optic_input, manager=manager, nscf_node=flow[0].nscf_task, ddk_nodes=flow[1])
        work.register_task(optic_task)

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
