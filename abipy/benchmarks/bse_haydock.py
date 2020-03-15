#!/usr/bin/env python
"""Analyze the parallel efficiency of the BSE code (Haydock method with model dielectric function)"""
import sys
import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_inputs(paw=False):
    # Crystalline silicon
    # Calculation of the GW correction to the direct band gap in Gamma
    # Dataset 1: ground state calculation
    # Dataset 2: BSE with haydock method and model dielectric function.
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")

    multi = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=2)

    ecut = 8
    multi.set_vars(
        ecut=ecut,
        pawecutdg=ecut*4 if paw else None,
        nsppol=1,
        timopt=-1,
        istwfk="*1",
        paral_kgb=0,
    )

    # This grid is the most economical, but does not contain the Gamma point.
    multi.set_kmesh(
        ngkpt=[8, 8, 8],
        shiftk=[0.0, 0.0, 0.0],
    )

    gs, bse = multi.split_datasets()

    gs.set_vars(tolvrs=1e-6,
                nband=8,
                )

    # Dataset 6 BSE equation with Model dielectric function and Haydock (only resonant + W + v)
    # Note that SCR file is not needed here
    bse.set_vars(
        optdriver=99,
        ecutwfn=ecut,
        ecuteps=4.0,
        inclvkb=2,
        bs_algorithm=2,        # Haydock
        bs_haydock_niter=60,   # No. of iterations for Haydock
        bs_exchange_term=1,
        bs_coulomb_term=21,    # Use model W and full W_GG.
        mdf_epsinf=12.0,
        bs_calctype=1,         # Use KS energies and orbitals to construct L0
        mbpt_sciss="0.8 eV",
        bs_coupling=0,
        bs_loband=2,
        nband=8,
        bs_freq_mesh="0 10 0.1 eV",
        bs_hayd_term=0,        # No terminator
        #gwmem=01               # Compute the model-dielectric function on-the-fly.
    )

    return gs, bse


def build_flow(options):
    """
    Build an `AbinitWorkflow` used for benchmarking ABINIT.
    """
    gs_inp, bse_inp = make_inputs(paw=options.paw)
    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)

    gs_work = flowtk.Work()
    gs_work.register_scf_task(gs_inp)
    flow.register_work(gs_work)
    flow.exclude_from_benchmark(gs_work)

    mpi_list = options.mpi_list

    if options.mpi_list is None:
        nkpt = len(gs_inp.abiget_ibz().points)
        ntrans = (2*2*nkpt)**2
        mpi_list = [p for p in range(1, 1 + ntrans) if ntrans % p == 0]
    if options.verbose: print("Using mpi_list:", mpi_list)

    bse_work = flowtk.Work()
    for mpi_procs, omp_threads in product(mpi_list, options.omp_list):
        if not options.accept_mpi_omp(mpi_procs, omp_threads): continue
        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        bse_work.register_bse_task(bse_inp, manager=manager, deps={gs_work[0]: "WFK"})
    flow.register_work(bse_work)

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
