#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import abipy.abilab as abilab
import abipy.data as abidata  

from abipy.benchmarks import bench_main, BenchmarkFlow


def make_inputs(paw=False):
    # Crystalline silicon
    # Calculation of the GW correction to the direct band gap in Gamma
    # Dataset 1: ground state calculation 
    # Dataset 2: NSCF calculation 
    # Dataset 3: calculation of the screening 
    # Dataset 4: calculation of the Self-Energy matrix elements (GW corrections)
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")

    multi = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=4)

    ecut = 8
    multi.set_vars(
        ecut=ecut,
        pawecutdg=ecut*4,
        timopt=-1,
        istwfk="*1",
        paral_kgb=0,
    )

    gs, nscf, scr, sigma = multi.split_datasets()

    # This grid is the most economical, but does not contain the Gamma point.
    ngkpt=[2, 2, 2],
    gs_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap. 
    gw_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.0, 0.0, 0.0,  
                0.0, 0.5, 0.5,  
                0.5, 0.0, 0.5,  
                0.5, 0.5, 0.0]
    )

    # Dataset 1 (GS run)
    gs.set_kmesh(**gs_kmesh)
    gs.set_vars(tolvrs=1e-6,
                nband=8,
                )

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    nscf.set_kmesh(**gw_kmesh)
    nscf.set_vars(iscf=-2,
                  tolwfr=1e-8,
                  nband=300,
                  nbdbuf=50,
                  )

    # Dataset3: Calculation of the screening.
    scr.set_kmesh(**gw_kmesh)
    scr.set_vars(
        optdriver=3,   
        gwpara=2,
        nband=25,
        ecutwfn=ecut,   
        symchi=1,
        awtr=2,
        inclvkb=0,
        ecuteps=4.0,    
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    sigma.set_kmesh(**gw_kmesh)
    sigma.set_vars(
        optdriver=4,
        gwpara=2,
        nband=35,
        ecutwfn=ecut,
        ecuteps=4.0,
        ecutsigx=6.0,
        symsigma=1,
        gw_qprange=1,
    )

    return gs, nscf, scr, sigma


def scr_benchmark(options):
    """
    Build an `AbinitWorkflow` used for benchmarking ABINIT.
    """
    gs_inp, nscf_inp, scr_inp, sigma_inp = make_inputs(paw=options.paw)
    flow = BenchmarkFlow(workdir="bench_scr")

    bands = abilab.BandStructureWork(gs_inp, nscf_inp)
    flow.register_work(bands)
    flow.exclude_from_benchmark(bands)

    print("Using mpi_range:", options.mpi_range)
    if options.mpi_range is None:
	raise RuntimeError("This benchmark requires --mpi-range")

    omp_threads = 1
    for nband in [100, 200, 300]:
	scr_work = abilab.Work()
	for mpi_procs in options.mpi_range:
            if not options.accept_mpi_omp(mpi_procs, omp_threads): continue
	    manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
	    scr_work.register(scr_inp, manager=manager, deps={bands.nscf_task: "WFK"})
	flow.register_work(scr_work)

    return flow.allocate()


@bench_main
def main(options):
    flow = scr_benchmark(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
