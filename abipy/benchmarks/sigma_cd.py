#!/usr/bin/env python
"""Analyze the parallel efficiency of the SIGMA code (one shot G0W0 with contour deformation and gwpara==2)"""
import sys
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from itertools import product
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

    ecut = 14
    multi.set_vars(
        ecut=ecut,
        pawecutdg=ecut*4 if paw else None,
        timopt=-1,
        istwfk="*1",
        paral_kgb=0,
    )

    multi.set_kmesh(
        ngkpt=[6,6,6],
        shiftk=[0.0, 0.0, 0.0],
    )

    gs, nscf, scr, sigma = multi.split_datasets()

    gs.set_vars(tolvrs=1e-6,
                nband=4,)

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    nscf.set_vars(iscf=-2,
                  tolwfr=1e-8,
                  nband=300,
                  nbdbuf=50,
                  )

    # Dataset3: Calculation of the screening.
    scr.set_vars(
        optdriver=3,
        gwcalctyp=9,
        gwpara=2,
        nband=25,
        ecutwfn=ecut,
        symchi=1,
        awtr=2,
        inclvkb=0,
        ecuteps=4.0,
        spmeth=1,        # Use Hilbert transform : Im chi0 --> chi0.
        nomegasf=250,    # Number of points for Imchi0
        nfreqre=50,
        nfreqim=10,
        freqremax="25 eV",
        freqremin="0 eV",
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    sigma.set_vars(
        optdriver=4,
        gwcalctyp=9,
        gwpara=2,
        nband=35,
        ecutwfn=ecut,
        ecuteps=4.0,
        ecutsigx=ecut,
        symsigma=1,
        gw_qprange=1,
    )

    return gs, nscf, scr, sigma


def build_flow(options):
    """
    Build an `AbinitWorkflow` used for benchmarking ABINIT.
    """
    gs_inp, nscf_inp, scr_inp, sigma_inp = make_inputs(paw=options.paw)
    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)

    bands = flowtk.BandStructureWork(gs_inp, nscf_inp)
    flow.register_work(bands)
    flow.exclude_from_benchmark(bands)

    scr_work = flowtk.Work()
    scr_work.register_scr_task(scr_inp, deps={bands.nscf_task: "WFK"})
    flow.register_work(scr_work)
    flow.exclude_from_benchmark(scr_work)

    mpi_list = options.mpi_list

    for nband in [100, 200, 300]:
        sigma_work = flowtk.Work()

        if options.mpi_list is None:
            # Cannot call autoparal here because we need a WFK file.
            print("Using hard coded values for mpi_list")
            mpi_list = [np for np in range(1, nband+1) if abs((nband - 4) % np) < 1]
        if options.verbose: print("Using nband %d and mpi_list: %s" % (nband, mpi_list))

        for mpi_procs, omp_threads in product(mpi_list, options.omp_list):
            if not options.accept_mpi_omp(mpi_procs, omp_threads): continue
            inp = sigma_inp.new_with_vars(nband=nband)
            manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
            sigma_work.register_sigma_task(inp, manager=manager,
                                           deps={bands.nscf_task: "WFK", scr_work[0]: "SCR"})
        flow.register_work(sigma_work)

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
