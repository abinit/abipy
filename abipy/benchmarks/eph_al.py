#!/usr/bin/env python
"""Benchmark for electron-phonon calculations."""
import sys
import numpy as np
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_flow_ephinp(options):
    # Preparatory run for E-PH calculations.
    # The sequence of datasets makes the ground states and
    # all of the independent perturbations of the single Al atom
    # for the irreducible qpoints in a 4x4x4 grid.
    # Note that the q-point grid must be a sub-grid of the k-point grid (here 8x8x8)
    pseudos = abidata.pseudos("Al.oncvpsp") if not options.paw else \
              abidata.pseudos("Al.GGA_PBE-JTH-paw.xml")

    structure = abilab.Structure.from_abivars(
        acell=3*[7.5],
        rprim=[0.0, 0.5, 0.5,
               0.5, 0.0, 0.5,
               0.5, 0.5, 0.0],
        typat=1,
        xred=[0.0, 0.0, 0.0],
        ntypat=1,
        znucl=13,
    )

    gs_inp = abilab.AbinitInput(structure, pseudos)

    gs_inp.set_vars(
        nsppol=1,
        prtpot=1,
        istwfk="*1",
        ecut=12.0,
        nband=5,
        occopt=7,    # include metallic occupation function with a small smearing
        tsmear=0.04,
        tolvrs=1e-7,
        timopt=-1,
    )

    # The kpoint grid is minimalistic to keep the calculation manageable.
    gs_inp.set_kmesh(
        ngkpt=[8, 8, 8],
        kptopt=3,
        shiftk=[0.0, 0.0, 0.0],
    )

    # Phonon calculation with 4x4x4
    qpoints = np.reshape([
         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
         2.50000000e-01,  0.00000000e+00,  0.00000000e+00,
         5.00000000e-01,  0.00000000e+00,  0.00000000e+00,
         2.50000000e-01,  2.50000000e-01,  0.00000000e+00,
         5.00000000e-01,  2.50000000e-01,  0.00000000e+00,
        -2.50000000e-01,  2.50000000e-01,  0.00000000e+00,
         5.00000000e-01,  5.00000000e-01,  0.00000000e+00,
        -2.50000000e-01,  5.00000000e-01,  2.50000000e-01,
        ], (-1,3))

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)
    work0 = flow.register_task(gs_inp, task_class=flowtk.ScfTask)
    flow.exclude_from_benchmark(work0)

    ph_work = flowtk.PhononWork.from_scf_task(work0[0], qpoints)
    flow.register_work(ph_work)
    flow.exclude_from_benchmark(ph_work)

    # Build input file for E-PH run.
    eph_inp = gs_inp.new_with_vars(
        optdriver=7,
        #ddb_ngqpt=[1, 1, 1],  # q-mesh used to produce the DDB file (must be consisten with DDB data)
        ddb_ngqpt=[4, 4, 4],   # q-mesh used to produce the DDB file (must be consisten with DDB data)
        eph_intmeth=2,         # Tetra
        eph_fsewin="0.8 eV",   # Energy window around Ef
        eph_mustar=0.12,       # mustar parameter
        # q-path for phonons and phonon linewidths.
        ph_ndivsm=20,
        ph_nqpath=3,
        ph_qpath=[
          0  , 0  , 0,
          0.5, 0  , 0,
          0.5, 0.5, 0,],
        # phonon DOS obtained via Fourier interpolation
        ph_intmeth=2,            # Tetra for phonon DOS and A2F
        ph_smear="0.001 eV",
        ph_wstep="0.0001 eV",
        ph_ngqpt=[16, 16, 16],   # q-mesh for Fourier interpolatation of IFC and a2F(w)
        ph_nqshift=1,
        ph_qshift=[0, 0, 0],
    )

    return flow, eph_inp


def build_flow(options):
    flow, eph_inp = make_flow_ephinp(options)

    mpi_list = options.mpi_list
    if mpi_list is None:
        nkpt = len(eph_inp.abiget_ibz().points)
        nks = nkpt * eph_inp["nsppol"]
        mpi_list = [p for p in range(1, nks+1) if nks % p == 0]
        if options.verbose: print("Using mpi_list:", mpi_list)
    else:
        print("Using mpi_list from cmd line:", mpi_list)

    eph_work = flowtk.Work()
    for mpi_procs, omp_threads in product(mpi_list, options.omp_list):
        if not options.accept_mpi_omp(mpi_procs, omp_threads): continue
        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        eph_work.register_eph_task(eph_inp, manager=manager, deps={flow[0][0]: "WFK", flow[1]: ["DDB", "DVDB"]})

    flow.register_work(eph_work)
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
