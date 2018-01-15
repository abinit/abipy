#!/usr/bin/env python
"""Benchmark for phonon calculation with DFPT."""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_inputs(paw=False):
    """
    This function constructs the input files for the phonon calculation:
    GS input + the input files for the phonon calculation.
    """
    # Crystalline AlAs: computation of the second derivative of the total energy
    structure = abidata.structure_from_ucell("AlAs")
    pseudos = abidata.pseudos("13al.981214.fhi", "33as.pspnc")
    if paw:
        raise NotImplementedError("PAW")

    # List of q-points for the phonon calculation.
    qpoints = [
             0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
             2.50000000E-01,  0.00000000E+00,  0.00000000E+00,
             5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
             2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
             5.00000000E-01,  2.50000000E-01,  0.00000000E+00,
            -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
             5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
            -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
            ]
    import numpy as np
    qpoints = np.reshape(qpoints, (-1,3))

    # Global variables used both for the GS and the DFPT run.
    global_vars = dict(
        nband=12,
        ecut=12.0,
        pawecutdg=24.0 if paw else None,
        ngkpt=[8, 8, 8],
        nshiftk=4,
        shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
                0.0, 0.5, 0.0,
                0.5, 0.0, 0.0,
                0.5, 0.5, 0.5],
        #shiftk=[0, 0, 0],
        paral_kgb=0,
        diemac=9.0,
    )

    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)
    gs_inp.set_vars(global_vars, tolvrs=1.0e-18)

    # Get the qpoints in the IBZ. Note that here we use a q-mesh with ngkpt=(4,4,4) and shiftk=(0,0,0)
    # i.e. the same parameters used for the k-mesh in gs_inp.
    #qpoints = gs_inp.abiget_ibz(ngkpt=(4,4,4), shiftk=(0,0,0), kptopt=1).points
    #print("get_ibz", qpoints)

    ph_inp = abilab.AbinitInput(structure, pseudos)

    # Response-function calculation for phonons.
    qpt = [0, 0, 0]
    qpt = [0.25, 0, 0]
    ph_inp.set_vars(
        global_vars,
        rfphon=1,        # Will consider phonon-type perturbation
        nqpt=1,          # One wavevector is to be considered
        qpt=qpt,         # This wavevector is q=0 (Gamma)
        toldfe=1.0e10,
        kptopt=3,
        timopt=-1,
    )
        #rfatpol   1 1   # Only the first atom is displaced
        #rfdir   1 0 0   # Along the first reduced coordinate axis
        #kptopt   2      # Automatic generation of k points, taking

    return gs_inp, ph_inp


def build_flow(options):
    gs_inp, ph_inp = make_inputs()

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)
    gs_work = flowtk.Work()
    gs_work.register_scf_task(gs_inp)
    flow.register_work(gs_work)
    flow.exclude_from_benchmark(gs_work)

    # Get the list of possible parallel configurations from abinit autoparal.
    max_ncpus, min_eff = options.max_ncpus, options.min_eff
    print("Getting all autoparal confs up to max_ncpus:", max_ncpus, "with efficiency >=", min_eff)

    pconfs = ph_inp.abiget_autoparal_pconfs(max_ncpus, autoparal=1)
    if options.verbose: print(pconfs)

    work = flowtk.Work()
    for conf, omp_threads in product(pconfs, options.omp_list):
        mpi_procs = conf.mpi_ncpus
        if not options.accept_conf(conf, omp_threads): continue

        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        inp = ph_inp.new_with_vars(conf.vars)
        work.register_phonon_task(inp, manager=manager, deps={gs_work[0]: "WFK"})

    print("Found %d configurations" % len(work))
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
