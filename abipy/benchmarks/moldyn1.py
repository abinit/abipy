#!/usr/bin/env python
"""
Benckmark for Molecular dynamics (based on tutoparal/Input/tmoldyn[1-3].in)
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import numpy as np
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_input(paw=False):
    """
    Build and return an input file for MD calculations with paral_kgb=1
    """
    pseudos = abidata.pseudos("Al.oncvpsp") if not paw else \
              abidata.pseudos("Al.GGA_PBE-JTH-paw.xml")

    # Atomic Positions.
    xred = np.fromstring("""
0.00000000E+00  0.00000000E+00  0.00000000E+00
0.00000000E+00  0.25000000E+00  0.25000000E+00
0.25000000E+00  0.00000000E+00  0.25000000E+00
0.25000000E+00  0.25000000E+00  0.00000000E+00
0.00000000E+00  0.00000000E+00  0.50000000E+00
0.00000000E+00  0.25000000E+00  0.75000000E+00
0.25000000E+00  0.00000000E+00  0.75000000E+00
0.25000000E+00  0.25000000E+00  0.50000000E+00
0.00000000E+00  0.50000000E+00  0.00000000E+00
0.00000000E+00  0.75000000E+00  0.25000000E+00
0.25000000E+00  0.50000000E+00  0.25000000E+00
0.25000000E+00  0.75000000E+00  0.00000000E+00
0.00000000E+00  0.50000000E+00  0.50000000E+00
0.00000000E+00  0.75000000E+00  0.75000000E+00
0.25000000E+00  0.50000000E+00  0.75000000E+00
0.25000000E+00  0.75000000E+00  0.50000000E+00
0.50000000E+00  0.00000000E+00  0.00000000E+00
0.50000000E+00  0.25000000E+00  0.25000000E+00
0.75000000E+00  0.00000000E+00  0.25000000E+00
0.75000000E+00  0.25000000E+00  0.00000000E+00
0.50000000E+00  0.00000000E+00  0.50000000E+00
0.50000000E+00  0.25000000E+00  0.75000000E+00
0.75000000E+00  0.00000000E+00  0.75000000E+00
0.75000000E+00  0.25000000E+00  0.50000000E+00
0.50000000E+00  0.50000000E+00  0.00000000E+00
0.50000000E+00  0.75000000E+00  0.25000000E+00
0.75000000E+00  0.50000000E+00  0.25000000E+00
0.75000000E+00  0.75000000E+00  0.00000000E+00
0.50000000E+00  0.50000000E+00  0.50000000E+00
0.50000000E+00  0.75000000E+00  0.75000000E+00
0.75000000E+00  0.50000000E+00  0.75000000E+00
0.75000000E+00  0.75000000E+00  0.50000000E+00
""", sep=" ").reshape((-1,3))

    # Crystal structure (32 Al atoms)
    structure = abilab.Structure.from_abivars(
        acell=3*[12.81],
        rprim=np.eye(3),
        ntypat=1,
        typat=32*[1],
        znucl=13.0,
        xred=xred,
    )

    inp = abilab.AbinitInput(structure, pseudos)

    inp.set_vars(
        # parallelization
        paral_kgb=1,
        #npband=10,
        #npfft=3,
        #npkpt=4,
        #bandpp=2,

        ecut=3.0,
        pawecutdg=6.0 if paw else None,
        nband=80,
        nsppol=1,

        # SCF cycle parameters
        tolvrs=1.e-3,
        nstep=50,

        # K-points and sym
        occopt=3,
        nsym=1,

        # Molecular Dynamics parameters
        ionmov=12,
        ntime=100,
        dtion=100,
        mdtemp=[3000, 3000],
        tsmear=0.009500446,

        # IO
        prtden=0,
        prtwf=0,
        prteig=0,
        timopt=-1,
    )

    inp.set_kmesh(ngkpt=[1,1,1], shiftk=[0,0,0])

#%% nprocs_to_test = 120

    return inp


def build_flow(options):
    template = make_input()
    #template.abivalidate()

    # Get the list of possible parallel configurations from abinit autoparal.
    max_ncpus, min_eff = options.max_ncpus, options.min_eff
    if max_ncpus is None:
        nkpt = len(template.abiget_ibz().points)
        max_ncpus = nkpt * template["nsppol"] * template["nband"] * 4
    print("Getting all autoparal confs up to max_ncpus:", max_ncpus, "with efficiency >=", min_eff)

    pconfs = template.abiget_autoparal_pconfs(max_ncpus, autoparal=1, verbose=options.verbose)
    if options.verbose: print(pconfs)

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)

    work = flowtk.Work()
    for conf, omp_threads in product(pconfs, options.omp_list):
        mpi_procs = conf.mpi_ncpus
        if not options.accept_conf(conf, omp_threads): continue

        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        inp = template.new_with_vars(conf.vars)
        work.register_scf_task(inp, manager=manager)

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
