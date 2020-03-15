#!/usr/bin/env python
"""Benchmark for phonon calculation with DFPT."""
import sys
import numpy as np
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_inputs(paw=False):
    pseudos = abidata.pseudos("56ba.psp_mod", "22ti.psp_mod", "8o.psp_mod")
    #pseudos = abidata.pseudos("56ba.psp_mod", "22ti.psp_mod", "8o.pspnc")
    if paw: raise NotImplementedError("PAW")

    # SLAB ending TiO2 double layer. paralelectric configuration
    # N=9
    # Supercell and atoms

    xcart = np.fromstring("""
0.0000000000E+00  0.0000000000E+00 -4.2633349730E+00
3.7794522658E+00  3.7794522658E+00 -3.2803418097E+00
0.0000000000E+00  3.7794522658E+00 -3.6627278067E+00
3.7794522658E+00  0.0000000000E+00  6.5250113947E-01
0.0000000000E+00  3.7794522658E+00 -1.0555036964E-01
3.7794522658E+00  3.7794522658E+00  3.2682166278E-01
0.0000000000E+00  0.0000000000E+00  3.9815918094E+00
3.7794522658E+00  3.7794522658E+00  4.0167907030E+00
3.7794522658E+00  0.0000000000E+00  7.7541444349E+00
0.0000000000E+00  3.7794522658E+00  7.6664087705E+00
3.7794522658E+00  3.7794522658E+00  7.7182324796E+00
0.0000000000E+00  0.0000000000E+00  1.1412913350E+01
3.7794522658E+00  3.7794522658E+00  1.1416615533E+01
3.7794522658E+00  0.0000000000E+00  1.5117809063E+01
0.0000000000E+00  3.7794522658E+00  1.5117809063E+01
3.7794522658E+00  3.7794522658E+00  1.5117809063E+01
0.0000000000E+00  0.0000000000E+00  1.8822704777E+01
3.7794522658E+00  3.7794522658E+00  1.8819002593E+01
3.7794522658E+00  0.0000000000E+00  2.2481473692E+01
0.0000000000E+00  3.7794522658E+00  2.2569209355E+01
3.7794522658E+00  3.7794522658E+00  2.2517385647E+01
0.0000000000E+00  0.0000000000E+00  2.6254026317E+01
3.7794522658E+00  3.7794522658E+00  2.6218827423E+01
3.7794522658E+00  0.0000000000E+00  2.9583116986E+01
0.0000000000E+00  3.7794522658E+00  3.0341168496E+01
3.7794522658E+00  3.7794522658E+00  2.9908796464E+01
0.0000000000E+00  0.0000000000E+00  3.4498953099E+01
3.7794522658E+00  3.7794522658E+00  3.3515959935E+01
0.0000000000E+00  3.7794522658E+00  3.3898345933E+01
""", sep=" ").reshape((-1,3))

    # Crystal structure.
    structure = abilab.Structure.from_abivars(
        #acell="4.0 4.0 28.0 Angstrom",
        acell=abilab.ArrayWithUnit([4.0, 4.0, 28], "ang").to("bohr"),
        rprim=np.eye(3),
        typat=[int(i) for i in "3 3 2 3 3 2 1 3 3 3 2 1 3 3 3 2 1 3 3 3 2 1 3 3 3 2 3 3 2".split()],
        znucl=[56, 22, 8],
        xcart=xcart,
    )

    multi = abilab.MultiDataset(structure, pseudos, ndtset=2)

    # Global variables used both for the GS and the DFPT run.
    multi.set_vars(
        #electronic structure
        nband=120,
        ecut=15.0,
        pawecutdg=30.0 if paw else None,
        nstep=80,
        ngkpt=[4,4,1],
        shiftk=[0,0,0],
        paral_kgb=1,
        timopt=-1,
        prtden=0,
    )

    gs_inp, ph_inp = multi.split_datasets()

    gs_inp.set_vars(tolvrs=1e-6, kptopt=1)

    ph_inp.set_vars(
        kptopt=3,
        rfphon=1,
        irdwfk=1,
        rfatpol=[1, 1],
        rfdir=[0, 0, 1],
        nqpt=1,
        qpt=[0.0, 0.25, 0.0],
        prtwf=0,
        #tolwfr=1.0e-22,
        toldfe=1.0e-9,
        tolrde=0.0, # This is a development input variable, used in the present case to avoid load unbalance
                    # when studying the scaling with respect to the number of cores. Do not define it in
                    # your production runs
    )

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
    print("Getting all autoparal confs up to max_ncpus:", max_ncpus," with efficiency >=", min_eff)

    pconfs = ph_inp.abiget_autoparal_pconfs(max_ncpus, autoparal=1)
    if options.verbose: print(pconfs)

    work = flowtk.Work()
    for conf, omp_threads in product(pconfs, options.omp_list):
        if not options.accept_conf(conf, omp_threads): continue

        manager = options.manager.new_with_fixed_mpi_omp(conf.mpi_procs, omp_threads)
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
