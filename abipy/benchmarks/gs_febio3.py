#!/usr/bin/env python
"""GS+NSCF calculation for FeBiO3"""

import sys
from itertools import product

import abipy.data as abidata
from abipy import abilab, flowtk
from abipy.benchmarks import BenchmarkFlow, bench_main
from abipy.flowtk import ParalHints

unit_cell = dict(
    acell=3 * [1.0385008112e01],
    natom=10,
    ntypat=3,
    typat=[1, 1, 2, 2, 3, 3, 3, 3, 3, 3],
    znucl=[26, 83, 8],
    xred=[
        [1.4071110772e-03, 1.4071110772e-03, 1.4071110772e-03],
        [5.0140711108e-01, 5.0140711108e-01, 5.0140711108e-01],
        [2.7037366934e-01, 2.7037366934e-01, 2.7037366934e-01],
        [7.7037366934e-01, 7.7037366934e-01, 7.7037366934e-01],
        [3.1190076695e-01, 1.6941360694e-01, 7.1602835037e-01],
        [1.6941360694e-01, 7.1602835037e-01, 3.1190076695e-01],
        [7.1602835037e-01, 3.1190076695e-01, 1.6941360694e-01],
        [2.1602835037e-01, 6.6941360694e-01, 8.1190076695e-01],
        [8.1190076695e-01, 2.1602835037e-01, 6.6941360694e-01],
        [6.6941360694e-01, 8.1190076695e-01, 2.1602835037e-01],
    ],
    rprim=[
        [5.7864276271e-01, 0.0000000000e00, 8.1558111379e-01],
        [-2.8932138135e-01, 5.0111933222e-01, 8.1558111379e-01],
        [-2.8932138135e-01, -5.0111933222e-01, 8.1558111379e-01],
    ],
)

global_vars = dict(
    paral_kgb=0,
    ecut=50,
    nstep=500,
    spinat=[
        [0.0000000000e00, 0.0000000000e00, 3.5716762600e00],
        [0.0000000000e00, 0.0000000000e00, -3.5716762600e00],
        [0.0000000000e00, 0.0000000000e00, 0.0000000000e00],
        [0.0000000000e00, 0.0000000000e00, 0.0000000000e00],
        [0.0000000000e00, 0.0000000000e00, 0.0000000000e00],
        [0.0000000000e00, 0.0000000000e00, 0.0000000000e00],
        [0.0000000000e00, 0.0000000000e00, 0.0000000000e00],
        [0.0000000000e00, 0.0000000000e00, 0.0000000000e00],
        [0.0000000000e00, 0.0000000000e00, 0.0000000000e00],
        [0.0000000000e00, 0.0000000000e00, 0.0000000000e00],
    ],
    nsppol=2,
    nspden=2,
    diemac=5,
    diemix=0.6,
    ixc=7,
    chksymbreak=0,
    dilatmx=1.05,
    ecutsm=0.5,
    nband=60,
    nbdbuf=5,
    ngkpt=[8, 8, 8],
    shiftk=[0.0, 0.0, 0.0],
    nsym=1,
    # iomode=3
)


def make_inputs(options):
    structure = abilab.Structure.from_abivars(unit_cell)

    if options.paw:
        raise RuntimeError("PAW is not implemented")
    pseudos = abidata.pseudos("26fe.pspnc", "83-Bi.GGA.fhi", "8o.pspnc")
    # pseudos = ["fe.pot", "bi.pot", 'o.pot']

    gs_inp = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=2)
    gs_inp.set_vars(global_vars)

    gs_inp.set_vars(timopt=-1, kptopt=3, mem_test=0)

    gs_inp[0].set_vars(tolvrs=1.0e-18)

    gs_inp[1].set_vars(
        iscf=-2,
        tolwfr=1.0e-22,
        nband=100,
        nbdbuf=10,
    )

    return gs_inp.split_datasets()


def build_flow(options):
    gs_inp, nscf_inp = make_inputs(options)

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)

    mpi_list = options.mpi_list
    if mpi_list is None:
        # Get the list of possible parallel configurations from abinit autoparal.
        max_ncpus, min_eff = options.max_ncpus, options.min_eff
        print("Getting all autoparal confs up to max_ncpus:", max_ncpus, "with efficiency >=", min_eff)

        pconfs = gs_inp.abiget_autoparal_pconfs(max_ncpus, autoparal=1)

    else:
        print("Initializing autoparal from command line options")
        pconfs = ParalHints.from_mpi_omp_lists(mpi_list, options.omp_list)
        if options.verbose:
            print(pconfs)

    work = flowtk.Work()
    for conf, omp_threads in product(pconfs, options.omp_list):
        mpi_procs = conf.mpi_ncpus
        # if not options.accept_mpi_omp(mpi_procs,omp_threads): continue
        if not options.accept_conf(conf, omp_threads):
            continue

        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        inp = gs_inp.new_with_vars(conf.vars)
        scf_task = work.register_scf_task(inp, manager=manager)

        inp2 = nscf_inp.new_with_vars(conf.vars)
        work.register_nscf_task(inp2, manager=manager, deps={scf_task: "DEN"})

    print("Found %d configurations" % len(work))
    flow.register_work(work)

    return flow.allocate()


@bench_main
def main(options):
    if options.info:
        # print doc string and exit.
        print(__doc__)
        return None
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
