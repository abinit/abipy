#!/usr/bin/env python
"""Benchmark for DDK calculations."""
import sys
import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk

from itertools import product
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_inputs(paw=False):
    """
    This function constructs the input files for the DDK in SiC.
    Return: GS input, NSCF input, DDK input.
    """
    pseudos = abidata.pseudos("14si.pspnc", "6c.pspnc") if not paw else \
              abidata.pseudos("Si.GGA_PBE-JTH-paw.xml", "6c.lda.atompaw")

    structure = abidata.structure_from_ucell("SiC")

    multi = abilab.MultiDataset(structure, pseudos=pseudos, ndtset=3)

    # Global variables.
    global_vars = dict(
        ecut=32,
        pawecutdg=32 if paw else None,
        nband=24,
        nbdbuf=4,
        # Definition of the SCF procedure
        nstep=100,
        diemac=9.0,
    )

    multi.set_vars(global_vars)

    # Definition of the k-point grid
    multi.set_kmesh(
        ngkpt=[16, 16, 16],
        #ngkpt=[2, 2, 2],
        shiftk=[0.0, 0.0, 0.5,
                0.0, 0.5, 0.0,
                0.5, 0.0, 0.0,
                0.5, 0.5, 0.5]
    )

    multi[0].set_vars(
      kptopt=1,          # Automatic generation of k points with symmetries.
      tolvrs=1.0e-6,
      nband=12,
      nbdbuf=2,
    )

    multi[1].set_vars(
      iscf=-2,
      #getwfk=1,
      #getden=1,
      kptopt=3,
      tolwfr=1.0e-6,
    )

    multi[2].set_vars(
      #getwfk=2,
      #getden=1,
      iscf=-3,        # Need this non-self-consistent option for d/dk
      rfelfd=2,       # Calculate d/dk wave function only
      tolwfr=1.0e-18, # Use wave function residual criterion instead
      kptopt=3,
      rfdir=[1, 0, 0],
      nstep=1,
      nline=0,
    )

    gs_inp, nscf_inp, ddk_inp = multi.split_datasets()

    return gs_inp, nscf_inp, ddk_inp


def build_flow(options):
    gs_inp, nscf_inp, ddk_inp = make_inputs()

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)

    ebands_work = flowtk.BandStructureWork(gs_inp, nscf_inp)
    flow.register_work(ebands_work)
    flow.exclude_from_benchmark(ebands_work)

    # Get the list of possible parallel configurations from abinit autoparal.
    max_ncpus, min_eff = options.max_ncpus, options.min_eff
    print("Getting all autoparal confs up to max_ncpus: ",max_ncpus," with efficiency >= ",min_eff)

    pconfs = ddk_inp.abiget_autoparal_pconfs(max_ncpus, autoparal=1)
    if options.verbose: print(pconfs)

    work = flowtk.Work()
    for conf, omp_threads in product(pconfs, options.omp_list):
        mpi_procs = conf.mpi_ncpus
        if not options.accept_conf(conf, omp_threads): continue

        manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
        inp = ddk_inp.new_with_vars(conf.vars)
        work.register_ddk_task(inp, manager=manager, deps={ebands_work[1]: "WFK"})

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
