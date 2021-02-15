#!/usr/bin/env python
"""
Benchmark paral_kgb=1 algorithm with wfoptalg in [default, 1].
default corresponds to the LOBPCG algorithm, 1 enables the Chebyschev solver.
"""
import sys
import operator
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from itertools import product
from functools import reduce
from abipy.benchmarks import bench_main, BenchmarkFlow


def make_input(paw=False):
    """
    Build and return an input file for GS calculations with paral_kgb=1
    """
    pseudos = abidata.pseudos("14si.pspnc", "8o.pspnc")
              #if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    structure = abidata.structure_from_ucell("SiO2-alpha")

    inp = abilab.AbinitInput(structure, pseudos)
    inp.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0,0,0.5])

    # Global variables
    ecut = 40
    inp.set_vars(
        ecut=ecut,
        pawecutdg=ecut*4 if paw else None,
        nsppol=1,
        nband=32,
        paral_kgb=1,
        kptopt=3,
        istwfk="*1",
        timopt=-1,
        chksymbreak=0,
        nsym=1,
        prtwf=0,
        prtden=0,
        tolvrs=1e-8,
        nstep=20,
    )

    return inp


def build_flow(options):
    template = make_input()

    # Get the list of possible parallel configurations from abinit autoparal.
    #max_ncpus, min_eff = options.max_ncpus, options.min_eff
    #print("Getting all autoparal configurations up to max_ncpus: ",max_ncpus," with efficiency >= ",min_eff)
    #pconfs = template.abiget_autoparal_pconfs(max_ncpus, autoparal=1, verbose=options.verbose)
    #if options.verbose: print(pconfs)

    # Processor distribution.
    pconfs = [
       dict(npkpt=64, npband=1, npfft=2), # 128
       dict(npkpt=64, npband=2, npfft=2), # 256
       dict(npkpt=64, npband=2, npfft=4), # 512
       dict(npkpt=64, npband=4, npfft=4), # 1024
    ]

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)

    for wfoptalg in [None, 1]:
        work = flowtk.Work()
        for conf, omp_threads in product(pconfs, options.omp_list):
            #if not options.accept_conf(conf, omp_threads): continue
            mpi_procs = omp_threads * reduce(operator.mul, conf.values(), 1)

            manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
            inp = template.new_with_vars(conf, wfoptalg=wfoptalg)
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
