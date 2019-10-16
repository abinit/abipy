#!/usr/bin/env python
"""
This benchmark compares GS calculations at the Gamma point done with istwfk in [1,2]
"""
import sys
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata

from abipy.benchmarks import bench_main, BenchmarkFlow


def make_input(paw=False):
    """
    Build and return an input file for GS calculations with paral_kgb=1
    """
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    structure = abidata.structure_from_ucell("Si")

    inp = abilab.AbinitInput(structure, pseudos)
    inp.set_kmesh(ngkpt=[1,1,1], shiftk=[0,0,0])

    # Global variables
    ecut = 20
    inp.set_vars(
        ecut=ecut,
        pawecutdg=ecut*4 if paw else None,
        nsppol=1,
        nband=20,
        paral_kgb=1,
        istwfk="*1",
        timopt=-1,
        chksymbreak=0,
        prtwf=0,
        prtden=0,
        tolvrs=1e-8,
        nstep=50,
    )

    return inp


def build_flow(options):
    template = make_input()

    # Get the list of possible parallel configurations from abinit autoparal.
    max_ncpus, min_eff = options.max_ncpus, options.min_eff
    print("Getting all autoparal confs up to max_ncpus: ",max_ncpus," with efficiency >= ",min_eff)

    pconfs = template.abiget_autoparal_pconfs(max_ncpus, autoparal=1, verbose=options.verbose)
    if options.verbose: print(pconfs)

    flow = BenchmarkFlow(workdir=options.get_workdir(__file__), remove=options.remove)

    omp_threads = 1
    for istwfk in [1, 2]:
        work = flowtk.Work()
        for conf in pconfs:
            mpi_procs = conf.mpi_ncpus
            if not options.accept_conf(conf, omp_threads): continue

            manager = options.manager.new_with_fixed_mpi_omp(mpi_procs, omp_threads)
            inp = template.new_with_vars(conf.vars, istwfk=istwfk)
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
