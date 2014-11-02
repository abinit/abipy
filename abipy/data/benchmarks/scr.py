#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals, absolute_import

import abipy.abilab as abilab
import abipy.data as abidata  

from abipy.data.benchmarks import bench_main


def make_inputs(paral_kgb=1, paw=False):
    # Crystalline silicon
    # Calculation of the GW correction to the direct band gap in Gamma
    # Dataset 1: ground state calculation 
    # Dataset 2: NSCF calculation 
    # Dataset 3: calculation of the screening 
    # Dataset 4: calculation of the Self-Energy matrix elements (GW corrections)
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc") if not paw else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")

    ecut = 6
    global_vars = dict(
        ecut=ecut,
        pawecutdg=ecut*4,
        gwpara=2,
        timopt=-1,
        istwfk="*1",
    )

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=4)
    inp.set_structure(structure)
    inp.set_variables(**global_vars)

    gs, nscf, scr, sigma = inp.split_datasets()

    # This grid is the most economical, but does not contain the Gamma point.
    gs_kmesh = dict(
        ngkpt=[2, 2, 2],
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap. 
    gw_kmesh = dict(
        ngkpt=[2, 2, 2],
        shiftk=[0.0, 0.0, 0.0,  
                0.0, 0.5, 0.5,  
                0.5, 0.0, 0.5,  
                0.5, 0.5, 0.0]
    )

    # Dataset 1 (GS run)
    gs.set_kmesh(**gs_kmesh)
    gs.set_variables(tolvrs=1e-6,
                     nband=4,
                     )

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    nscf.set_kmesh(**gw_kmesh)
    nscf.set_variables(iscf=-2,
                       tolwfr=1e-12,
                       nband=35,
                       nbdbuf=5,
                       )

    # Dataset3: Calculation of the screening.
    scr.set_kmesh(**gw_kmesh)
    scr.set_variables(
        optdriver=3,   
        nband=25,
        ecutwfn=ecut,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    sigma.set_kmesh(**gw_kmesh)
    sigma.set_variables(
        optdriver=4,
        nband=35,
        ecutwfn=ecut,
        ecuteps=4.0,
        ecutsigx=6.0,
        symsigma=1,
        gw_qprange=1,
    )

    return gs, nscf, scr, sigma


def g0w0_benchmark(options):
    """
    Build an `AbinitWorkflow` used for benchmarking ABINIT.
    """
    gs_inp, nscf_inp, scr_inp, sigma_inp = make_inputs(paw=options.paw, paral_kgb=options.paral_kgb)
    flow = abilab.AbinitFlow(workdir="bench_scr")
    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)

    bands = abilab.BandStructureWorkflow(gs_inp, nscf_inp)
    flow.register_work(bands)

    mpi_list = [1, 2]
    scr_work = abilab.Workflow()

    for mpi_procs in mpi_list:
        manager.set_autoparal(0)
        manager.set_mpi_procs(mpi_procs)
        scr_work.register(scr_inp, manager=manager, deps={bands.nscf_task: "WFK"})
    flow.register_work(scr_work)

    return flow.allocate()


@bench_main
def main(options):
    flow = g0w0_benchmark(options)
    flow.build_and_pickle_dump()
    flow.make_scheduler().start()
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
