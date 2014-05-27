#!/usr/bin/env python
from __future__ import division, print_function

import abipy.abilab as abilab
import abipy.data as data  


def make_inputs():
    # Crystalline silicon
    # Calculation of the GW correction to the direct band gap in Gamma
    # Dataset 1: ground state calculation 
    # Dataset 2: NSCF calculation 
    # Dataset 3: calculation of the screening 
    # Dataset 4: calculation of the Self-Energy matrix elements (GW corrections)
    structure = abilab.Structure.from_file(data.cif_file("si.cif"))
    pseudos = data.pseudos("14si.pspnc")

    ecut = 6

    global_vars = dict(
        ecut=ecut,
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


def g0w0_benchmark():
    """
    Build an `AbinitWorkflow` used for benchmarking ABINIT.
    """
    #max_ncpus = manager.policy.max_ncpus
    manager = abilab.TaskManager.from_user_config()
    flow = abilab.AbinitFlow(workdir="g0w0_benchmark", manager=manager)

    gs, nscf, scr, sigma = make_inputs()

    bands = abilab.BandStructureWorkflow(gs, nscf)
    flow.register_work(bands)

    ncpu_list = [1, 2]
    scr_work = abilab.Workflow()

    for ncpu in ncpu_list:
        manager.set_mpi_ncpus(ncpu)
        scr_work.register(scr, manager=manager, deps={bands.nscf_task: "WFK"})
    flow.register_work(scr_work)

    return flow.allocate()


def main():
    flow = g0w0_benchmark()
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    import sys
    sys.exit(main())
