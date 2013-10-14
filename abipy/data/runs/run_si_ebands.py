#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import Tester, decorate_main

def bands_flow()
    # Create the task defining the calculation and run.
    tester = Tester()
    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1)
    #manager = abilab.TaskManager.from_user_config()

    flow  = abilab.AbinitFlow(workdir=tester.workdir, manager=manager)

    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=2)
    structure = inp.set_structure_from_file(data.cif_file("si.cif"))

    # Global variables
    global_vars = dict(ecut=6,
                       nband=8,
                       timopt=-1,
                       accesswff=3,
                       istwfk="*1",
                    )

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0], dtset=1)
    inp[1].set_variables(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    inp[2].set_kpath(ndivsm=6, kptbounds=kptbounds)
    inp[2].set_variables(tolwfr=1e-12,
                         getden=-1
                        )
    
    # Initialize the workflow.
    gs_inp, nscf_inp = inp.split_datasets()
    work = abilab.BandStructureWorkflow(gs_inp, nscf_inp)

    flow.register_work(work)

    return flow.allocate()

    #tester.set_work_and_run(work)
    #if tester.retcode != 0:
    #    return tester.retcode

    ## Remove all files except those matching these regular expression.
    ##work.rmtree(exclude_wildcard="*.abi|*.abo|*_WFK*|*_GSR.nc|*DEN-etsf.nc")

    #work[0].rename("out_DS1_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
    #work[0].rename("out_DS1_DEN-etsf.nc", "si_DEN-etsf.nc")
    #work[0].rename("out_DS1_GSR.nc", "si_scf_GSR.nc")

    #work[0].rename("out_DS2_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")
    #work[0].rename("out_DS2_GSR.nc", "si_nscf_GSR.nc")

    #work[0].remove_files("out_DS2_DEN-etsf.nc")

    #tester.finalize()
    #return tester.retcode 


@decorate_main
def main():
    flow = bands_flow()
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    import sys
    sys.exit(main())

