#!/usr/bin/env python
from __future__ import division, print_function

import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import decorate_main


def relax_flow():
    cif_file = data.cif_file("si.cif")
    structure = abilab.Structure.from_file(cif_file)

    pseudos_list = [
        data.pseudos("14si.pspnc"),
        #data.pseudos("14si.fhi"),
    ]

    global_vars = dict(
        ecut=4,  # Assuming that all the pseudos require the same cutoff.
        ngkpt=[8,8,8], 
        shiftk=[0,0,0],
        nshiftk=1,
    )

    workdir = "IONCELL"
    manager = abilab.TaskManager.simple_mpi()
    #manager = abilab.TaskManager.from_user_config()

    flow = abilab.AbinitFlow(workdir, manager)

    for pseudos in pseudos_list:
        inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
        inp.set_structure(structure)

        # Global variables
        inp.set_variables(**global_vars)

        # Dataset 1 (Atom Relaxation)
        inp[1].set_variables(
            optcell=0,
            ionmov=1,
            tolvrs=1e-6
        )

        # Dataset 2 (Atom + Cell Relaxation)
        inp[2].set_variables(
            optcell=1,
            ionmov=1,
            dilatmax=1.1,
            tolvrs=1e-6,
            #getxred=-1,
            )

        ion_inp, ioncell_inp = inp.split_datasets()
        work = abilab.RelaxWorkflow(ion_inp, ioncell_inp)

        flow.register_work(work)

    return flow.allocate()


@decorate_main
def main():
    flow = relax_flow()
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    import sys
    sys.exit(main())

