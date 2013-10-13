#!/usr/bin/env python
from __future__ import division, print_function

import os
import sys
import abipy.data as data  
import abipy.abilab as abilab

from collections import namedtuple
from abipy.data.runs import Tester, decorate_main

def relax_flow():

    pseudo_dirname = data.pseudo_dir
    pp_names = ["14si.pspnc"]

    cif_names = ["si.cif"]
    cif_dirname = data._CIF_DIRPATH

    ksampling_for_cif = {
        "si.cif": dict(ngkpt=[8,8,8], shiftk=[0,0,0]),
        #"foo.cif": dict(ngkpt=[12,12,8], shiftk=[0,0,0]),
    }

    pseudos = [os.path.join(pseudo_dirname, pp_name) for pp_name in pp_names]

    works = []
    for cif, ksampling in ksampling_for_cif.items():

        structure = abilab.Structure.from_file(os.path.join(cif_dirname, cif))

        inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
        inp.set_structure(structure)

        # Global variables
        #inp.set_variables(**conf.vars)

        inp.set_kmesh(**ksampling)

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
            getxred=-1,
            )

        print(inp)

        # Initialize the workflow.
        tester = Tester()
        manager = tester.make_manager()

        workdir = structure.formula + "_natom" + str(len(structure))
        assert workdir not in [work.workdir for work in works]
        work = abilab.Workflow(workdir, manager)

        # Register the input.
        work.register(inp)

        works.append(work)

    return works

@decorate_main
def main():
    flow = relax_flow()
    return flow.build_and_pickle_dump()

if __name__ == "__main__":
    sys.exit(main())

