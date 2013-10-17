#!/usr/bin/env python
from __future__ import division, print_function

import sys
import numpy as np
import os

import abipy.abilab as abilab
import abipy.data as data

from abipy.data.runs import Tester, enable_logging

def relax_flow():

    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=1)
    #manager = abilab.TaskManager.from_user_config()
    workdir = "relax_SiC"
                                                         
    flow = abilab.AbinitFlow(workdir, manager)

    pseudos = data.pseudos("14si.pspnc", "6c.pspnc")
    structure = data.structure_from_ucell("sic")

    global_vars = dict(
        chksymbreak=0,
        ecut = 20
        #istwfk="*1",
    )

    ngkpt = [4,4,4]
    shiftk = [[0.5,0.5,0.5],
              [0.5,0.0,0.0],
              [0.0,0.5,0.0],
              [0.0,0.0,0.5]]


    inp = abilab.AbiInput(pseudos=pseudos, ndtset=2)
    inp.set_structure(structure)
    inp.set_variables(**global_vars)

    relax_inp, nscf_inp = inp[1:]

    relax_inp.set_kmesh(ngkpt=ngkpt, shiftk=shiftk)
    relax_inp.set_variables(
        toldff=1e-6,
        tolmxf=1e-5,
        strfact=100,
        ecutsm=0.5,
        dilatmx=1.15,
        ntime=100,
        ionmov=2,
        optcell=1,
    )

    nscf_inp.set_kpath(ndivsm=20)
    nscf_inp.tolwfr = 1e-22

    relax_inp, nscf_inp = inp.split_datasets()

    # Initialize the workflow.
    relax_task = flow.register_task(relax_inp, task_class=abilab.RelaxTask)

    nscf_task = flow.register_task(nscf_inp, deps={relax_task: "DEN"}, task_class=abilab.NscfTask)

    return flow.allocate()

@enable_logging
def main():
    flow = relax_flow()
    return flow.build_and_pickle_dump()


if __name__=="__main__":
    sys.exit(main())

