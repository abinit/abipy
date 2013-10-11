#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import Tester, decorate_main
from pseudo_dojo.dojo.deltaworks import DeltaFactory


def delta_flow():
    # Path of the pseudopotential to test.
    pseudo = data.pseudos("14si.pspnc")[0]

    # Manager used to submit the jobs.
    manager = abilab.TaskManager.simple_mpi(mpi_ncpus=2)

    # Use this for manneback and edit the YAML file according to your platform
    #manager = abilab.TaskManager.from_file("taskmanager.yaml") 

    # Initialize the flow.
    # Don't know why protocol=-1 does not work here.
    flow = abilab.AbinitFlow(workdir="DELTAFACTOR", manager=manager, pickle_protcol=0)

    # Build the workflow for the computation of the deltafactor.
    # The calculation is done with the paramenters and the cif files
    # used in the original paper. We only have to specify 
    # the cutoff energy ecut (Ha) for the pseudopotential.
    factory = DeltaFactory()

    kppa = 6750  # Use this to have the official k-point sampling
    kppa = 50    # this value is for testing purpose.

    work = factory.work_for_pseudo(pseudo, accuracy="normal", kppa=kppa, 
                                   ecut=8, toldfe=1.e-8, smearing="fermi_dirac:0.0005")

    # Register the workflow.
    flow.register_work(work)
    return flow.allocate()


@decorate_main
def main():
    flow = delta_flow()

    flow.build_and_pickle_dump()

if __name__ == "__main__":
    import sys
    sys.exit(main())
