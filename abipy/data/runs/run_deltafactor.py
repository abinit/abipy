#!/usr/bin/env python
"""Compute the deltafactor for a given pseudopotential."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import sys
import abipy.data as data
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import DeltaFactory


def build_flow(options):
    # Path of the pseudopotential to test.
    #pseudo = data.pseudo("14si.pspnc")
    pseudo = data.pseudo("Si.GGA_PBE-JTH-paw.xml")

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize the flow.
    flow = abilab.Flow(workdir=workdir, manager=options.manager, remove=options.remove)

    # Build the workflow for the computation of the deltafactor.
    # The calculation is done with the parameters and the cif files
    # used in the original paper. We only have to specify
    # the cutoff energy ecut (Ha) for the pseudopotential.
    # The workflow will produce a pdf file with the equation of state
    # and a file deltafactor.txt with the final results in the
    # outdir directory DELTAFACTOR/work_0/outdir.
    factory = DeltaFactory("PBE")

    kppa = 6750  # Use this to have the official k-point sampling
    kppa = 50    # this value is for testing purpose.

    ecut = 8
    pawecutdg = ecut * 2 if pseudo.ispaw else None

    work = factory.work_for_pseudo(pseudo, accuracy="normal", kppa=kppa,
                                   ecut=ecut, pawecutdg=pawecutdg,
                                   toldfe=1.e-8, smearing="fermi_dirac:0.0005")

    # Register the workflow.
    flow.register_work(work)
    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
