#!/usr/bin/env python
"""
Flow for quasi-harmonic calculations
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy.flows.qha import QhaFlow


def build_flow(options):
    """
    Create a `QhaFlow` for quasi-harmonic calculations.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize structure and pseudos
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    # Build input for GS calculation.
    gsinp = abilab.AbinitInput(structure, pseudos)
    gsinp.set_vars(ecut=4, nband=4, toldff=1.e-6)
    gsinp.set_autokmesh(nksmall=2)

    volumes = [gsinp.structure.volume]
    flow = QhaFlow.from_gsinp(workdir, gsinp, volumes, ngqpt=[2,2,2])

    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
