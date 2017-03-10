#!/usr/bin/env python
"""Band structure of silicon with the HT interface."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata  
import abipy.flowapi as flowapi

from abipy import abilab

def build_flow(options):

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Initialize structure and pseudos.
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    # Initialize the flow.
    flow = flowapi.Flow(workdir=workdir, manager=options.manager, remove=options.remove)

    # Use ebands_input factory function to build inputs.
    multi = abilab.ebands_input(structure, pseudos, kppa=40, nscf_nband=6, ndivsm=10, ecut=6)
    work = flowapi.BandStructureWork(scf_input=multi[0], nscf_input=multi[1])

    flow.register_work(work)
    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
