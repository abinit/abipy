#!/usr/bin/env python
"""Compute phonon frequencies with phonopy (supercells and finite-difference method)."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk
from abipy.flowtk.abiphonopy import PhonopyWork


def build_flow(options):
    """
    Create a `Flow` for phonon calculations with phonopy:
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
    # This gives ngkpt = 4x4x4 with 4 shifts for the initial unit cell.
    # The k-point sampling will be rescaled when we build the supercell in PhonopyWork.
    gsinp.set_autokmesh(nksmall=4)
    #gsinp.set_vars(ngkpt=[4, 4, 4])

    flow = flowtk.Flow(workdir=workdir)

    # Use a 2x2x2 supercell to compute phonons with phonopy
    work = PhonopyWork.from_gs_input(gsinp, scdims=[2,2,2])
    flow.register_work(work)

    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
