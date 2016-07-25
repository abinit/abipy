#!/usr/bin/env python
"""Compute phonon frequencies with phonopy (supercells and finite-difference method)."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy.dfpt.abiphonopy import PhonopyWork


def build_flow(options):
    """
    Create a `Flow` for phonon calculations with phonopy:
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")
    #gsinp = abilab.gs_input(structure, pseudos,
    #         kppa=None, ecut=4, pawecutdg=None, scf_nband=6, accuracy="normal", spin_mode="unpolarized",
    #         smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)
    #gsinp.pop("spinat", None)

    gsinp = abilab.AbinitInput(structure, pseudos)
    gsinp.set_vars(ecut=4, nband=5, toldff=1.e-6)
    gsinp.set_autokmesh(nksmall=4)

    flow = abilab.Flow(workdir=workdir)
    work = PhonopyWork.from_gs_input(gsinp, [2,2,2])
    flow.register_work(work)
    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
