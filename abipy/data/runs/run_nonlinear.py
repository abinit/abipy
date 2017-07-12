#!/usr/bin/env python
"""Flow to compute non-linear optical properties with optic."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import abipy.flowtk as flowtk
import abipy.data as abidata

from abipy import abilab


def make_scf_input(ecut=10, ngkpt=(8, 8, 8)):
    """
    This function constructs an `AbinitInput` for performing a
    GS-SCF calculation in crystalline AlAs.

    Args:
        ecut: cutoff energy in Ha.
        ngkpt: 3 integers specifying the k-mesh for the electrons.

    Return:
        `AbinitInput` object
    """
    # Initialize the AlAs structure from an internal database. Use the pseudos shipped with AbiPy.
    gs_inp = abilab.AbinitInput(structure=abidata.structure_from_ucell("AlAs"),
                                pseudos=abidata.pseudos("13al.981214.fhi", "33as.pspnc"))

    # Set the value of the Abinit variables needed for GS runs.
    gs_inp.set_vars(
        nband=4,
        ecut=ecut,
        ngkpt=ngkpt,
        nshiftk=1,
        shiftk=[0.0, 0.0, 0.0],
        ixc=7,
        nstep=500,
        iscf=7,
        diemac=5.0,
        toldfe=1.0e-22,
        nbdbuf=0,
        kptopt=1,
    )

    gs_inp.set_mnemonics(True)
    return gs_inp


def build_flow(options):
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    scf_input = make_scf_input(ecut=10, ngkpt=(6, 6, 6))
    flow = flowtk.NonLinearCoeffFlow.from_scf_input(workdir, scf_input)
    flow.minimum_abinit_version = "8.5.2"

    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
