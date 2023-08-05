#!/usr/bin/env python
"""
Gruneisen with Phonopy and AbiPy
================================

Compute Gruneisedn parameters with phonopy (supercells and finite-difference method).
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.flowtk.abiphonopy import PhonopyGruneisenWork


def build_flow(options):
    """
    Create a `Flow` for phonon calculations with phonopy:
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        __file__ = os.path.join(os.getcwd(), "run_gruneisen_phonopy_si.py")
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize structure and pseudos
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    # Build input for GS calculation.
    gsinp = abilab.AbinitInput(structure, pseudos)
    gsinp.set_vars(ecut=8, nband=4, toldff=1.e-6)

    # This gives ngkpt = 4x4x4 with 4 shifts for the initial unit cell.
    # The k-point sampling will be rescaled when we build the supercell in PhonopyWork.
    gsinp.set_autokmesh(nksmall=4)

    flow = flowtk.Flow(workdir=options.workdir)

    # Use a 2x2x2 supercell to compute phonons with phonopy
    work = PhonopyGruneisenWork.from_gs_input(gsinp, voldelta=0.01, scdims=[2, 2, 2])
    flow.register_work(work)

    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).graphviz_imshow()


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
