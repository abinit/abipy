#!/usr/bin/env python
r"""
Flow for quasi-harmonic calculations
====================================

Warning: This code is still under development.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata
from abipy import flowtk

from abipy.flowtk.qha import QhaFlow


def build_flow(options):
    """
    Create a `QhaFlow` for quasi-harmonic calculations.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize structure and pseudos
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    # Build input for GS calculation.
    gsinp = abilab.AbinitInput(structure, pseudos)
    gsinp.set_vars(ecut=4, nband=4, toldff=1.e-6)
    gsinp.set_autokmesh(nksmall=2)

    volumes = [gsinp.structure.volume]
    return QhaFlow.from_gsinp(options.workdir, gsinp, volumes, ngqpt=[2,2,2])


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("GENERATE_SPHINX_GALLERY", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).plot_networkx(tight_layout=True)


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
