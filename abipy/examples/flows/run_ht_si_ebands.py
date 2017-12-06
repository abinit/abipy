#!/usr/bin/env python
r"""
Band structure Flow with factory functions
==========================================

Band structure of silicon with the HT interface.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy import abilab

def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_")

    # Initialize structure and pseudos.
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    # Initialize the flow.
    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    # Use ebands_input factory function to build inputs.
    multi = abilab.ebands_input(structure, pseudos, kppa=40, nscf_nband=6, ndivsm=10, ecut=6)
    work = flowtk.BandStructureWork(scf_input=multi[0], nscf_input=multi[1])

    flow.register_work(work)
    return flow


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
