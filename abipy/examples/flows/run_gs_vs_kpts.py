#!/usr/bin/env python
r"""
Convergence study with different k-meshes
==========================================

In this example,
"""
import sys
import os

import abipy.abilab as abilab
import abipy.flowtk as flowtk
import abipy.data as abidata


def build_flow(options):
    """
    Build and return a flow for GS properties with different k-meshes
    """
    # Set working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    structure = abilab.Structure.from_file(abidata.cif_file("gan2.cif"))
    pseudos = abidata.pseudos("Ga.oncvpsp", "N.oncvpsp")

    scf_input = abilab.AbinitInput(structure, pseudos)

    scf_input.set_vars(
        ecut=15,       # Too low.
        nstep=50,      # Increase default
        tolvrs=1e-8,
    )

    flow = flowtk.Flow(workdir=options.workdir)

    from abipy.flowtk.gs_works import GsKmeshConvWork
    nksmall_list = [2, 4, 6, 8]
    flow.register_work(GsKmeshConvWork.from_scf_input(scf_input, nksmall_list))

    return flow



# This block generates the thumbnails in the AbiPy gallery.
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

