#!/usr/bin/env python
r"""
Flow for quasi-harmonic calculations
====================================

Warning: This code is still under development.
"""
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
        __file__ = os.path.join(os.getcwd(), "run_qha.py")
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize structure and pseudos
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    # Build input for GS calculation.
    scf_input = abilab.AbinitInput(structure, pseudos)
    scf_input.set_vars(ecut=12, nband=8, tolvrs=1e-8)
    scf_input.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0, 0, 0])

    v0 = scf_input.structure.volume
    volumes = [0.08 * v0, v0, v0 * 1.02]
    return QhaFlow.from_scf_input(options.workdir, scf_input, volumes,
                                  ngqpt=[2, 2, 2], with_becs=False, edos_ngkpt=(4, 4, 4))


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
