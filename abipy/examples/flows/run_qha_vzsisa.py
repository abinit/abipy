#!/usr/bin/env python
r"""
Flow for QHA calculations v-ZSISA-QHA
=====================================

See [Phys. Rev. B 110, 014103](https://doi.org/10.1103/PhysRevB.110.014103)
"""
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk
from abipy.flowtk.vzsisa import VzsisaFlow


def build_flow(options):
    """
    Create a `VzsisaFlow`.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        __file__ = os.path.join(os.getcwd(), "run_qha_vzsisa.py")
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize structure from cif file.
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    # Get NC pseudos from pseudodojo.
    from abipy.flowtk.psrepos import get_oncvpsp_pseudos
    pseudos = get_oncvpsp_pseudos(xc_name="PBE", version="0.4")

    # Select k-mesh for electrons and q-mesh for phonons.
    ngkpt = [2, 2, 2]; ngqpt = [1, 1, 1]
    #ngkpt = [4, 4, 4]; ngqpt = [4, 4, 4]

    with_becs = False
    with_quad = False
    #with_quad = not structure.has_zero_dynamical_quadrupoles

    #bo_vol_scales = [0.96, 0.98, 1.0, 1.02, 1.04, 1.06]
    #ph_vol_scales = [0.98, 1.0, 1.02, 1.04, 1.06] # EinfVib4(D)
    bo_vol_scales = [0.96, 0.98, 1, 1.02, 1.04]    # EinfVib4(S)
    ph_vol_scales = [1, 1.02, 1.04]                # EinfVib2(D)

    scf_input = abilab.AbinitInput(structure, pseudos)

    # Set other important variables
    # All the other DFPT runs will inherit these parameters.
    scf_input.set_vars(
        nband=4,
        nline=10,
        nbdbuf=0,
        nstep=100,
        ecut=8.0,
        ecutsm=1.0,
        occopt=1,
        tolvrs=1.0e-18,   # SCF stopping criterion (modify default)
    )

    scf_input.set_kmesh(ngkpt=ngkpt, shiftk=[0, 0, 0])
    print("scf_input\n", scf_input)

    return VzsisaFlow.from_scf_input(options.workdir, scf_input, bo_vol_scales, ph_vol_scales, ngqpt,
                                     with_becs, with_quad, edos_ngkpt=None)


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
