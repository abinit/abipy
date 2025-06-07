#!/usr/bin/env python
r"""
Flow for quasi-harmonic calculations under development
======================================================

Warning: This code is still under development.
"""
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk
from abipy.flowtk.zsisa import ZsisaFlow


def build_flow(options):
    """
    Create a `QhaFlow` for quasi-harmonic calculations.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        __file__ = os.path.join(os.getcwd(), "run_qha_zsisa.py")
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize structure and pseudos for ZnO.
    structure = abilab.Structure.from_abistring("""
natom 4
ntypat 2
typat
1 1 2
2
znucl 30 8
xred
   0.0000000000    0.0000000000   -0.0025137620
   0.3333333333    0.6666666667    0.4974862380
   0.0000000000    0.0000000000    0.3835201241
   0.3333333333    0.6666666667    0.8835201241
acell    1.0    1.0    1.0
rprim
   6.3016720624    0.0000000000    0.0000000000
  -3.1508360312    5.4574080923    0.0000000000
   0.0000000000    0.0000000000    9.7234377918
""")

    # Initialize structure and pseudos
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    # Use NC PBE pseudos from pseudodojo v0.4
    from abipy.flowtk.psrepos import get_oncvpsp_pseudos
    pseudos = get_oncvpsp_pseudos(xc_name="PBEsol", version="0.4")

    # Select k-mesh for electrons and q-mesh for phonons.
    #ngkpt = [6, 6, 4]; ngqpt = [1, 1, 1]
    ngkpt = [2, 2, 2]; ngqpt = [1, 1, 1]

    with_becs = False
    with_quad = False
    #with_quad = not structure.has_zero_dynamical_quadrupoles

    scf_input = abilab.AbinitInput(structure, pseudos)

    # Set other important variables
    scf_input.set_vars(
        nband=scf_input.num_valence_electrons // 2,
        nline=10,
        nbdbuf=0,
        nstep=100,
        ecutsm=1.0,
        #tolvrs=1.0e-18,      # SCF stopping criterion (modify default)
        tolvrs=1.0e-6,      # SCF stopping criterion (modify default)
    )

    scf_input.set_kmesh(ngkpt=ngkpt, shiftk=[0, 0, 0])

    # mode: 'TEC' for thermal expansion only, 'ECs' to include elastic constants
    mode = 'TEC'
    eps = 0.005

    flow = ZsisaFlow.from_scf_input(options.workdir, scf_input, eps, mode, ngqpt,
                                    with_becs, with_quad, edos_ngkpt=None)

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
