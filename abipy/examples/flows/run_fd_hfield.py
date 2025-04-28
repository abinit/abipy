#!/usr/bin/env python
r"""
Dynamical magnetic charges
==========================

Crystalline Cr2O3 with magnetic field (zeemanfield)
Flow to compute dynamical magnetic charges with finite differences.

Z_jv^m=Ω_0 (∂M_v)/(∂u_j ) = (∂F_j)/(∂H_v ) = Ω_0 (∂^2 E)/(∂H_β ∂u_i ).
"""

import sys
import os
import abipy.flowtk as flowtk

from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput
from abipy.flowtk.finitediff_works import FiniteHfieldWork


def build_flow(options):
    # Set working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    structure = Structure.from_abistring("""
natom 10
ntypat 2
typat
1 1 1
1 2 2
2 2 2
2
znucl 24 8
xred      3.4574472807E-01  3.4574472807E-01  3.4574472807E-01
          1.5425527193E-01  1.5425527193E-01  1.5425527193E-01
          6.5425527193E-01  6.5425527193E-01  6.5425527193E-01
          8.4574472807E-01  8.4574472807E-01  8.4574472807E-01
          5.5995881675E-01  2.5000000000E-01  9.4004118325E-01
          2.5000000000E-01  9.4004118325E-01  5.5995881675E-01
          9.4004118325E-01  5.5995881675E-01  2.5000000000E-01
          4.4004118325E-01  7.5000000000E-01  5.9958816751E-02
          5.9958816751E-02  4.4004118325E-01  7.5000000000E-01
          7.5000000000E-01  5.9958816751E-02  4.4004118325E-01

acell    3*1.0223825450E+01
rprim     5.2802747870E-01  0.0000000000E+00  8.4922728509E-01
         -2.6401373935E-01  4.5728521045E-01  8.4922728509E-01
         -2.6401373935E-01 -4.5728521045E-01  8.4922728509E-01
""")

    # Get NC pseudos from pseudodojo.
    from abipy.flowtk.psrepos import get_oncvpsp_pseudos
    pseudos = get_oncvpsp_pseudos(xc_name="PBEsol", version="0.4",
                                  relativity_type="FR", accuracy="standard")
                                  #relativity_type="SR", accuracy="standard")
    #nspinor = 1
    #nsppol, nspden = 1, 4
    #if nspinor == 1:
    #    nsppol, nspden  = 2, 2

    scf_input = AbinitInput(structure, pseudos)
    num_ele = scf_input.num_valence_electrons
    #scf_input.set_spin_mode(self, spin_mode)
    # AFM
    #nspinor, nsppol, nspden  = 1, 1, 2
    nspinor, nsppol, nspden  = 2, 1, 4
    metallic = True

    nband = int(1.1 * num_ele) if metallic else num_ele
    nband += nband % 2

    scf_input.set_vars(
        #ecut=43,
        ecut=12,
        #nband=60,          # Cr.psp8:14; O.psp8:6; (14*6+6*4)/2=54;54(occupied)+6(unoccupied)=60
        nband=nband,
        tolvrs=1e-3,       # SCF stopping criterion
        #tolvrs=1e-8,       # SCF stopping criterion
        nspinor=nspinor,
        nsppol=nsppol,
        nspden=nspden,
        nstep=800,         # Maximal number of SCF cycles
        #nstep=0,         # Maximal number of SCF cycles
        diemac=12.0,
        occopt=7,
        tsmear=0.01,
        ixc=7,
        paral_kgb=0,       # The system is too small to use paral_kgb = 1
    )

    # kptopt=4,          # NO TR.
    scf_input.set_kmesh(ngkpt=[1, 1, 1], shiftk=[0, 0, 0], kptopt=4)
    #scf_input.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0, 0, 0], kptopt=4)

    scf_input["spinat"] = """
        0.0 0.0  -4.0
        0.0 0.0   4.0
        0.0 0.0   4.0
        0.0 0.0  -4.0
        18*0
    """

    # Initialize the flow
    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    work = FiniteHfieldWork.from_scf_input(
        scf_input,
        fd_accuracy=4,
        step_au=0.01,
        relax_ions=False,
        relax_ions_opts=None,
    )

    # Add the work to the flow.
    flow.register_work(work)

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
