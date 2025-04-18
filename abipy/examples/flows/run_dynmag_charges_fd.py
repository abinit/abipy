#!/usr/bin/env python
r"""
Dynamical magnetic charges with finite difference
==========================

Flow to compute dynamical magnetic charges with finite difference
"""

import sys
import os
#import abipy.data as abidata
#import abipy.abilab as abilab
import abipy.flowtk as flowtk

from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput
from abipy.flowtk.finitediff_works import FdDynMagneticChargeWork


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
                                  relativity_type="SR", accuracy="standard")
                                  #relativity_type="FR", accuracy="standard")


    # This is a calculation with spin-up and spin-down wavefunctions,         ... nsppol=  2
    # in which the occupation numbers are to be determined automatically.     ... occopt=  1
    # However, in this case, the target total spin magnetization
    # must be specified, while the default value is observed.                 ... spinmagntarget= -99.99
    # Action: if you are doing an antiferromagnetic calculation, please use nsppol=1 with nspden=2;
    # on the other hand, if you are doing a ferromagnetic calculation, either specify your own spinmagntarget,
    # or let the code determine the total spin-polarization, by using a metallic value for occopt (e.g. 7 or 4 ...).

    nspinor = 1
    nsppol, nspden = 1, 4
    if nspinor == 1:
        nsppol, nspden  = 2, 2

    scf_input = AbinitInput(structure, pseudos)
    #scf_input.set_spin_mode(self, spin_mode)
    # AFM
    nspinor, nsppol, nspden  = 1, 1, 2

    scf_input.set_vars(
        ecut=43,
        #ecut=12,
        nband=60,          #Cr.psp8:14; O.psp8:6; (14*6+6*4)/2=54;54(occupied)+6(unoccupied)=60
        tolvrs=1e-8,       # SCF stopping criterion
        #tolvrs=1e-2,       # SCF stopping criterion
        nspinor=nspinor,
        nsppol=nsppol,
        nspden=nspden,
        nstep=800,         # Maximal number of SCF cycles
        diemac=12.0,       # Although this is not mandatory, it is worth to
        #so_psp=??
        paral_kgb=0,
    )

    scf_input.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0, 0, 0])

    scf_input["spinat"] = """
        0.0 0.0  -4.0
        0.0 0.0   4.0
        0.0 0.0   4.0
        0.0 0.0  -4.0
        18*0
    """

    # Initialize the flow
    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    work = FdDynMagneticChargeWork.from_scf_input(
        scf_input,
        berryopt=-1,
        num_points=3,
        delta_h=0.01,
        relax_opts=None
    )

    # Add the work to the flow
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
