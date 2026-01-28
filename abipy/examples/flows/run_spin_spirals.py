#!/usr/bin/env python
r"""
Spin spirals in iron with GBT
=============================

"""
import os
import sys
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    # Fe in bcc structure.
    structure = abilab.Structure.from_abistring("""
natom 1
ntypat 1
typat 1
znucl 26
xred    0.0000000000    0.0000000000    0.0000000000
acell    1.0    1.0    1.0
rprim
   4.2302629971    0.0000000000    2.4423434801
   1.4100876657    3.9883302019    2.4423434801
   0.0000000000    0.0000000000    4.8846869602
""")

    # gamma-iron. fcc structure with Cu lattice parameters
    a = 6.82
    #a = 6.833
    structure = abilab.Structure.fcc(a, ["Fe"], units="bohr")
    print(structure)
    print(structure.abi_string)

    # Use relativistic NC PBEsol pseudos from pseudodojo v0.4.
    from abipy.flowtk.psrepos import get_oncvpsp_pseudos
    pseudos = get_oncvpsp_pseudos(xc_name="PBEsol", version="0.4", relativity_type="FR")

    scf_input = abilab.AbinitInput(structure=structure, pseudos=pseudos)

    connect = False

    scf_input.set_vars(
        # SPIN
        nspinor=2,
        nspden=4,
        so_psp=0,
        spinat=[4.0, 0.0, 0.0],
        ixc=7,        # Use LDA instead of PBEsol
        ecut=30,
        nband=24,
        nline=12,         # To facilitate convergence.
        nstep=100,        # Maximal number of SCF cycles
        toldfe=1e-6,
        #
        nsym=1,           # Disable spatial symmetries
        #
        occopt=7,
        tsmear=0.01,
        paral_kgb=0,
        prtwf=-1 if connect else 1,
    )

    # Create the Flow.
    flow = flowtk.Flow(options.workdir, manager=options.manager)

    # Define q-path. Two modes are available
    # If qnames is None, a standardized path is automatically generated from the structure
    # else the labels are taken from qnames.
    qnames = None
    qnames = ["G", "X"]
    line_density = 10

    # This to show how to perform convergence studies.
    # Here we compute E(q) for different ngkpt
    ngkpt_list = [
        (2, 2, 2),
        #(4, 4, 4),
    ]

    from abipy.flowtk.gs_works import SpinSpiralWork
    for ngkpt in ngkpt_list:
        new_input = scf_input.deepcopy()
        new_input.set_kmesh(ngkpt=ngkpt, shiftk=[0.0, 0.0, 0.0], kptopt=4)
        work = SpinSpiralWork.from_scf_input(new_input, qnames=qnames, line_density=line_density, connect=connect)
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
