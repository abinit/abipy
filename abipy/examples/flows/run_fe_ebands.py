#!/usr/bin/env python
r"""
Band structure w/wo magnetization
=================================

Calculation of the band structure of Fe with and without magnetization,
including L-projected (FATBANDS and FATDOS)
See also <~abinit/tutorial/Input/tspin_1.in>
"""
import os
import sys
import abipy.data as data
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def make_scf_input(nsppol, paral_kgb=1):
    """
    Generate input file for GS and given `nsppol`.
    """
    # Fe normal bcc structure for test of a ferromagnetic calculation
    scf_input = abilab.AbinitInput(structure=data.structure_from_ucell("Fe-fm"),
                                   pseudos=data.pseudos("26fe.pspnc"))

    # Global variables
    global_vars = dict(
        nsppol=nsppol,
        ecut=18,
        nband=8,
        occopt=3,
        tsmear=0.01,
        paral_kgb=paral_kgb,
    )

    if nsppol == 2:
        global_vars.update(spinat=[0.0, 0.0, 4.0])

    scf_input.set_vars(global_vars)

    scf_input.set_kmesh(ngkpt=[4, 4, 4], shiftk=[0.5, 0.5, 0.5])
    scf_input.set_vars(tolvrs=1e-6)

    return scf_input


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    # Create the Flow.
    flow = flowtk.Flow(options.workdir, manager=options.manager)

    for nsppol in [1, 2]:
        # Build a BandStructureWork from the scf_input with the given nsppol and add it to the flow
        # L-projection (prtdos 3) is used by default.
        scf_input = make_scf_input(nsppol)
        work = flowtk.BandStructureWork.from_scf_input(scf_input, dos_ngkpt=(8, 8, 8))
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

############################################################################
#
# Run the script with:
#
#     run_fe_ebands -s
#
# then use:
#
#    abirun.py flow_fe_ebands ebands -p
#
# to analyze all the band structures produced by the Flow and plot the data
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_fe_ebands.png?raw=true
#    :alt: Band structures of Fe with nsppol 1, 2.
#
