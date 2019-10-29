#!/usr/bin/env python
r"""
Effective masses with DFPT
==========================

Flow to compute the band structure of silicon.
"""

import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def make_scf_input(usepaw=0):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    #pseudos = abidata.pseudos("14si.pspnc") if usepaw == 0 else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")
    pseudos = abidata.pseudos("Si_r.psp8") if usepaw == 0 else abidata.pseudos("Si.GGA_PBE-JTH-paw.xml")

    #structure = abidata.cif_file("si.cif"),

    structure = dict(
         ntypat=1,
         natom=2,
         typat=[1, 1],
         znucl=14,
         #acell=3 * [10.26310667319252], # https://docs.abinit.org/tests/v7/Input/t82.in
         acell=3 * [10.2073557], # 5.4015 Ang
         rprim=[[0.0,  0.5,  0.5],
                [0.5,  0.0,  0.5],
                [0.5,  0.5,  0.0]],
         xred=[ [0.0 , 0.0 , 0.0],
                [0.25, 0.25, 0.25]],
    )

    # Get structure from cif file.
    scf_input = abilab.AbinitInput(structure=structure, pseudos=pseudos)

    # Global variables
    ecut = 8
    scf_input.set_vars(
        ecut=ecut,
        nband=8,
        #nband=16,
        #nspinor=2,
        nstep=100,
    )

    if scf_input.ispaw:
        scf_input.set_vars(pawecutdg=2 * ecut)

    # Dataset 1 (GS run)
    scf_input.set_kmesh(ngkpt=[8, 8, 8], shiftk=[0, 0, 0])
    scf_input.set_vars(tolvrs=1e-8)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0],  # L point
        [0.0, 0.0, 0.0],  # Gamma point
        [0.0, 0.5, 0.5],  # X point
    ]

    return scf_input


def build_flow(options):
    # Set working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        __file__ = os.path.join(os.getcwd(), "run_effmass_dfpt.py")
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Get the SCF and the NSCF input.
    scf_input = make_scf_input(usepaw=1)

    # Build the flow.
    from abipy.flowtk.effmass_works import EffMassDFPTWork, EffMassAutoDFPTWork

    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    work = EffMassDFPTWork.from_scf_input(scf_input, k0_list=(0, 0, 0), effmass_bands_f90=[1, 4],
                                          #red_dirs=[[1, 0, 0], [1, 1, 0]],
                                          #red_dirs=None,
                                          #cart_dirs=[[1, 0, 0], [1, 1, 1], [1, 1, 0]],
                                          #den_node=None

                                          )

    work = EffMassAutoDFPTWork.from_scf_input(scf_input, ndivsm=5, tolwfr=1e-12)

    flow.register_work(work)

    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    #build_flow(options).plot_networkx(with_edge_labels=True, tight_layout=True)
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
#     run_si_ebands.py -s
#
# then use:
#
#    abirun.py flow_si_ebands ebands --plot
#
# to analyze (and plot) the electronic bands produced by the Flow.
#
# .. code-block:: bash
#
#    KS electronic bands:
#           nsppol  nspinor  nspden  nkpt  nband  nelect  fermie formula  natom  \
#    w0_t0       1        1       1    29      8     8.0   5.598     Si2      2
#    w0_t1       1        1       1    14      8     8.0   5.598     Si2      2
#
#           angle0  angle1  angle2      a      b      c  volume abispg_num scheme  \
#    w0_t0    60.0    60.0    60.0  3.867  3.867  3.867  40.888        227   none
#    w0_t1    60.0    60.0    60.0  3.867  3.867  3.867  40.888        227   none
#
#           occopt  tsmear_ev  bandwidth_spin0  fundgap_spin0  dirgap_spin0  \
#    w0_t0       1      0.272           11.856          0.562         2.532
#    w0_t1       1      0.272           11.856          0.524         2.532
#
#          task_class                                   ncfile              status
#    w0_t0    ScfTask  flow_si_ebands/w0/t0/outdata/out_GSR.nc  Completed
#    w0_t1   NscfTask  flow_si_ebands/w0/t1/outdata/out_GSR.nc  Completed
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_si_ebands.png?raw=true
#    :alt: Band structure of Si in the IBZ and along a k-path
#
