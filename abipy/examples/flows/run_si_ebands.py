#!/usr/bin/env python
r"""
Band structure Flow
===================

Flow to compute the band structure of silicon.
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def make_scf_nscf_inputs(paral_kgb=0, usepaw=0):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    pseudos = abidata.pseudos("14si.pspnc") if usepaw == 0 else data.pseudos("Si.GGA_PBE-JTH-paw.xml")

    # Get structure from cif file.
    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"), pseudos=pseudos, ndtset=2)

    # Global variables
    ecut = 6
    multi.set_vars(
        ecut=ecut,
        nband=8,
        paral_kgb=paral_kgb,
        iomode=3,
        timopt=-1,
    )

    if multi.ispaw:
        multi.set_vars(pawecutdg=2 * ecut)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(ngkpt=[8, 8, 8], shiftk=[0, 0, 0])
    multi[0].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0],  # L point
        [0.0, 0.0, 0.0],  # Gamma point
        [0.0, 0.5, 0.5],  # X point
    ]

    multi[1].set_kpath(ndivsm=6, kptbounds=kptbounds)
    multi[1].set_vars(tolwfr=1e-12)

    # Generate two input files for the GS and the NSCF run
    scf_input, nscf_input = multi.split_datasets()
    return scf_input, nscf_input


def build_flow(options):
    # Set working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Get the SCF and the NSCF input.
    scf_input, nscf_input = make_scf_nscf_inputs()

    # Build the flow.
    return flowtk.bandstructure_flow(options.workdir, scf_input, nscf_input, manager=options.manager)


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
