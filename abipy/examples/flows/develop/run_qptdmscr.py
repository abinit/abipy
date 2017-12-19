#!/usr/bin/env python
r"""
Screening Flow with q-point parallelism
=======================================

This example shows how to compute the SCR file by splitting the calculation of the SCR file
over q-points with the input variables nqptdm and qptdm.
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as data

from abipy import flowtk


def all_inputs(paral_kgb=1):
    """
    Build the input files of the calculation.
    Returns: gs_input, nscf_input, scr_input, sigma_input
    """
    ecut = ecutwfn = 4
    global_vars = dict(
        ecut=ecut,
        istwfk="*1",
        paral_kgb=paral_kgb,
    )

    multi = abilab.MultiDataset(structure=data.cif_file("si.cif"),
                                pseudos=data.pseudos("14si.pspnc"), ndtset=4)
    multi.set_vars(global_vars)

    gs, nscf, scr, sigma = multi.split_datasets()

    # This grid is the most economical, but does not contain the Gamma point.
    gs_kmesh = dict(
        ngkpt=[2, 2, 2],
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap.
    gw_kmesh = dict(
        ngkpt=[2, 2, 2],
        shiftk=[0.0, 0.0, 0.0,
                0.0, 0.5, 0.5,
                0.5, 0.0, 0.5,
                0.5, 0.5, 0.0]
    )

    # Dataset 1 (GS run)
    gs.set_kmesh(**gs_kmesh)
    gs.set_vars(tolvrs=1e-6, nband=4)

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    nscf.set_kmesh(**gw_kmesh)
    nscf.set_vars(iscf=-2,
                  tolwfr=1e-10,
                  nband=15,
                  nbdbuf=5,
    )

    # Dataset3: Calculation of the screening.
    scr.set_kmesh(**gw_kmesh)

    scr.set_vars(
        optdriver=3,
        nband=6,
        ecutwfn=ecutwfn,
        symchi=1,
        inclvkb=0,
        ecuteps=2.0,
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    sigma.set_kmesh(**gw_kmesh)

    sigma.set_vars(
            optdriver=4,
            nband=8,
            ecutwfn=ecutwfn,
            ecuteps=2.0,
            ecutsigx=2.0,
            symsigma=1,
    )

    kptgw = [
         -2.50000000E-01, -2.50000000E-01,  0.00000000E+00,
         -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
          5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
         -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
          5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
          0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
    ]

    bdgw = [1, 8]

    sigma.set_kptgw(kptgw, bdgw)

    return multi.split_datasets()


def build_flow(options):
    """
    Construct the flow for G0W0 calculations.
    The calculation of the SCR file is split into nqptdm independent calculations.
    Partial SCR files are then merged with mrgddb before starting the sigma run.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Build the input files for GS, NSCF, SCR and SIGMA runs.
    gs, nscf, scr_input, sigma_input = all_inputs()

    # Construct the flow.
    return flowtk.G0W0WithQptdmFlow(options.workdir, gs, nscf, scr_input, sigma_input,
                                    manager=options.manager)


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("GENERATE_SPHINX_GALLERY", False):
    __name__ = None
    #import tempfile
    #options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    #build_flow(options).plot_networkx(with_edge_labels=True, tight_layout=True)


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
