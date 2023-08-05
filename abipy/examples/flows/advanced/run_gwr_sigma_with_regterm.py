#!/usr/bin/env python
r"""
GWR flow with convergence studies
=================================

This script shows how to compute the G0W0 corrections in silicon.
More specifically, we build a flow to analyze the convergence of the QP corrections
wrt to the number of bands in the self-energy. More complicated convergence studies
can be implemented on the basis of this example.
"""

import os
import sys
import abipy.data as data
import abipy.abilab as abilab

from abipy import flowtk


def build_flow(options):

    from abipy.data.gwr_structures import get_gwr_structure
    symbol = "Si"
    structure = get_gwr_structure(symbol)

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_","flow_")

    from abipy.flowtk.psrepos import get_repo_from_name
    pseudos = get_repo_from_name("ONCVPSP-PBE-SR-PDv0.4").get_pseudos("stringent")

    scf_input = abilab.AbinitInput(structure=structure, pseudos=pseudos)
    scf_input.set_cutoffs_for_accuracy("normal")
    ecut = scf_input["ecut"]
    scf_input.set_scf_nband_semicond()

    # Global variables.
    scf_input.set_vars(
        tolvrs=1e-8,
        paral_kgb=1,
    )
    scf_input.set_kmesh(
        ngkpt=[2, 2, 2],
        #ngkpt=[8, 8, 8],
        shiftk=[0.0, 0.0, 0.0], # IMPORTANT: k-grid for GWR must be Gamma-centered.
    )

    # Get max number of PWs.
    dims, _ = scf_input.abiget_dims_spginfo()
    mpw = dims["mpw"]

    flow = flowtk.Flow(workdir=options.workdir)

    # GS-SCF run to get the DEN, followed by direct diago to obtain green_nband bands.
    from abipy.flowtk.gwr_works import DirectDiagoWork, GWRSigmaConvWork
    green_nband = -1  # -1 this means full diago
    diago_work = DirectDiagoWork.from_scf_input(scf_input, green_nband)
    flow.register_work(diago_work)

    # Build template for GWR.
    ecuteps = 12
    gwr_template = scf_input.make_gwr_qprange_input(gwr_ntau=6, nband=int(mpw * 0.9), ecuteps=ecuteps)

    gwr_ntau_list = list(range(6, 34, 2))
    gwr_ntau_list = [6, 8]

    # 1) Change the value of one variable:
    varname_values = ("gwr_ntau", gwr_ntau_list)
    # or take the Cartesian product of two or more variables with e.g.:
    #
    #varname_values = [
    #   ("gwr_ntau", gwr_ntau_list),
    #   #("userra", [0.0, 1e-6),    # Compute QP corrections with/without regterm.
    #   #("ecuteps", [4, 6]),
    #]

    # Conpute QP corrections without/with regularization term.
    for userra in [0.0, 1e-6]:
        gwr_template["userra"] = userra
        gwr_work = GWRSigmaConvWork.from_varname_values(
                varname_values, gwr_template, den_node=diago_work[0], wfk_node=diago_work[1])
        flow.register_work(gwr_work)

    return flow


@flowtk.flow_main
def main(options):
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
