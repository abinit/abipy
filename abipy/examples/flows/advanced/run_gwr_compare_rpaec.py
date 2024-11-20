#!/usr/bin/env python
r"""
RPA flow with convergence studies
=================================

This script computes the RPA energy with the GWR code
and the quartic scaling algorithm (Adler-Wiser exact expression in frequency domain)
so that one can then compare the results.
"""

import os
import sys
import abipy.data as data
import abipy.abilab as abilab

from abipy import flowtk


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_","flow_")

    from abipy.data.gwr_structures import get_gwr_structure
    symbol = "Si"
    #symbol = "LiF"
    #symbol = "Si"
    #symbol = "C"
    #symbol = "BN"
    #symbol = "MgO"
    #symbol = "GaAs"
    #symbol = "ZnO"
    structure = get_gwr_structure(symbol)

    from abipy.flowtk.psrepos import get_repo_from_name
    #pseudos = get_repo_from_name("ONCVPSP-PBE-SR-PDv0.4").get_pseudos("stringent")
    pseudos = get_repo_from_name("ONCVPSP-PBE-SR-PDv0.4").get_pseudos("standard")

    scf_input = abilab.AbinitInput(structure=structure, pseudos=pseudos)
    scf_input.set_cutoffs_for_accuracy("normal")
    ecut = scf_input["ecut"]
    scf_input.set_scf_nband_semicond()

    # Global variables.
    scf_input.set_vars(
        tolvrs=1e-8,
        paral_kgb=0,
        npfft=1,
        timopt=-1,
    )
    scf_input.set_kmesh(
        ngkpt=[1, 1, 1],
        #ngkpt=[4, 4, 4],
        #ngkpt=[8, 8, 8],
        shiftk=[0.0, 0.0, 0.0], # IMPORTANT: k-grid for GWR must be Gamma-centered.
    )

    # Get max number of PWs.
    dims, _ = scf_input.abiget_dims_spginfo()
    mpw = dims["mpw"]

    flow = flowtk.Flow(workdir=options.workdir)
    small_manager = options.manager.new_with_fixed_mpi_omp(4, 1)

    # GS-SCF run to get the DEN, followed by direct diago to obtain green_nband bands.
    from abipy.flowtk.gwr_works import DirectDiagoWork, GWRRPAConvWork
    green_nband = -1  # -1 this means full diago
    diago_work = DirectDiagoWork.from_scf_input(scf_input, green_nband)
    wfk_node = diago_work[1]
    diago_work[0].set_manager(small_manager)
    flow.register_work(diago_work)

    gwr_ntau_list = list(range(6, 34, 2))
    gwr_ntau_list = [6, 8]
    nband = int(mpw*0.8)
    ecuteps = 12

    # GWR with different values of regterm
    for gwr_regterm in (0.0, 1e-6):
        work = GWRRPAConvWork.from_scf_input_ntaus(scf_input, gwr_ntau_list,
                                                   nband=nband, ecuteps=ecuteps,
                                                   den_node=diago_work[0], wfk_node=wfk_node,
                                                   gwr_kwargs=dict(gwr_regterm=gwr_regterm))
        flow.register_work(work)

    # RPA with quartic scaling. See also https://docs.abinit.org/tests/v67mbpt/Input/t19.abi
    work = flow.new_work()
    for gwr_ntau in gwr_ntau_list:
        chi_input = scf_input.new_with_vars(optdriver=3,       # Screening calculation.
                                            gwcalctyp=1,       # Gauss-Legendre frequency mesh on the imaginary axis.
                                            gwrpacorr=1,       # Exact integration over the coupling constant.
                                            nfreqim=gwr_ntau,  # No. of points along the imaginary axis for chi0.
                                            nband=nband,
                                            ecuteps=ecuteps,
                                           )
        work.register_scr_task(chi_input, deps={wfk_node: "WFK"})

    return flow


@flowtk.flow_main
def main(options):
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
