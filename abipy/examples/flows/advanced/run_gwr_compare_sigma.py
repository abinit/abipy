#!/usr/bin/env python
r"""
GWR flow with convergence studies
=================================

This script computes the G0W0 corrections with the GWR code
and the analytic continuation implemented in quartic scaling code in order to compare the results
"""

import os
import sys
import abipy.data as data
import abipy.abilab as abilab
import abipy.core.abinit_units as abu

from abipy import flowtk


def build_flow(options):

    from abipy.data.gwr_structures import get_gwr_structure
    symbol = "Si"
    #symbol = "LiF"
    #symbol = "Si"
    #symbol = "C"
    #symbol = "BN"
    #symbol = "MgO"
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
        paral_kgb=0,
        npfft=1,
        timopt=-1,
    )
    scf_input.set_kmesh(
        #ngkpt=[1, 1, 1],
        ngkpt=[2, 2, 2],
        #ngkpt=[4, 4, 4],
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
    #diago_work[0].with_fixed_mpi_omp(1, 1)
    flow.register_work(diago_work)

    # Build template for GWR.
    ecuteps = 12

    gwr_template = scf_input.make_gwr_qprange_input(gwr_ntau=6, nband=int(mpw * 0.9), ecuteps=ecuteps)

    # Define kptgw and bdgw
    nval = scf_input.num_valence_electrons // 2
    kptgw = [ # k-points in reduced coordinates
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0), # X
        #(0.5    0.000    0.000),
    ]

    nkptgw = len(kptgw)
    bdgw = (nval, nval+1) * nkptgw

    sigma_kcalc_dict = dict(
        nkptgw=nkptgw,
        kptgw=kptgw,
        bdgw=bdgw,
    )

    gwr_template.set_vars(**sigma_kcalc_dict)

    gwr_ntau_list = list(range(6, 34, 2))
    gwr_ntau_list = [6, 8]

    # Conpute QP corrections without/with regularization term.
    # 1) Change the value of one variable:
    varname_values = ("gwr_ntau", gwr_ntau_list)

    nband = nval * 100

    gwr_template["nband"] = nband
    gwr_template["userra"] = 1e-6
    #gwr_template["userra"] = 0.0

    wfk_node = diago_work[1]
    gwr_work = GWRSigmaConvWork.from_varname_values(
            varname_values, gwr_template, den_node=diago_work[0], wfk_node=wfk_node)
    flow.register_work(gwr_work)

    # Create work for AC calculation
    #nfreqim  25      # This is equal to gwr_ntau
    #nomegasi  10
    #omegasimax 10 eV
    work = flow.new_work()
    for gwr_ntau in gwr_ntau_list:
        chi_ac_input = scf_input.new_with_vars(optdriver=3,
                                               gwcalctyp=1, # Analytic continuation.
                                               nfreqim=gwr_ntau,
                                               ecuteps=ecuteps,
                                               nband=nband,
                                               )

        sigma_ac_input = chi_ac_input.new_with_vars(optdriver=4,
                                                    nomegasi=gwr_ntau,
                                                    omegasimax=0.2 * abu.eV_Ha * gwr_ntau,
                                                    ecutsigx=gwr_template["ecutsigx"],
                                                   )
        sigma_ac_input.set_vars(**sigma_kcalc_dict)
        scr_ac_task = work.register_scr_task(chi_ac_input, deps={wfk_node: "WFK"})
        work.register_sigma_task(sigma_ac_input, deps={wfk_node: "WFK", scr_ac_task: "SCR"})

    # Create work for plasmon-pole calculation
    chi_ppm_input = scf_input.new_with_vars(optdriver=3,
                                            ecuteps=ecuteps,
                                            nband=nband,
                                            )

    sigma_ppm_input = chi_ppm_input.new_with_vars(optdriver=4,
                                                  ecutsigx=gwr_template["ecutsigx"],
                                                 )
    sigma_ppm_input.set_vars(**sigma_kcalc_dict)

    work = flow.new_work()
    scr_ppm_task = work.register_scr_task(chi_ppm_input, deps={wfk_node: "WFK"})
    work.register_sigma_task(sigma_ppm_input, deps={wfk_node: "WFK", scr_ac_task: "SCR"})

    return flow


@flowtk.flow_main
def main(options):
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
