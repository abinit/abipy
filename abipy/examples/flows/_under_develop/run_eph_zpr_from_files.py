#!/usr/bin/env python
r"""
Convergence studies for e-ph ZPR from previous DFPT flow
========================================================

This example shows how to build a new Flow to converge the QP corrections due to e-ph
using the DDB/DVDB/POT files produced by a previous flow.

More specifically, we perform NSCF calculations with different dense k-meshes and empty states.
Then we use these dense WFK files to compute the e-ph self-energy by varying the value of zcut and nband.
"""
import sys
import os
import itertools
import numpy as np
import abipy.abilab as abilab

from abipy import flowtk


def build_flow(options):
    """
    Create a `Flow` for ZPR convergene studies
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Recostruct previous DFPT flow from the pickle file.
    # Assume flow in which:
    #   - w0/t0 is the SCF task producing DEN, POT files.
    #   - w1 is the PhononWork producing the DDB and the DVDB file in w1/outdir

    flow_workdir = options.extra
    print(f"Reconstructing Flow from {flow_workdir} ...")
    if not flow_workdir:
        raise ValueError("Please specify the location of the DFPT flow with `--extra FLOWDIR`")

    prev_flow = flowtk.Flow.pickle_load(flow_workdir)

    scf_task = prev_flow[0][0]
    gs_input = scf_task.input
    gs_input.pop_vars(["autoparal", "max_ncpus"])
    dfpt_work = prev_flow[1]

    # Paths to output files.
    ddb_filepath = dfpt_work.outdir.has_abiext("DDB")
    dvdb_filepath = dfpt_work.outdir.has_abiext("DVDB")
    den_filepath = scf_task.outdir.has_abiext("DEN")
    pot_filepath = scf_task.outdir.has_abiext("POT")

    print(f"""
Creating ZPR flow from the following files:

    ddb_filepath  = `{ddb_filepath}`
    dvdb_filepath = `{dvdb_filepath}`
    den_filepath  = `{den_filepath}`
    pot_filepath  = `{pot_filepath}`
""")

    if any(not f for f in (ddb_filepath, dvdb_filepath, den_filepath)):
        raise RuntimeError("DDB, DVDB and DEN file are required to build the ZPR flow.")

    if not pot_filepath:
        raise RuntimeError("POT file is higly recommended for the Sternheimer method")

    # Get ab-initio coarse q-mesh from the DDB file.
    with abilab.abiopen(ddb_filepath) as ddb:
        ddb_ngqpt = ddb.guessed_ngqpt
        #print(ddb)

    # Now we build NSCF tasks with different k-meshes and "enough" empty states.
    # CHANGE this list according to your needs
    ngkpt_fine_list = [
        [8, 8, 8],
        #[12, 12, 12],
        #[32, 32, 32],
    ]

    # By default, all dense k-meshes are Gamma-centered.
    num_kmesh = len(ngkpt_fine_list)
    shiftk_fine_list = np.reshape(num_kmesh * [0, 0, 0], (-1, 3))
    nshiftk_list = [shift.shape[0] for shift in shiftk_fine_list]

    flow = flowtk.Flow(workdir=options.workdir)

    # Compute nband and nbdbuf for the NSCF run.
    # Increase the number of GS-SCF bands by 50 + buffer
    # This should be enough for convergence studies as we're gonna use Sternheimer in the EPH code.
    gs_nband = gs_input["nband"]
    nscf_nband = gs_nband + 50
    nbdbuf = max(int(nscf_nband * 0.1), 4)
    nscf_nband += nbdbuf

    # Pack all NSCF calculations in nscf_work.
    nscf_work = flow.new_work()
    for i, ngkpt_fine in enumerate(ngkpt_fine_list):

        nscf_inp = gs_input.new_with_vars(
            ngkpt=ngkpt_fine,
            shiftk=shiftk_fine_list[i],
            nband=nscf_nband,
            nbdbuf=nbdbuf,      # Reduces considerably the time needed to converge empty states!
            tolwfr=1e-20,
            iscf=-2,
        )
        nscf_work.register_nscf_task(nscf_inp, deps={scf_task: "DEN"})

    # Build template for e-ph self-energy calculations (real + imag part).

    eph_template = gs_input.new_with_vars(
        optdriver=7,             # EPH driver.
        eph_task=4,              # EPH self-energy (Re + Im)
        ddb_ngqpt=ddb_ngqpt,     # Ab-initio q-mesh used for the DDB file.
        tmesh=[0, 100, 4],       # T-mesh in Kelvin (start, step, num)
        eph_stern=1,             # Use Sternheimer equation

        ######################
        # k-points in Sigma_nk

        gw_qprange=0,           # Compute the QP corrections only for the fundamental and the direct gap
                                # gw_qprange is handy as we don't need to specify kptgw and bdgw
                                # but keep in mind that the band edges are computed from the input WFK k-mesh
                                # so their position in the BZ MIGHT CHANGE

        #nkptgw=2,
        #kptgw=[0, 0, 0,        # NB: The k-points must be in the WFK file.
        #       0.5, 0.5, 0],
        #bdgw=[1, 8, 1, 8],

        #####################
        # Spectral function
        #nfreqsp=8000,
        #freqspmax="8.0 eV",

        #####################
        # Advanced options.
        mixprec=1,               # These two varables accelerats the FFT.
        boxcutmin=1.1,
        #symsigma=0,             # Deactivate symmetries in self-energy integration (BZ instead of IBZ_k)
    )

    # Set q-path for Fourier interpolation of phonons.
    eph_template.set_qpath(10)

    # Set q-mesh for phonons DOS.
    eph_template.set_phdos_qmesh(nqsmall=16, method="tetra")

    # Now we use the EPH template to perform a convergence study in which
    # we change the dense q-mesh used to integrate the self-energy and the number of bands.
    # The code will activate the Fourier interpolation of the DFPT potentials
    # if eph_ngqpt_fine != ddb_ngqpt

    nbsum_list = np.arange(gs_nband + 10, nscf_nband, step=10)
    zcut_list = [f"{zc} meV" for zc in [1, 10, 50, 100]]

    print("Performing convergenge studies wrt BZ integration, nbands in Sigma sum and zcut:")
    print("ngkpt_fine_list", ngkpt_fine_list)
    print("shiftk_fine_list", shiftk_fine_list)
    print("nbsum_list:", nbsum_list)
    print("zcut_list:", zcut_list)

    for nscf_task in nscf_work:
        # Create empty work to contain EPH tasks with this value of eph_ngqpt_fine
        eph_work = flow.new_work()
        kfine_dict = {k: nscf_task.input[k] for k in ("ngkpt", "nshiftk", "shiftk")}

        for nband_sum, zcut in itertools.product(nbsum_list, zcut_list):

            new_inp = eph_template.new_with_vars(
                nband=nband_sum,
                zcut=zcut,
                eph_ngqpt_fine=kfine_dict["ngkpt"],
                **kfine_dict,
            )

            # The EPH code requires the GS WFK, the DDB file with all perturbations
            # and the DVDB file with the DFPT potentials (already merged by the previous Flow)
            # The POT file is needed for the Sternheimer.
            deps = {
                nscf_task: "WFK",
                ddb_filepath: "DDB",
                dvdb_filepath: "DVDB",
                pot_filepath: "POT",
            }

            eph_work.register_eph_task(new_inp, deps=deps)

    return flow


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
