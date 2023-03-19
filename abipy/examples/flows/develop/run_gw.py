#!/usr/bin/env python
r"""
G0W0 convergence study
======================

G0W0 convergence study wrt ecuteps and the number of bands in W.
"""
import sys
import os
import numpy as np

import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk


def build_flow(options):

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    # Build the SCF input.
    scf_inp = abilab.AbinitInput(structure=abidata.structure_from_ucell("SiC"),
                                 pseudos=abidata.pseudos("14si.pspnc", "6c.pspnc"))

    ngkpt = [4, 4, 4]
    shiftk = [0, 0, 0]

    print("Number of valence electrons in the unit cell:", scf_inp.num_valence_electrons)

    scf_inp.set_vars(
        ecut=12,
        nband=10,
        ngkpt=ngkpt,
        shiftk=shiftk,
        tolvrs=1.e-4,
        #paral_kgb=1,
        #iomode=3,
    )

    #ebands_inp = scf_inp.make_ebands_input(ndivsm=15)

    # Build NSCF input on a k-mesh including empty states.
    nscf_inp = scf_inp.new_with_vars(
        nband=25,
        #nbdbuf=10
        tolwfr=1.e-8,
        iscf=-2
    )

    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    # Band structure work to produce the WFK file with a NSCF run + empty states.
    bands_work = flowtk.BandStructureWork(scf_inp, nscf_inp) # dos_inputs=[nscf_inp])
    flow.register_work(bands_work)

    # Compute screening by splitting all the q-points.
    ecuteps_list = np.arange(2, 8, 2)
    #ecuteps_list = [2]
    max_ecuteps = max(ecuteps_list)

    for ecuteps in ecuteps_list:
        scr_inp = nscf_inp.new_with_vars(
            optdriver=3,
            ecuteps=ecuteps,
            nband=20,
        )

        scr_work = flowtk.ScreeningWork.from_nscf_task(bands_work.nscf_task, scr_inp)

        # Ii you alredy have a WFK file and you want to skip the SCF + NSCF part
        # build the work using:

        #scr_work = ScreeningWork.from_wfk_file(wfk_filepath, scr_inp)

        scr_work.set_readme(f"SCR computation with ecuteps={ecuteps}")
        flow.register_work(scr_work)

    # Do a convergence study wrt ecuteps, each work is connected to a
    # different SCR file computed with a different value of nband.
    # Build a list of sigma inputs with different ecuteps

    sigma_work = flow.new_work()

    # Input file for SIGMA
    sigma_inp = scr_inp.new_with_vars(
        optdriver=4,
        nband=20,
        ecuteps=max_ecuteps,
        ecutsigx=scf_inp["ecut"],
        #ecutsigx=(3.5 * scf_input["ecut"]), ! This is problematic
        #gw_qprange
    )

    sigma_inp.set_kptgw(kptgw=[[0, 0, 0], [0.5, 0, 0]], bdgw=[1, 8])

    sigma_work.register_sigma_task(sigma_inp, deps={bands_work.nscf_task: "WFK", scr_work: "SCR"})

    #for scr_task in scr_work:
    #    sigma_conv = flowtk.SigmaConvWork(wfk_node=bands.nscf_task, scr_node=scr_task,
    #                                      sigma_inputs=sigma_inputs)
    #    flow.register_work(sigma_conv)

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
