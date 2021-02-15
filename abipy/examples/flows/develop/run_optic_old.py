#!/usr/bin/env python
r"""
Optic Flow
==========

Optical spectra with Optic.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def build_flow(options, paral_kgb=0):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    multi = abilab.MultiDataset(structure=abidata.structure_from_ucell("GaAs"),
                                pseudos=abidata.pseudos("31ga.pspnc", "33as.pspnc"), ndtset=5)

    # Global variables
    kmesh = dict(ngkpt=[4, 4, 4],
                 nshiftk=4,
                 shiftk=[[0.5, 0.5, 0.5],
                         [0.5, 0.0, 0.0],
                         [0.0, 0.5, 0.0],
                         [0.0, 0.0, 0.5]]
                )

    global_vars = dict(ecut=2, paral_kgb=paral_kgb)
    global_vars.update(kmesh)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_vars(
        tolvrs=1e-6,
        nband=4,
    )

    # NSCF run with large number of bands, and points in the the full BZ
    multi[1].set_vars(
        iscf=-2,
        nband=20,
        nstep=25,
        kptopt=1,
        tolwfr=1.e-9,
        #kptopt=3,
    )

    # Fourth dataset : ddk response function along axis 1
    # Fifth dataset : ddk response function along axis 2
    # Sixth dataset : ddk response function along axis 3
    for idir in range(3):
        rfdir = 3 * [0]
        rfdir[idir] = 1

        multi[2+idir].set_vars(
            iscf=-3,
            nband=20,
            nstep=1,
            nline=0,
            prtwf=3,
            kptopt=3,
            nqpt=1,
            qpt=[0.0, 0.0, 0.0],
            rfdir=rfdir,
            rfelfd=2,
            tolwfr=1.e-9,
        )

    scf_inp, nscf_inp, ddk1, ddk2, ddk3 = multi.split_datasets()

    # Initialize the flow.
    flow = flowtk.Flow(options.workdir, manager=options.manager)

    bands_work = flowtk.BandStructureWork(scf_inp, nscf_inp)
    flow.register_work(bands_work)

    ddk_work = flowtk.Work()
    for inp in [ddk1, ddk2, ddk3]:
        ddk_work.register_ddk_task(inp, deps={bands_work.nscf_task: "WFK"})

    flow.register_work(ddk_work)

    # Optic does not support MPI with ncpus > 1.
    optic_input = abilab.OpticInput(
        broadening=0.002,          # Value of the smearing factor, in Hartree
        domega=0.0003,             # Frequency mesh.
        maxomega=0.3,
        scissor=0.000,             # Scissor shift if needed, in Hartree
        tolerance=0.002,           # Tolerance on closeness of singularities (in Hartree)
        num_lin_comp=1,            # Number of components of linear optic tensor to be computed
        lin_comp=11,               # Linear coefficients to be computed (x=1, y=2, z=3)
        num_nonlin_comp=2,         # Number of components of nonlinear optic tensor to be computed
        nonlin_comp=(123, 222),    # Non-linear coefficients to be computed
    )

    # TODO
    # Check is the order of the 1WF files is relevant. Can we use DDK files ordered
    # in an arbitrary way or do we have to pass (x,y,z)?
    optic_task = flowtk.OpticTask(optic_input, nscf_node=bands_work.nscf_task, ddk_nodes=ddk_work)
    flow.register_task(optic_task)

    return flow


def optic_flow_from_files():
    # Optic does not support MPI with ncpus > 1.
    manager = flowtk.TaskManager.from_user_config()
    manager.set_mpi_procs(1)

    flow = flowtk.Flow(workdir="OPTIC_FROM_FILE", manager=manager)

    ddk_nodes = [
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_0/outdata/out_1WF",
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_1/outdata/out_1WF",
        "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_1/task_2/outdata/out_1WF",
    ]
    nscf_node = "/Users/gmatteo/Coding/abipy/abipy/data/runs/OPTIC/work_0/task_1/outdata/out_WFK"

    optic_task = flowtk.OpticTask(optic_input, nscf_node=nscf_node, ddk_nodes=ddk_nodes)
    flow.register_task(optic_task)
    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).plot_networkx(tight_layout=True)


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
