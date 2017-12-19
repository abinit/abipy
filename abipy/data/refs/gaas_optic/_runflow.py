#!/usr/bin/env python
"""Optical spectra with Optic."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def build_flow(options, paral_kgb=0):
    """
    Build flow for the calculation of optical properties with optic + band structure
    along high-symmetry k-path. DDK are computed with 3 k-meshes of increasing density
    to monitor the convergece of the spectra.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    multi = abilab.MultiDataset(structure=abidata.structure_from_ucell("GaAs"),
                                pseudos=abidata.pseudos("31ga.pspnc", "33as.pspnc"), ndtset=2)

    # Usa same shifts in all tasks.
    shiftk= [[0.5, 0.5, 0.5],
             [0.5, 0.0, 0.0],
             [0.0, 0.5, 0.0],
             [0.0, 0.0, 0.5]]

    # Global variables.
    multi.set_vars(ecut=2, paral_kgb=paral_kgb)

    # Dataset 1 (GS run)
    multi[0].set_vars(tolvrs=1e-8, nband=4)
    multi[0].set_kmesh(ngkpt=[4, 4, 4], shiftk=shiftk)

    # NSCF run on k-path with large number of bands
    multi[1].set_vars(iscf=-2, nband=20, tolwfr=1.e-9)
    multi[1].set_kpath(ndivsm=10)

    # Initialize the flow.
    flow = flowtk.Flow(workdir, manager=options.manager, remove=options.remove)

    # GS to get the density + NSCF along the path.
    scf_inp, nscf_inp = multi.split_datasets()
    bands_work = flowtk.BandStructureWork(scf_inp, nscf_inp)
    flow.register_work(bands_work)

    # Build OpticInput used to compute optical properties.
    optic_input = abilab.OpticInput(
        broadening=0.002,          # Value of the smearing factor, in Hartree
        domega=0.0003,             # Frequency mesh.
        maxomega=0.3,
        scissor=0.000,             # Scissor shift if needed, in Hartree
        tolerance=0.002,           # Tolerance on closeness of singularities (in Hartree)
        num_lin_comp=2,            # Number of components of linear optic tensor to be computed
        lin_comp=(11, 33),         # Linear coefficients to be computed (x=1, y=2, z=3)
        num_nonlin_comp=2,         # Number of components of nonlinear optic tensor to be computed
        nonlin_comp=(123, 222),    # Non-linear coefficients to be computed
        num_linel_comp=1,          # Number of components of LEO tensor to be computed
        linel_comp=(123),          # Non-linear coefficients to be computed
    )

    # ddk_nband is fixed here, in principle it depends on nelect and the frequency range in chi(w).
    ddk_nband = 20

    # Perform converge study wrt ngkpt (shiftk is constant).
    ngkpt_convergence = [[4, 4, 4], [8, 8, 8], [12, 12, 12]]

    for ddk_ngkpt in ngkpt_convergence:
        # Build work for NSCF from DEN produced by the first GS task + 3 DDKs.
        # All tasks use more bands and a denser k-mesh defined by ddk_ngkpt.
        ddks_work = flowtk.NscfDdksWork.from_scf_task(bands_work[0], ddk_ngkpt, shiftk, ddk_nband)
        flow.register_work(ddks_work)

        # Build optic task to compute chi with this value of ddk_ngkpt.
        optic_task = flowtk.OpticTask(optic_input, nscf_node=ddks_work.task_with_ks_energies,
                                      ddk_nodes=ddks_work.ddk_tasks, use_ddknc=False)
        ddks_work.register_task(optic_task)

    return flow


@flowtk.flow_main
def main(options):
    return build_flow(options)


if __name__ == "__main__":
    retcode = main()
    if retcode != 0: sys.exit(retcode)

    rename_table = [
        #  src, dest
        #"_runflow/w0/t1/outdata/out_GSR.nc",
        ("_runflow/w1/t4/outdata/out_OPTIC.nc", "gaas_444_OPTIC.nc"),
        ("_runflow/w2/t4/outdata/out_OPTIC.nc", "gaas_888_OPTIC.nc"),
        ("_runflow/w3/t4/outdata/out_OPTIC.nc", "gaas_121212_OPTIC.nc"),
        ("_runflow/w1/t1/outdata/out_DDK.nc", "gaas_444_dir1_DDK.nc"),
        ("_runflow/w1/t2/outdata/out_DDK.nc", "gaas_444_dir2_DDK.nc"),
        ("_runflow/w1/t3/outdata/out_DDK.nc", "gaas_444_dir3_DDK.nc"),
    ]
    import shutil
    for old, new in rename_table:
        shutil.copyfile(old, new)

    #shutil.rmtree("_runflow")
    sys.exit(0)
