#!/usr/bin/env python
"""Electron-phonon calculations."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import sys
import numpy as np
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk


def build_flow(options):
    """
    C in diamond structure. Very rough q-point mesh, low ecut, completely unconverged.
    The flow compute the ground state density a WFK file on a k-mesh used for DFPT phonons
    a WFK file on a k-mesh with empty states used for the self-energy.
    Then all of the independent perturbations for the irreducible qpoints in a 4x4x4 grid are obtained.
    Note that the q-point grid must be a sub-grid of the k-point grid (here 8x8x8)
    Finally, we enter the EPH driver to compute the Fan-Migdal self-energy.
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Define structure explicitly.
    structure = abilab.Structure.from_abivars(
        acell=3*[6.70346805],
        rprim=[0.0, 0.5, 0.5,
               0.5, 0.0, 0.5,
               0.5, 0.5, 0.0],
        typat=[1, 1],
        xred=[0.0, 0.0, 0.0, 0.25, 0.25, 0.25],
        ntypat=1,
        znucl=6,
    )

    gs_inp = abilab.AbinitInput(structure, pseudos=abidata.pseudos("C.oncvpsp"))

    gs_inp.set_vars(
        istwfk="*1",
        ecut=12.0,
        nband=4,
        tolvrs=1e-10,
    )

    # The kpoint grid is minimalistic to keep the calculation manageable.
    gs_inp.set_kmesh(
        ngkpt=[4, 4, 4],
        shiftk=[0.0, 0.0, 0.0],
        #kptopt=3,
    )

    # NSCF run with k-path (just for plotting purpose)
    nscf_kpath_inp = gs_inp.new_with_vars(
        nband=8,
        tolwfr=1e-16,
    )
    #nscf_kpath_inp.pop_vars(["tolvrs"])
    nscf_kpath_inp.set_kpath(ndivsm=10)

    # NSCF run with k-mesh to get WFK with empty states.
    nscf_empty_kmesh_inp = gs_inp.new_with_vars(
        nband=54,
        nbdbuf=4,
        tolwfr=1e-16,
        iscf=-2,
    )

    flow = flowtk.Flow(workdir, manager=options.manager, remove=options.remove)

    # GS run to get the WFK
    work0 = flowtk.BandStructureWork(gs_inp, nscf_kpath_inp, dos_inputs=[nscf_empty_kmesh_inp])
    flow.register_work(work0)

    # Phonon work with 4x4x4 q-mesh
    ddb_ngqpt = [2, 2, 2]
    ph_work = flowtk.PhononWork.from_scf_task(work0[0], ddb_ngqpt, is_ngqpt=True)
    flow.register_work(ph_work)

    # Build input file for E-PH run. See v8/Input/t44.in
    eph_inp = gs_inp.new_with_vars(
        optdriver=7,               # EPH driver.
        eph_task=4,                # For electronic self-energy due to phonon
        nband=54,
        ddb_ngqpt=ddb_ngqpt,       # q-mesh used to produce the DDB file (must be consistent with DDB data)
        symsigma=0,
        gw_qprange=2,
        #eph_intmeth=2,            # Tetra method
        #gw_qprange -2
    )

    # Set q-path for phonons and phonon linewidths.
    eph_inp.set_qpath(10)

    # Set q-mesh for phonons DOS.
    eph_inp.set_phdos_qmesh(nqsmall=16, method="tetra")

    # EPH part requires the GS WFK, the DDB file with all perturbations
    # and the database of DFPT potentials (already merged by PhononWork)
    eph_work = flow.register_work(flowtk.Work())
    deps = {work0[2]: "WFK", ph_work: ["DDB", "DVDB"]}

    eph_work.register_eph_task(eph_inp, deps=deps)

    # Activate Fourier interpolation of DFPT potentials.
    eph_work.register_eph_task(eph_inp.new_with_vars(eph_ngqpt_fine=[4, 4, 4]), deps=deps)
    #eph_work.register_eph_task(eph_inp.new_with_vars(eph_ngqpt_fine=[12, 12, 12]), deps=deps)

    flow.allocate()

    return flow


@flowtk.flow_main
def main(options):
    return build_flow(options)


if __name__ == "__main__":
    retcode = main()
    if retcode != 0: sys.exit(retcode)

    rename_table = [
        # src, dest
        #("_runflow/w0/t1/outdata/out_GSR.nc", "diamond_kpath_GSR.nc"),
        #("_runflow/w1/outdata/out_DDB", "diamond_444q_DDB"),
        ("_runflow/w2/t1/outdata/out_SIGEPH.nc", "diamond_444q_full_SIGEPH.nc"),
    ]

    import shutil
    for old, new in rename_table:
        shutil.copyfile(old, new)
    #shutil.rmtree("_runflow")
    sys.exit(0)
