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
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Preparatory run for E-PH calculations.
    # The sequence of datasets makes the ground states and
    # all of the independent perturbations of the single Al atom
    # for the irreducible qpoints in a 4x4x4 grid.
    # Note that the q-point grid must be a sub-grid of the k-point grid (here 8x8x8)
    pseudos = abidata.pseudos("Al.oncvpsp")

    structure = abilab.Structure.from_abivars(
        acell=3*[7.5],
        rprim=[0.0, 0.5, 0.5,
               0.5, 0.0, 0.5,
               0.5, 0.5, 0.0],
        typat=1,
        xred=[0.0, 0.0, 0.0],
        ntypat=1,
        znucl=13,
    )

    gs_inp = abilab.AbinitInput(structure, pseudos)

    gs_inp.set_vars(
        istwfk="*1",
        ecut=12.0,
        nband=5,
        occopt=7,    # include metallic occupation function with a small smearing
        tsmear=0.04,
        tolvrs=1e-7,
        timopt=-1,
    )

    # The kpoint grid is minimalistic to keep the calculation manageable.
    gs_inp.set_kmesh(
        ngkpt=[12, 12, 12],
        shiftk=[0.0, 0.0, 0.0],
        #kptopt=3,
    )

    flow = flowtk.Flow(workdir, manager=options.manager, remove=options.remove)

    # GS run to get the WFK
    work0 = flow.register_task(gs_inp, task_class=flowtk.ScfTask)

    # Phonon work with 4x4x4 q-mesh
    ddb_ngqpt = [4, 4, 4]
    ph_work = flowtk.PhononWork.from_scf_task(work0[0], ddb_ngqpt, is_ngqpt=True)
    flow.register_work(ph_work)

    # Build input file for E-PH run.
    eph_inp = gs_inp.new_with_vars(
        optdriver=7,
        ddb_ngqpt=ddb_ngqpt,       # q-mesh used to produce the DDB file (must be consistent with DDB data)
        eph_intmeth=2,             # Tetra method
        eph_fsewin="0.8 eV",       # Energy window around Ef
        eph_mustar=0.12,           # mustar parameter
    )

    # Set q-path for phonons and phonon linewidths.
    eph_inp.set_qpath(10)

    # Set q-mesh for phonons DOS and a2F(w)
    eph_inp.set_phdos_qmesh(nqsmall=16, method="tetra")

    # EPH part requires the GS WFK, the DDB file with all perturbations
    # and the database of DFPT potentials (already merged by PhononWork)
    eph_work = flow.register_work(flowtk.Work())
    deps = {work0[0]: "WFK", ph_work: ["DDB", "DVDB"]}

    eph_work.register_eph_task(eph_inp, deps=deps)

    # Activate Fourier interpolation of DFPT potentials.
    eph_work.register_eph_task(eph_inp.new_with_vars(eph_ngqpt_fine=[6, 6, 6]), deps=deps)
    eph_work.register_eph_task(eph_inp.new_with_vars(eph_ngqpt_fine=[12, 12, 12]), deps=deps)

    flow.allocate()

    # EPH does not support autoparal (yet)
    #for eph_task in eph_work:
    #    eph_task.with_fixed_mpi_omp(1, 1)

    return flow


@flowtk.flow_main
def main(options):
    return build_flow(options)


if __name__ == "__main__":
    retcode = main()
    if retcode != 0: sys.exit(retcode)

    rename_table = [
        #  src, dest
        ("_runflow/w2/t0/outdata/out_A2F.nc", "al_888k_161616q_A2F.nc"),
        ("_runflow/w2/t1/outdata/out_A2F.nc", "al_888k_161616q_vinterp_A2F.nc"),
	("_runflow/w2/t2/outdata/out_PHDOS.nc", "al_161616q_PHDOS.nc"),
    ]
    import shutil
    for old, new in rename_table:
        shutil.copyfile(old, new)
    #shutil.rmtree("_runflow")
    sys.exit(0)
