#!/usr/bin/env python
r"""
G0W0 convergence study
======================

G0W0 convergence study wrt ecuteps and the number of bands in W.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import numpy as np

import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk

# system aopt SR DFT FR DFT SR G0W0 FR G0W0 aexpt157 SR DFT SR G0W0 expt.157
alatt = {
"Al-Sb": 6.23, # 1.22 1.02 1.89 1.65 6.136 1.20 2.23 1.60
"Cd-S": 5.93,  # 1.02 1.01 1.61 1.59 5.832 1.14 1.77 2.42
"Cd-Te": 6.63, # 0.56 0.29 1.07 0.76 6.477 0.76 1.29 1.44
"Ga-As": 5.75, # 0.15 0.01 0.25 0.13 5.653 0.50 0.62 1.35
"Ga-P": 5.51,  # 1.54 1.51 1.95 1.91 5.450 1.63 2.17 2.24
"In-P": 5.96,  # 0.42 0.39 0.61 0.57 5.869 0.68 0.90 1.27
"Si-Si": 5.48,   # 0.72 0.71 1.38 1.36 5.431 0.67 1.29 1.12
"Zn-S": 5.45,  # 1.99 1.97 3.08 3.06 5.409 2.09 3.26 3.54
"Zn-Se": 5.73, # 1.13 1.01 1.89 1.74 5.668 1.26 2.06 2.58
"Zn-Te": 6.18, # 1.08 0.79 1.87 1.49 6.101 1.24 2.08 2.26
}


def make_inputs(nspinor, paral_kgb=1):
    """
    Returns a tuple of 4 input files for SCF, NSCF, SCR, SIGMA calculations.
    These files are then used as templates for the convergence study
    wrt ecuteps and the number of bands in W.
    """
    formula = "Si-Si"
    species = formula.split("-")
    assert len(species) == 2
    a = alatt[formula]
    pseudos = "Si_r.psp8"
    #pseudos = abilab.PseudoTable.from_dir()
    structure = abilab.Structure.fcc(a, species)
    structure = abidata.structure_from_ucell("Si")
    multi = abilab.MultiDataset(structure, pseudos, ndtset=5)
    #pseudos = multi[0].pseudos

    #kptopt=1 if nspinor == 1 else 4
    kptopt=4
    ngkpt = [4, 4, 4]
    shiftk = [0, 0, 0]

    ecut = 4
    ecuteps = 4
    nband_gs = 6 * nspinor
    nband_mbpt = 10 * nspinor
    nbdbuf = 8

    # non-collinear GS (activate SOC allow for magnetic solution although Si is not)
    global_vars = dict(
        nspinor=nspinor,
        nspden=1 if nspinor == 1 else 4,
        so_psp="*0" if nspinor == 1 else "*1",
        ecut=ecut,
        istwfk="*1",
        paral_kgb=paral_kgb,
        #nsym=1,
        #iomode=1,
    )

    multi.set_vars(global_vars)
    multi.set_kmesh(ngkpt=ngkpt, shiftk=shiftk, kptopt=kptopt)

    # SCF
    multi[0].set_vars(
        nband=nband_gs,
        tolvrs=1.e-10,
    )

    # NSCF on k-path
    multi[1].set_vars(
        nband=nband_gs + nbdbuf,
        nbdbuf=nbdbuf,
        tolwfr=1.e-18,
        iscf=-2,
    )
    multi[1].set_kpath(ndivsm=10)

    # NSCF on k-mesh
    multi[2].set_vars(
        nband=nband_gs + nbdbuf,
        nbdbuf=nbdbuf,
        tolwfr=1.e-18,
        iscf=-2
    )

    # SCR
    multi[3].set_vars(
        optdriver=3,
        ecutwfn=ecut,
        ecuteps=ecuteps,
        nband=nband_mbpt,
        symchi=1,
        awtr=1 if kptopt == 1 else 0,
        gwpara=2 if kptopt == 1 else 1,
        inclvkb=0,
    )

    # SIGMA
    multi[4].set_vars(
        optdriver=4,
        nband=nband_mbpt,
        ecutwfn=ecut,
        ecuteps=ecuteps,
        ecutsigx=ecut,
        #ecutsigx=(4*ecut), ! This is problematic
        symsigma=1,
        gw_qprange=0,
        )

    #multi[4].set_kptgw(kptgw=[[0,0,0], [0.5, 0, 0]], bdgw=[1, 8])

    return multi.split_datasets()


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    for nspinor in (1, 2):
    #for nspinor in (2,):
        # Get our templates
        scf_inp, bands_inp, nscf_inp, scr_inp, sig_inp = make_inputs(nspinor)

        # Band structure work to produce the WFK file
        bands_work = flowtk.BandStructureWork(scf_inp, bands_inp, dos_inputs=[nscf_inp])
        flow.register_work(bands_work)

        # Build a work made of two SCR runs with different value of nband
        gw_work = flowtk.Work()
        scr_task = gw_work.register_scr_task(scr_inp, deps={bands_work[2]: "WFK"})
        gw_work.register_sigma_task(sig_inp, deps={bands_work[2]: "WFK", scr_task: "SCR"})
        flow.register_work(gw_work)

    return flow


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)


if __name__=="__main__":
    sys.exit(main())
