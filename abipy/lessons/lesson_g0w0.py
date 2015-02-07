#!/usr/bin/env python
"""This script shows how to compute the G0W0 corrections in silicon."""
from __future__ import division, print_function, unicode_literals

import os
import sys
import abipy.data as data  
import abipy.abilab as abilab

def help(stream=sys.stdout):
    """
    Display the tutorial text.
    """
    stream.write(__doc__)


def get_local_copy():
    """
    Copy this script to the current working dir to explore and edit
    """
    dst = os.path.basename(__file__[:-1])
    if os.path.exists(dst):
        raise RuntimeError("file %s already exists. Remove it before calling get_local_copy" % dst)
    shutil.copyfile(__file__[:-1], dst)


def make_inputs(ngkpt, paral_kgb=1):
    # Crystalline silicon
    # Calculation of the GW correction to the direct band gap in Gamma
    # Dataset 1: ground state calculation 
    # Dataset 2: NSCF calculation 
    # Dataset 3: calculation of the screening 
    # Dataset 4: Self-Energy matrix elements (GW corrections) with different values of nband

    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=6)
    inp.set_structure(data.cif_file("si.cif"))

    # This grid is the most economical, but does not contain the Gamma point.
    scf_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    dos_kmesh = dict(
        ngkpt=(6, 6, 6),
        shiftk=[0.0, 0.0, 0.0])

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap. 
    gw_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.0, 0.0, 0.0,  
                0.0, 0.5, 0.5,  
                0.5, 0.0, 0.5,  
                0.5, 0.5, 0.0]
    )

    # Global variables. gw_kmesh is used in all datasets except DATASET 1.
    ecut = 6
       
    inp.set_vars(
        ecut=ecut,
        istwfk="*1",
        paral_kgb=paral_kgb,
        gwpara=2,
    )
    #inp.set_kmesh(**gw_kmesh)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(**scf_kmesh)
    inp[1].set_vars(
        tolvrs=1e-6,
        nband=4,
    )
    inp[1].set_kmesh(**scf_kmesh)

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    inp[2].set_vars(iscf=-2,
                    tolwfr=1e-12,
                    nband=8,
                    nbdbuf=5,
                   )
    inp[2].set_kpath(ndivsm=5)

    # Dataset 3 (DOS NSCF)
    inp[3].set_vars(iscf=-2,
                    tolwfr=1e-12,
                    nband=35,
                    nbdbuf=5,
                   )
    inp[3].set_kmesh(**dos_kmesh)

    # Dataset 4 (NSCF run for GW)
    inp[4].set_vars(iscf=-2,
                    tolwfr=1e-12,
                    nband=35,
                    nbdbuf=5,
                   )
    inp[4].set_kmesh(**gw_kmesh)

    # Dataset3: Calculation of the screening.
    inp[5].set_vars(
        optdriver=3,   
        nband=25,    
        ecutwfn=ecut,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
    )
    inp[5].set_kmesh(**gw_kmesh)

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    #kptgw = [
    #     -2.50000000E-01, -2.50000000E-01,  0.00000000E+00,
    #     -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
    #      5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
    #     -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
    #      5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
    #      0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
    #  ]

    #bdgw = [1,8]

    inp[6].set_vars(
            optdriver=4,
            nband=10,      
            ecutwfn=ecut,
            ecuteps=4.0,
            ecutsigx=6.0,
            symsigma=1,
            gw_qprange=-4,  # All occupied states and 4 empty states
            #nkptgw=0,
        )
    inp[6].set_kmesh(**gw_kmesh)
    #inp[6].set_kptgw(kptgw, bdgw)

    # Return scf, bands_nscf, dos_nscf, gw_nscf, scr, sig
    return inp.split_datasets()


def build_flow():
    # Change the value of ngkpt below to perform a GW calculation with a different k-mesh.
    scf, bands_nscf, dos_nscf, gw_nscf, scr, sig = make_inputs(ngkpt=[2,2,2])

    flow = abilab.Flow(workdir="lesson_g0w0")
    work0 = abilab.BandStructureWork(scf, bands_nscf, dos_inputs=dos_nscf)
    flow.register_work(work0)

    work1 = abilab.Work()
    gw_nscf_task = work1.register_nscf_task(gw_nscf, deps={work0[0]: "DEN"})
    scr_task = work1.register_scr_task(scr, deps={gw_nscf_task: "WFK"})
    sigma_task = work1.register_sigma_task(sig, deps={gw_nscf_task: "WFK", scr_task: "SCR"})
    flow.register_work(work1)

    return flow.allocate()
    #return abilab.g0w0_flow(workdir, scf, nscf, scr, [sig1, sig2, sig3])


def analyze_flow(flow):
    sigma_task = flow[1][2]
    builder = sigma_task.get_scissors_builder()
    builder.plot_qpe_vs_e0()
    builder.build(domains_spin=[[-10, 6.02], [6.1, 20]])

    builder.plot_fit()

    bands_task = flow[0][1]
    bands_gsr_path = bands_task.outdir.has_abiext("GSR")

    builder.plot_qpbands(bands_gsr_path, title="Silicon Bands (KS and KS+scissors)")

    #dos_task = flow[0][2]
    #dos_gsr_path = dos_task.outdir.has_abiext("GSR")
    #builder.plot_qpbands(bands_gsr_path,
    #                     dos_filepath=dos_gsr_path,
    #                     title="Silicon Bands and DOS (KS and KS+scissors)")


if __name__ == "__main__":
    #flow = build_flow()
    #flow.build_and_pickle_dump()
    #flow.show_inputs()
    #flow.make_scheduler().start()

    flow = abilab.Flow.pickle_load("lesson_g0w0")
    analyze_flow(flow)

