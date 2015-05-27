#!/usr/bin/env python
"""This script shows how to compute the G0W0 corrections in silicon."""
from __future__ import division, print_function, unicode_literals

import os
import sys
import abipy.data as data  
import abipy.abilab as abilab

def make_inputs(ngkpt, paral_kgb=1):
    # Crystalline silicon
    # Calculation of the GW correction to the direct band gap in Gamma
    # Dataset 1: ground state calculation 
    # Dataset 2: NSCF calculation 
    # Dataset 3: calculation of the screening 
    # Dataset 4-5-6: Self-Energy matrix elements (GW corrections) with different values of nband

    multi = abilab.MultiDataset(structure=data.cif_file("si.cif"), pseudos=data.pseudos("14si.pspnc"), ndtset=6)

    # This grid is the most economical, but does not contain the Gamma point.
    scf_kmesh = dict(
        ngkpt=ngkpt,
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

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
       
    multi.set_vars(
        ecut=ecut,
        timopt=-1,
        istwfk="*1",
        paral_kgb=paral_kgb,
        gwpara=2,
    )
    multi.set_kmesh(**gw_kmesh)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(**scf_kmesh)
    multi[0].set_vars(
        tolvrs=1e-6,
        nband=4,
    )

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax multi[1]
    multi[1].set_vars(iscf=-2,
                    tolwfr=1e-12,
                    nband=35,
                    nbdbuf=5,
                   )

    # Dataset3: Calculation of the screening.
    multi[2].set_vars(
        optdriver=3,   
        nband=25,    
        ecutwfn=ecut,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
        ppmfrq="16.7 eV",
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    kptgw = [
         -2.50000000E-01, -2.50000000E-01,  0.00000000E+00,
         -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
          5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
         -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
          5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
          0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
      ]

    bdgw = [1,8]

    for idx, nband in enumerate([10, 20, 30]):
        multi[3+idx].set_vars(
            optdriver=4,
            nband=nband,      
            ecutwfn=ecut,
            ecuteps=4.0,
            ecutsigx=6.0,
            symsigma=1,
            #gw_qprange=0,
            #nkptgw=0,
        )
        multi[3+idx].set_kptgw(kptgw, bdgw)

    return multi.split_datasets()


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Change the value of ngkpt below to perform a GW calculation with a different k-mesh.
    scf, nscf, scr, sig1, sig2, sig3 = make_inputs(ngkpt=[2,2,2])

    return abilab.g0w0_flow(workdir, scf, nscf, scr, [sig1, sig2, sig3], manager=options.manager)


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
