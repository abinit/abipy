#!/usr/bin/env python
"""This script shows how to compute the G0W0 corrections in silicon."""
from __future__ import division, print_function

import os
import sys
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import AbipyTest, MixinTest


class SiG0W0FlowTest(AbipyTest, MixinTest):
    """
    Unit test for the flow defined in this module.  
    Users who just want to learn how to use this flow can ignore this section.
    """
    def setUp(self):
        super(SiG0W0FlowTest, self).setUp()
        self.init_dirs()
        self.flow = build_flow()

    #task.rename("out_DS4_SIGRES.nc", "si_g0w0ppm_nband10_SIGRES.nc")
    #task.rename("out_DS5_SIGRES.nc", "si_g0w0ppm_nband20_SIGRES.nc")
    #task.rename("out_DS6_SIGRES.nc", "si_g0w0ppm_nband30_SIGRES.nc")


def make_inputs(ngkpt):
    # Crystalline silicon
    # Calculation of the GW correction to the direct band gap in Gamma
    # Dataset 1: ground state calculation 
    # Dataset 2: NSCF calculation 
    # Dataset 3: calculation of the screening 
    # Dataset 4-5-6: Self-Energy matrix elements (GW corrections) with different values of nband

    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=6)

    inp.set_structure_from_file(data.cif_file("si.cif"))

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
       
    inp.set_variables(
        ecut=ecut,
        timopt=-1,
        istwfk="*1",
        paral_kgb=0,
        gwpara=2,
    )
    inp.set_kmesh(**gw_kmesh)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(**scf_kmesh)
    inp[1].set_variables(
        tolvrs=1e-6,
        nband=4,
    )

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    inp[2].set_variables(iscf=-2,
                         tolwfr=1e-12,
                         nband=35,
                         nbdbuf=5,
                        )

    # Dataset3: Calculation of the screening.
    inp[3].set_variables(
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
        inp[4+idx].set_variables(
            optdriver=4,
            nband=nband,      
            ecutwfn=ecut,
            ecuteps=4.0,
            ecutsigx=6.0,
            symsigma=1,
        )
        inp[4+idx].set_kptgw(kptgw, bdgw)

    return inp.split_datasets()


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed)
    workdir = os.path.basename(os.path.abspath(__file__).replace(".py", "")) if not options.workdir else options.workdir

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else options.manager

    # Change the value of ngkpt below to perform a GW calculation with a different k-mesh.
    scf, nscf, scr, sig1, sig2, sig3 = make_inputs(ngkpt=[2,2,2])

    flow = abilab.g0w0_flow(workdir, manager, scf, nscf, scr, [sig1, sig2, sig3])
    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    return flow.build_and_pickle_dump()

if __name__ == "__main__":
    sys.exit(main())
