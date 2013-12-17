#!/usr/bin/env python
"""
This example shows how to compute the SCR file by splitting the calculation
over q-points with the input variables nqptdm and qptdm.
"""
from __future__ import division, print_function

import sys
import os
import abipy.abilab as abilab
import abipy.data as data  

from abipy.data.runs import AbipyTest, MixinTest
from abipy.data.runs.qptdm_workflow import *


class QptdmFlowTest(AbipyTest, MixinTest):
    """
    Unit test for the flow defined in this module.  
    Users who just want to learn how to use this flow can ignore this section.
    """
    def setUp(self):
        super(QptdmFlowTest, self).setUp()
        self.init_dirs()
        self.flow = qptdm_flow()


def all_inputs():
    structure = abilab.Structure.from_file(data.cif_file("si.cif"))
    pseudos = data.pseudos("14si.pspnc")

    ecut = ecutwfn = 6

    global_vars = dict(
        ecut=ecut,
        timopt=-1,
        istwfk="*1",
        paral_kgb=0,
    )

    inp = abilab.AbiInput(pseudos=pseudos, ndtset=4)
    print("pseudos",inp.pseudos)
    inp.set_structure(structure)
    inp.set_variables(**global_vars)

    gs, nscf, scr, sigma = inp[1:]

    # This grid is the most economical, but does not contain the Gamma point.
    gs_kmesh = dict(
        ngkpt=[2,2,2],
        shiftk=[0.5, 0.5, 0.5,
                0.5, 0.0, 0.0,
                0.0, 0.5, 0.0,
                0.0, 0.0, 0.5]
    )

    # This grid contains the Gamma point, which is the point at which
    # we will compute the (direct) band gap. 
    gw_kmesh = dict(
        ngkpt=[2,2,2],
        shiftk=[0.0, 0.0, 0.0,  
                0.0, 0.5, 0.5,  
                0.5, 0.0, 0.5,  
                0.5, 0.5, 0.0]
    )

    # Dataset 1 (GS run)
    gs.set_kmesh(**gs_kmesh)
    gs.set_variables(tolvrs=1e-6,
                     nband=4,
                    )

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    nscf.set_kmesh(**gw_kmesh)

    nscf.set_variables(iscf=-2,
                       tolwfr=1e-12,
                       nband=35,
                       nbdbuf=5,
                       )

    # Dataset3: Calculation of the screening.
    scr.set_kmesh(**gw_kmesh)

    scr.set_variables(
        optdriver=3,   
        nband=25,    
        ecutwfn=ecutwfn,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
        ppmfrq="16.7 eV",
    )

    # Dataset4: Calculation of the Self-Energy matrix elements (GW corrections)
    sigma.set_kmesh(**gw_kmesh)

    sigma.set_variables(
            optdriver=4,
            nband=25,      
            ecutwfn=ecutwfn,
            ecuteps=4.0,
            ecutsigx=6.0,
            #symsigma=1,
            gwcalctyp=20,
        )

    kptgw = [
         -2.50000000E-01, -2.50000000E-01,  0.00000000E+00,
         -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
          5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
         -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
          5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
          0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
      ]

    bdgw = [1,8]

    sigma.set_kptgw(kptgw, bdgw)

    return inp.split_datasets()


def gw_flow(options):
    #workdir = "GW"
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir 
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 
                                                                                                                         
    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else options.manager

    gs, nscf, scr_input, sigma_input = all_inputs()
                                                                        
    manager = abilab.TaskManager.from_user_config()

    flow = g0w0_flow_with_qptdm(workdir, manager, gs, nscf, scr_input, sigma_input)

    #from pymatgen.io.abinitio.tasks import G_Task
    #flow[2][0].__class__ = G_Task

    return flow
    
def qptdm_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir 
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else options.manager

    gs, nscf, scr_input, sigma_input = all_inputs()

    return g0w0_flow_with_qptdm(workdir, manager, gs, nscf, scr_input, sigma_input)


@abilab.flow_main
def main(options):
    # GW Works
    #flow = gw_flow(options)
                        
    # QPTDM
    flow = qptdm_flow(options)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
