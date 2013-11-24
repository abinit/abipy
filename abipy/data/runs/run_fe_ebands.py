#!/usr/bin/env python
"""
Calculation of the band structure of Fe with and without magnetization.
See tutorial/Input/tspin_1.in
"""
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import enable_logging, AbipyTest, MixinTest

class SpinEbandsFlowTest(AbipyTest, MixinTest):
    """
    Unit test for the flow defined in this module.  
    Users who just want to learn how to use this flow can ignore this section.
    """
    def setUp(self):
        super(SpinEbandsFlowTest, self).setUp()
        self.init_dirs()
        self.flow = build_bands_flow(workdir=self.workdir)


def make_scf_nscf_inputs(nsppol):
    inp = abilab.AbiInput(pseudos=data.pseudos("26fe.pspnc"), ndtset=2)

    # Fe normal bcc structure for test of a ferromagnetic calculation
    structure = data.structure_from_ucell("Fe-fm")
    inp.set_structure(structure)

    # Global variables
    global_vars = dict(nsppol=nsppol,
                       ecut=18,
                       nband=8,
                       occopt=3,
                       tsmear=0.01,
                       paral_kgb=0,
                    )
    if nsppol == 2:
        global_vars.update(spinat=[0.0, 0.0, 4.0])

    inp.set_variables(**global_vars)

    # Dataset 1 (GS run)
    inp[1].set_kmesh(ngkpt=[4,4,4], shiftk=[0.5,0.5,0.5])
    inp[1].set_variables(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    inp[2].set_kpath(ndivsm=4)
    inp[2].set_variables(tolwfr=1e-8)
    
    # Generate two input files for the GS and the NSCF run 
    scf_input, nscf_input = inp.split_datasets()

    return scf_input, nscf_input

def build_bands_flow(workdir="tmp_fe_ebands"):
    manager = abilab.TaskManager.from_user_config()

    flow = abilab.AbinitFlow(workdir, manager)

    # Create the task defining the calculation and run.
    for nsppol in [1,2]:
        scf_input, nscf_input = make_scf_nscf_inputs(nsppol)
        work = abilab.BandStructureWorkflow(scf_input, nscf_input)
        flow.register_work(work)

    return flow.allocate()


@enable_logging
def main():
    flow = build_bands_flow()
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    import sys
    sys.exit(main())
