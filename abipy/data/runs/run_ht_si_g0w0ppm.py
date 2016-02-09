#!/usr/bin/env python
"""G0W0 corrections with the HT interface."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.data as abidata  

from abipy import abilab


def build_flow(options):
    # Init structure and pseudos.
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Initialize the flow.
    flow = abilab.Flow(workdir, manager=options.manager, remove=options.remove) 

    scf_kppa = 10
    nscf_nband = 10
    #nscf_ngkpt = [4,4,4]
    #nscf_shiftk = [0.0, 0.0, 0.0]
    ecut, ecuteps, ecutsigx = 4, 2, 3
    #scr_nband = 50
    #sigma_nband = 50

    multi = abilab.g0w0_with_ppmodel_inputs(
        structure, pseudos, 
        scf_kppa, nscf_nband, ecuteps, ecutsigx,
        ecut=ecut, pawecutdg=None,
        accuracy="normal", spin_mode="unpolarized", smearing=None,
        #ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None,
        #sigma_nband=None, gw_qprange=1):

    )
    #multi.set_vars(paral_kgb=1)

    scf_input, nscf_input, scr_input, sigma_input = multi.split_datasets()
    work = abilab.G0W0Work(scf_input, nscf_input, scr_input, sigma_input)

    flow.register_work(work)
    return flow
    

@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
