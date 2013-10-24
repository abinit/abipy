#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import Tester

def main():
    structure = abilab.Structure.from_file(data.cif_file("si.cif"))

    # Global variables
    global_vars = dict(
        ecut=6,
        #accesswff=3,
        #istwfk = "*1",
    )

    # GS run
    scf_inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))
    scf_inp.set_structure(structure)
    scf_inp.set_variables(**global_vars)
    scf_inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0])
    scf_inp.nband  = 8
    scf_inp.tolvrs = 1e-6

    # NSCF run
    kptbounds = [
        [0.5, 0.0,  0.0],  # L point
        [0.0, 0.0,  0.0],  # Gamma point
        [0.0, 0.5,  0.5],  # X point
    ]

    nscf_inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"))

    nscf_inp.set_structure(structure)
    nscf_inp.set_kpath(ndivsm=6, kptbounds=kptbounds)
    nscf_inp.set_variables(**global_vars)
    #nscf_inp.set_kpath(ndivsm=5)
    nscf_inp.tolwfr = 1e-12

    # Create the task defining the calculation and run.
    tester = Tester()
    manager = abilab.TaskManager.from_user_config()
    
    # Initialize the workflow.
    work = abilab.Workflow(tester.workdir, manager)

    # Register the input for the SCF calculation. 
    # scf_task is the object that describes this node of the workflow.
    scf_task = work.register(scf_inp)

    # Register the input for the NSCF calculation and tell the workflow
    # that this step depens on the SCF run 
    # Iin this case, the nscf run requires the DEN file produced in the SCF run.
    work.register(nscf_inp, links={scf_task: "DEN"})

    #work.build()
    #tester.set_work_and_run(work)

    work.show_inputs()
    #from abipy.gui.wxapps import wxapp_showfiles
    #wxapp_showfiles(dirpath=work.workdir, walk=True, wildcard="*.abo").MainLoop()

    if tester.retcode != 0:
        return tester.retcode

    # Remove all files except those matching these regular expressions.
    work.rmtree(exclude_wildcard="*.abi|*.abo")

    #work.rename("out_DS1_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
    #work.rename("out_DS1_DEN-etsf.nc", "si_DEN-etsf.nc")
    #work.rename("out_DS1_GSR.nc", "si_scf_GSR.nc")

    #work.rename("out_DS2_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")
    #work.rename("out_DS2_GSR.nc", "si_nscf_GSR.nc")
    #work.remove_files("out_DS2_DEN-etsf.nc")

    tester.finalize()
    return tester.retcode 

if __name__ == "__main__":
    import sys
    sys.exit(main())
