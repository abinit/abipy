#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab

from pymatgen.io.abinitio.task import RunMode
from abipy.data.runs import RunManager

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

    print(scf_inp)

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

    print(nscf_inp)

    # Create the task defining the calculation and run.
    runmode = RunMode.sequential()
    manager = RunManager()
    
    # Initialize the workflow.
    work = abilab.Workflow(manager.workdir, runmode)

    # Register the input for the SCF calculation and receive an object 
    # that describes this node of the worflow.
    scf_link = work.register_input(scf_inp)

    # Register the input for the NSCF calculation and tell the workflow
    # that this step depens on the SCF run (requires the DEN file produced in the SCF run).
    work.register_input(nscf_inp, links=scf_link.produces_exts("_DEN"))

    work.show_inputs()

    manager.set_work_and_run(work)

    if manager.retcode != 0:
        return manager.retcode

    # Remove all files except those matching these regular expressions.
    work.rmtree(exclude_wildcard="*.abi|*.abo")

    #work.rename("out_DS1_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
    #work.rename("out_DS1_DEN-etsf.nc", "si_DEN-etsf.nc")
    #work.rename("out_DS1_GSR.nc", "si_scf_GSR.nc")

    #work.rename("out_DS2_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")
    #work.rename("out_DS2_GSR.nc", "si_nscf_GSR.nc")
    #work.remove_files("out_DS2_DEN-etsf.nc")

    manager.finalize()
    return manager.retcode 

if __name__ == "__main__":
    import sys
    sys.exit(main())
