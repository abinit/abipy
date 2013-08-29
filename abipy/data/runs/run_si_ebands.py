#!/usr/bin/env python
from __future__ import division, print_function

import sys
import os
import abipy.data as data  
import abipy.abilab as abilab

from pymatgen.io.abinitio.task import RunMode
from abipy.data.runs import RunManager

def main():
    structure = abilab.Structure.from_file(data.cif_file("si.cif"))

    inp = abilab.AbiInput(pseudos=data.pseudos("14si.pspnc"), ndtset=2)
    inp.set_structure_from_file(data.cif_file("si.cif"))

    # Global variables
    inp.ecut = 6
    inp.nband  = 8
    inp.timopt = -1
    inp.accesswff = 3
    inp.istwfk = "*1"

    # Dataset 1 (GS run)
    inp.set_kmesh(ngkpt=[8,8,8], shiftk=[0,0,0], dtset=1)
    inp.tolvrs1 = 1e-6

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0,  0.0],  # L point
        [0.0, 0.0,  0.0],  # Gamma point
        [0.0, 0.5,  0.5],  # X point
    ]

    inp.set_kpath(ndivsm=6, kptbounds=kptbounds, dtset=2)
    #inp.set_kpath(ndivsm=5, dtset=2)
    inp.tolwfr2 = 1e-12
    inp.getden2 = -1

    print(inp)

    # Create the task defining the calculation and run.
    runmode = RunMode.sequential()
    manager = RunManager()

    task = abilab.AbinitTask.from_input(inp, manager.workdir, runmode)

    manager.set_work_and_run(task)

    if manager.retcode != 0:
        return manager.retcode

    # Remove all files except those matching these regular expression.
    task.rmtree(exclude_wildcard="*.abi|*_WFK*|*_GSR.nc|*DEN-etsf.nc")

    task.rename("out_DS1_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
    task.rename("out_DS1_DEN-etsf.nc", "si_DEN-etsf.nc")
    task.rename("out_DS1_GSR.nc", "si_scf_GSR.nc")

    task.rename("out_DS2_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")
    task.rename("out_DS2_GSR.nc", "si_nscf_GSR.nc")

    task.remove_files("out_DS2_DEN-etsf.nc")

    manager.finalize()
    return manager.retcode 

if __name__ == "__main__":
    sys.exit(main())
