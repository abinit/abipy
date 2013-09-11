#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  
import abipy.abilab as abilab

from abipy.data.runs import Tester

def main():

    # Change the value of ngkpt below to perform a GW calculation with a different ngkpt.
    inp = make_input(ngkpt=[2,2,2])
    print(inp)

    # Create the task defining the calculation and run it.
    tester = Tester()
    manager = tester.make_manager()

    task = abilab.AbinitTask.from_input(inp, tester.workdir, manager)

    tester.set_work_and_run(task)

    if tester.retcode != 0:
        return tester.retcode

    # Remove all files except those matching these regular expression.
    task.rmtree(exclude_wildcard="*.abi|*.abo|*SIGRES.nc")

    task.rename("out_DS4_SIGRES.nc", "si_g0w0ppm_nband10_SIGRES.nc")
    task.rename("out_DS5_SIGRES.nc", "si_g0w0ppm_nband20_SIGRES.nc")
    task.rename("out_DS6_SIGRES.nc", "si_g0w0ppm_nband30_SIGRES.nc")

    tester.finalize()
    return tester.retcode 

def make_input(ngkpt):

    # Crystalline silicon
    # Calculation of the GW correction to the direct band gap in Gamma
    # Dataset 1: ground state calculation 
    # Dataset 2: NSCF calculation 
    # Dataset 3: calculation of the screening 
    # Dataset 4: calculation of the Self-Energy matrix elements (GW corrections)

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
    inp.ecut = 6
    inp.timopt = -1
    inp.istwfk = "*1"
    inp.set_kmesh(dtset=0, **gw_kmesh)

    # Dataset 1 (GS run)
    inp.set_kmesh(dtset=1, **scf_kmesh)
    inp.tolvrs1 = 1e-6
    inp.nband1 = 4

    # Dataset 2 (NSCF run)
    # Here we select the second dataset directly with the syntax inp[2]
    inp[2].set_variables(iscf=-2,
                         getden=1,
                         tolwfr=1e-12,
                         nband=35,
                         nbdbuf=5,
                        )

    # Dataset3: Calculation of the screening.
    screening_run = dict(
        optdriver=3,   
        getkss=2,      
        nband=25,    
        ecutwfn=inp.ecut,   
        symchi=1,
        inclvkb=0,
        ecuteps=4.0,    
        ppmfrq="16.7 eV",
    )

    inp.set_variables(dtset=3, **screening_run)

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
        sigma_run = dict(
            optdriver=4,
            getkss=2,      
            getscr=3,     
            nband=nband,      
            ecutwfn=inp.ecut,
            ecuteps=4.0,
            ecutsigx=6.0,
            symsigma=1,
        )

        inp.set_variables(dtset=4+idx, **sigma_run)
        inp.set_kptgw(kptgw, bdgw, dtset=4+idx)

    return inp


if __name__ == "__main__":
    import sys
    sys.exit(main())
