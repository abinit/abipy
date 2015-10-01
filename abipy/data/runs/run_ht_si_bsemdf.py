#!/usr/bin/env python
"""Calculation of the BSE spectrum with the HT interface."""
from __future__ import division, print_function, unicode_literals

import sys
import os
import abipy.data as abidata  

from abipy import abilab


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_") 

    # Initialize pseudos and Structure.
    pseudos = abidata.pseudos("14si.pspnc")
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    kppa = scf_kppa = 1
    nscf_nband = 6
    nscf_ngkpt = [4,4,4]
    nscf_shiftk = [0.1, 0.2, 0.3]
    bs_loband = 2
    bs_nband = nscf_nband
    soenergy = 0.7
    mdf_epsinf = 12
    ecuteps = 2
    ecut = 12

    #extra_abivars = dict(
    #    ecut=12, 
    #    istwfk="*1",
    #)

    flow = abilab.Flow(workdir=workdir, manager=options.manager)

    # BSE calculation with model dielectric function.
    multi = abilab.bse_with_mdf_inputs(
        structure, pseudos, 
        scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk, 
        ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf, 
        ecut=ecut,#$ pawecutdg=None, 
        exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="unpolarized", 
        smearing=None)
        #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

    work = abilab.BseMdfWork(scf_input=multi[0], nscf_input=multi[1], bse_inputs=multi[2:])

    #from pymatgen.io.abinit.calculations import bse_with_mdf_work
    #work = bse_with_mdf_work(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
    #                         ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf,
    #                         accuracy="normal", spin_mode="unpolarized", smearing=None,
    #                         charge=0.0, scf_solver=None, **extra_abivars)

    flow.register_work(work)
    return flow


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    flow.build_and_pickle_dump()
    return flow


if __name__ == "__main__":
    sys.exit(main())
