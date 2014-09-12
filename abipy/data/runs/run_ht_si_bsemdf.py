#!/usr/bin/env python
"""Calculation of the BSE spectrum with the HT interface."""
from __future__ import division, print_function

import sys
import os
from abipy import abilab
import abipy.data as abidata  

from pymatgen.io.abinitio.abiobjects import AbiStructure
from pymatgen.io.abinitio.calculations import bse_with_mdf


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    workdir = options.workdir
    if not options.workdir:
        workdir = os.path.basename(__file__).replace(".py", "").replace("run_","flow_") 

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else \
              abilab.TaskManager.from_file(options.manager)

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
    max_ncpus = 1
    ecuteps = 2

    extra_abivars = dict(
        ecut=12, 
        istwfk="*1",
    )

    flow = abilab.AbinitFlow(workdir=workdir, manager=manager)

    # BSE calculation with model dielectric function.
    work = bse_with_mdf(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
                       ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf,
                       accuracy="normal", spin_mode="unpolarized", smearing=None,
                       charge=0.0, scf_solver=None, **extra_abivars)

    flow.register_work(work)
    return flow.allocate()


@abilab.flow_main
def main(options):
    flow = build_flow(options)
    return flow.build_and_pickle_dump()


if __name__ == "__main__":
    sys.exit(main())
