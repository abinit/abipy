#!/usr/bin/env python
"""Calculation of the BSE spectrum with the HT interface."""
from __future__ import division, print_function

import sys
from abipy import abilab
import abipy.data as data  

from pymatgen.io.abinitio.abiobjects import AbiStructure
from pymatgen.io.abinitio.calculations import bse_with_mdf
from abipy.data.runs import Tester, enable_logging


@enable_logging
def main():
    #return 0
    pseudos = data.pseudos("14si.pspnc")
    structure = abilab.Structure.from_file(data.cif_file("si.cif"))

    tester = Tester()
    manager = abilab.TaskManager.from_user_config()

    kppa = scf_kppa = 1
    nscf_nband = 6
    nscf_ngkpt = [4,4,4]
    nscf_shiftk = [0.1, 0.2, 0.3]
    bs_loband = 2
    soenergy = 0.7
    mdf_epsinf = 12
    max_ncpus = 1
    ecuteps = 2

    extra_abivars = dict(
        ecut=12, 
        istwfk="*1",
    )

    flow = abilab.AbinitFlow(workdir=tester.workdir, manager=manager)

    # BSE calculation with model dielectric function.
    work = bse_with_mdf(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
                       ecuteps, bs_loband, soenergy, mdf_epsinf,
                       accuracy="normal", spin_mode="unpolarized", smearing=None,
                       charge=0.0, scf_solver=None, **extra_abivars)

    flow.register_work(work)
    flow.allocate()

    #tester.set_flow_and_run(flow)
    #tester.set_work_and_run(work)
    #if tester.retcode !=0:
    #    return tester.retcode

    # Remove all files except those matching these regular expression.
    #work[0].rename("out_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
    #work[0].rename("out_DEN-etsf.nc", "si_DEN-etsf.nc")
    #work[0].rename("out_GSR.nc", "si_scf_GSR.nc")

    #work[1].rename("out_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")
    #work[1].rename("out_GSR.nc", "si_nscf_GSR.nc")

    #work[3].rename("out_SIGRES.nc", "si_g0w0ppm_SIGRES.nc")

    #work.rmtree(exclude_wildcard="*.abin|*.about|*_SIGRES.nc")

    #tester.finalize()
    #return tester.retcode 

if __name__ == "__main__":
    sys.exit(main())
