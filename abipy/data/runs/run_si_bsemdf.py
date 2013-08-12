#!/usr/bin/env python
from __future__ import division, print_function

import sys
import abipy.data as data  

from pymatgen.io.abinitio.abiobjects import AbiStructure
from pymatgen.io.abinitio.task import RunMode
from pymatgen.io.abinitio.pseudos import PseudoTable
from pymatgen.io.abinitio.calculations import bse_with_mdf
from abipy.data.runs import RunManager


def main():
    return 0
    structure = AbiStructure.asabistructure(data.cif_file("si.cif"))

    pseudos = PseudoTable(data.pseudos("14si.pspnc"))
    runmode = RunMode.sequential()

    manager = RunManager()

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


    work = bse_with_mdf(manager.workdir, runmode, structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
                       ecuteps, bs_loband, soenergy, mdf_epsinf,
                       accuracy="normal", spin_mode="unpolarized", smearing=None,
                       charge=0.0, scf_solver=None, **extra_abivars)

    manager.set_workflow_and_run(work)

    if manager.retcode !=0:
        return manager.retcode

    # Remove all files except those matching these regular expression.
    #work[0].rename("out_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
    #work[0].rename("out_DEN-etsf.nc", "si_DEN-etsf.nc")
    #work[0].rename("out_GSR.nc", "si_scf_GSR.nc")

    #work[1].rename("out_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")
    #work[1].rename("out_GSR.nc", "si_nscf_GSR.nc")

    #work[3].rename("out_SIGRES.nc", "si_g0w0ppm_SIGRES.nc")

    #work.rmtree(exclude_wildcard="*.abin|*.about|*_SIGRES.nc")

    manager.finalize()

    return manager.retcode 

if __name__ == "__main__":
    sys.exit(main())
