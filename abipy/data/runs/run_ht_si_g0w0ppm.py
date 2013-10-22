#!/usr/bin/env python
from __future__ import division, print_function

import os
import abipy.data as data  

from pymatgen.io.abinitio.abiobjects import AbiStructure
from pymatgen.io.abinitio.calculations import g0w0_with_ppmodel
from abipy.data.runs import Tester, enable_logging


@enable_logging
def main():
    structure = AbiStructure.asabistructure(data.cif_file("si.cif"))

    pseudos = data.pseudos("14si.pspnc")

    tester = Tester()
    manager = tester.make_manager()

    scf_kppa = 40
    nscf_nband = 100
    #nscf_ngkpt = [4,4,4]
    #nscf_shiftk = [0.0, 0.0, 0.0]
    ecuteps, ecutsigx = 6, 8
    #scr_nband = 50
    #sigma_nband = 50

    extra_abivars = dict(
        ecut=8, 
        istwfk="*1",
        timopt=-1,
    )

    work = g0w0_with_ppmodel(structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                             accuracy="normal", spin_mode="unpolarized", smearing=None, 
                             ppmodel="godby", charge=0.0, inclvkb=2, sigma_nband=None, scr_nband=None,
                             workdir=tester.workdir, manager=manager, **extra_abivars)

    tester.set_work_and_run(work)

    if tester.retcode != 0:
        return tester.retcode

    # Remove all files except those matching these regular expression.
    work[3].rename("out_SIGRES.nc", "si_g0w0ppm_SIGRES.nc")

    work.rmtree(exclude_wildcard="*.abin|*.about|*_SIGRES.nc")

    tester.finalize()

    return tester.retcode


if __name__ == "__main__":
    import sys
    sys.exit(main())
