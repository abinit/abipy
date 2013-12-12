#!/usr/bin/env python
"""Calculation of the BSE spectrum with the HT interface."""
from __future__ import division, print_function

import sys
import os
from abipy import abilab
import abipy.data as abidata  

from pymatgen.io.abinitio.abiobjects import AbiStructure
from pymatgen.io.abinitio.calculations import bse_with_mdf
from abipy.data.runs import AbipyTest, MixinTest

class HtBseMdfFlowTest(AbipyTest, MixinTest):
    """
    Unit test for the flow defined in this module.  
    Users who just want to learn how to use this flow can ignore this section.
    """
    def setUp(self):
        super(HtBseMdfFlowTest, self).setUp()
        self.init_dirs()
        self.flow = build_flow()

    # Remove all files except those matching these regular expression.
    #work[0].rename("out_WFK_0-etsf.nc", "si_scf_WFK-etsf.nc")
    #work[0].rename("out_DEN-etsf.nc", "si_DEN-etsf.nc")
    #work[0].rename("out_GSR.nc", "si_scf_GSR.nc")
                                                                       
    #work[1].rename("out_WFK_0-etsf.nc", "si_nscf_WFK-etsf.nc")
    #work[1].rename("out_GSR.nc", "si_nscf_GSR.nc")
                                                                       
    #work[3].rename("out_SIGRES.nc", "si_g0w0ppm_SIGRES.nc")
                                                                       
    #work.rmtree(exclude_wildcard="*.abin|*.about|*_SIGRES.nc")


def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed)
    workdir = os.path.basename(os.path.abspath(__file__).replace(".py", "")) if not options.workdir else options.workdir

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config() if not options.manager else options.manager

    pseudos = abidata.pseudos("14si.pspnc")
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    kppa = scf_kppa = 1
    nscf_nband = 6
    nscf_ngkpt = [4,4,4]
    nscf_shiftk = [0.1, 0.2, 0.3]
    bs_loband = 2
    bs_nban = nscf_nband
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
