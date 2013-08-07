#!/usr/bin/env python
from __future__ import division, print_function

import sys
import abipy.data as data  

from pymatgen.io.abinitio.abiobjects import AbiStructure
from pymatgen.io.abinitio.task import RunMode
from pymatgen.io.abinitio.pseudos import PseudoTable
from pymatgen.io.abinitio.launcher import SimpleResourceManager
from pymatgen.io.abinitio.calculations import g0w0_with_ppmodel


def main():
    structure = AbiStructure.asabistructure(data.cif_file("si.cif"))

    pseudos = PseudoTable(data.pseudos("14si.pspnc"))
    runmode = RunMode.sequential()

    kppa = scf_kppa = 1
    nscf_nband = 6
    nscf_ngkpt = [4,4,4]
    nscf_shiftk = [0.1, 0.2, 0.3]
    ecuteps, ecutsigx = 2, 2
    max_ncpus = 1

    extra_abivars = dict(
        ecut=12, 
        istwfk="*1",
    )

    g0w0 = g0w0_with_ppmodel("test_g0w0", runmode, structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                             accuracy="normal", spin_mode="unpolarized", smearing=None, 
                             ppmodel="godby", charge=0.0, scf_solver=None, inclvkb=0, sigma_nband=None, scr_nband=None,
                             **extra_abivars)

    retcodes = SimpleResourceManager(g0w0, max_ncpus).run()
    return max(retcodes)

    #bse = bse_with_mdf("hello_bse", runmode, structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
    #                   ecuteps, bs_loband, soenergy, mdf_epsinf,
    #                   accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
    #                   charge=0.0, scf_solver=None)


if __name__ == "__main__":
    sys.exit(main())
