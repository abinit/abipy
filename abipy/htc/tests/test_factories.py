# coding: utf-8
from __future__ import unicode_literals, division, print_function

import sys
import abipy.data as abidata  
import abipy.abilab as abilab

from abipy.htc.factories import *
from abipy.core.testing import AbipyTest


class FactoryTest(AbipyTest):

    def validate_inp(self, inp):
        # Hack neede because ecut is not in the pseudos.
        inp.set_vars(ecut=3)

        v = inp.validate()
        if v.retcode != 0:
            raise RuntimeError(v.err)
        else:
            print("Valid input!")

        # Test validity of individual datasets.
        for dtset in inp.split_datasets():
            v = dtset.validate()

    def test_factories(self):
        pseudos = abidata.pseudos("14si.pspnc")
        structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        scf_kppa, scf_nband, nscf_nband, dos_kppa = 1000, 10, 10, 4

        inp = ebands_input(structure, pseudos, scf_kppa)
        print(inp)
        self.validate_inp(inp)
        sys.exit(0)

        inp = ebands_input(structure, pseudos, scf_kppa, nscf_nband, dos_kppa=dos_kppa)
        print(inp)
        self.validate_inp(inp)
        #sys.exit(0)

        inp = ion_ioncell_relax_input(structure, pseudos, scf_kppa, scf_nband)
                            #accuracy="normal", spin_mode="polarized",
                            #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

        self.validate_inp(inp)

        ecuteps, ecutsigx = 3, 2
        inp = g0w0_with_ppmodel_input(structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx)
                                      #accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
                                      #ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None,
                                      #sigma_nband=None, gw_qprange=1):

        self.validate_inp(inp)

        nscf_ngkpt, nscf_shiftk = [2,2,2], [[0,0,0]]
        inp = bse_with_mdf_input(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
                                 ecuteps=2, bs_loband=1, bs_nband=2, soenergy="0.1 eV", mdf_epsinf=12)
                                 #exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized", 
                                 #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

        print(inp)
        self.validate_inp(inp)
        
        inps = scf_phonons_inputs(structure, pseudos, scf_kppa)
                                 #ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
                                 #smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
        print(inps[0])
        print(inps[1])
        pprint(inps[1].get_irred_perts())

        for inp in inps:
            self.validate_inp(inp)


if __name__ == '__main__':
    import unittest
    unittest.main()
