# coding: utf-8
import sys
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.abio.decorators as ideco

from abipy.core.testing import AbipyTest
from abipy.abio.factories import *


class DecoratorTest(AbipyTest):

    def setUp(self):
        # Si ebands
        si_structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        self.si_ebands = ebands_input(si_structure, abidata.pseudos("14si.pspnc"), ecut=2, kppa=10)

        # Reference input string. Used to test if decorators do not change the initial Input.
        self.si_ebands_inpstr = str(self.si_ebands)

        # NiO bands with PAW
        nio_structure = abidata.structure_from_ucell("NiO")
        self.nio_ebands = ebands_input(nio_structure, abidata.pseudos("28ni.paw", "8o.2.paw"),
                                       ecut=2, pawecutdg=4, kppa=10)

        self.nio_ebands_inpstr = str(self.nio_ebands)

    def tearDown(self):
        """Testing if initial inputs are unchanged."""
        assert all(not inp.decorators for inp in self.si_ebands)
        assert self.si_ebands_inpstr == str(self.si_ebands)

        assert all(not inp.decorators for inp in self.nio_ebands)
        assert self.nio_ebands_inpstr == str(self.nio_ebands)

    def validate_inp(self, inp, ndec=1):
        # Hack needed because ecut is not in the pseudos.
        inp.set_vars(ecut=3)

        #v = inp.validate()
        #if v.retcode != 0:
        #    raise RuntimeError(v.err)
        #else:
        #    print("Valid input!")

        # Test validity of individual datasets.
        for dtset in inp.split_datasets():
            v = dtset.abivalidate()
            #assert dtset.decorators == inp.decorators
            #assert len(dtset.decorators) == ndec

            if v.retcode != 0:
                raise RuntimeError("Wrong input. See {0}".format(v))
            else:
                print("Valid input!")

    def test_spin_decorator(self):
        """Testing spin decorator."""
        spinor_deco = ideco.SpinDecorator("spinor")
        self.assertMSONable(spinor_deco)
        print(spinor_deco)

        new_inp = spinor_deco(self.si_ebands)
        print(new_inp)

        # kptopt is set to 4 if non-collinear magnetism and kptopt == 3 is not specified.
        for dt in new_inp:
            assert dt["nsppol"] == 1 and dt["nspinor"] == 2 and dt["kptopt"] == 4
        #self.validate_inp(new_inp)

        # kptopt should not be changes if it's set to 3 and non-collinear magnetism
        inp_with_kpt3 = self.si_ebands.deepcopy()
        inp_with_kpt3.kptopt = 3

        # FIXME: Here there's a bug because get should check the global variables!
        #for dt in spinor_deco(inp_with_kpt3):
        #    assert dt["nsppol"] == 1 and dt["nspinor"] == 2 and dt["kptopt"] == 3

    def test_smearing_decorator(self):
        """Testing electronic smearing decorator."""
        smearing_deco = ideco.SmearingDecorator("fermi_dirac:0.1 eV")
        self.assertMSONable(smearing_deco)

        new_inp = smearing_deco(self.si_ebands)
        self.validate_inp(new_inp)

    def test_xcdecorator(self):
        """Testing XCdecorator."""
        xc_deco = ideco.XcDecorator(17)
        self.assertMSONable(xc_deco)

        new_inp = xc_deco(self.si_ebands)
        self.validate_inp(new_inp)

    def test_ldau_decorators(self):
        """Testing LdaUDecorator."""
        symbols_luj = dict(Ni=dict(l=2, u=5.0, j=0.5))

        ldau_deco = ideco.LdaUDecorator(symbols_luj, usepawu=1, unit="eV")
        self.assertMSONable(ldau_deco)

        new_inp = ldau_deco(self.nio_ebands)
        new_inp.set_vars(chkprim=0, ecut=3, pawecutdg=3)
        print(new_inp)
        self.validate_inp(new_inp)
        #assert 0

        # LDA+U only if PAW
        with self.assertRaises(ldau_deco.Error):
            ldau_deco(self.si_ebands)

    def test_lexx_decorators(self):
        """Testing LexxDecorator."""
        lexx_deco = ideco.LexxDecorator({"Ni": 2})
        self.assertMSONable(lexx_deco)

        new_inp = lexx_deco(self.nio_ebands)
        new_inp.set_vars(chkprim=0, ecut=3, pawecutdg=3)
        print(new_inp)
        self.validate_inp(new_inp)
        #assert 0

    def test_new_with_decorators(self):
        """Testing AbinitInput.new_with_decorators."""
        spinor_deco = ideco.SpinDecorator("spinor")
        smearing_deco = ideco.SmearingDecorator("nosmearing")
        new_inp = self.si_ebands.new_with_decorators(spinor_deco)
        new_inp = self.si_ebands.new_with_decorators([spinor_deco, smearing_deco])


if __name__ == '__main__':
    import unittest
    unittest.main()
