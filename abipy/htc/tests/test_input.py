"""Tests for htc.FilesFile."""
from __future__ import print_function, division

import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.htc.input import *

class AbiInputTest(AbipyTest):

    def test_input(self):
        """Testing AbiInput."""
        aequal = self.assertEqual
        inp = AbiInput(pseudos=abidata.pseudos("28ni.paw", "8o.2.paw"), ndtset=2, comment="NiO calculation")
        inp.set_structure(abidata.structure_from_ucell("NiO"))
        print(inp)

        aequal(inp.ndtset, 2)
        aequal(inp.ispaw, True)

        # Set global variables.
        inp.set_variables(ecut=10)

        # Setting an unknown variable should raise an error.
        with self.assertRaises(inp.Error):
            inp.set_variables(foobar=10)


class LdauLexxTest(AbipyTest):

    def test_nio(self):
        """Test LdauParams and LexxParams."""
        structure = abidata.structure_from_ucell("NiO")
        pseudos = abidata.pseudos("28ni.paw", "8o.2.paw")

        u = 8.0
        luj_params = LdauParams(usepawu=1, structure=structure)
        luj_params.luj_for_symbol("Ni", l=2, u=u, j=0.1*u, unit="eV")
        vars = luj_params.to_abivars()

        self.serialize_with_pickle(luj_params, test_eq=False)

        self.assertTrue(vars["usepawu"] == 1),
        self.assertTrue(vars["lpawu"] ==  "2 -1"),
        self.assertTrue(vars["upawu"] == "8.0 0.0 eV"),
        self.assertTrue(vars["jpawu"] == "0.8 0.0 eV"),

        # Cannot add UJ for non-existent species.
        with self.assertRaises(ValueError):
            luj_params.luj_for_symbol("Foo", l=2, u=u, j=0.1*u, unit="eV")

        # Cannot overwrite UJ.
        with self.assertRaises(ValueError):
            luj_params.luj_for_symbol("Ni", l=1, u=u, j=0.1*u, unit="eV")

        lexx_params = LexxParams(structure)
        lexx_params.lexx_for_symbol("Ni", l=2)
        vars = lexx_params.to_abivars()

        self.serialize_with_pickle(lexx_params, test_eq=False)

        self.assertTrue(vars["useexexch"] == 1),
        self.assertTrue(vars["lexexch"] ==  "2 -1"),

        # Cannot add LEXX for non-existent species.
        with self.assertRaises(ValueError):
            lexx_params.lexx_for_symbol("Foo", l=2)
                                                                            
        # Cannot overwrite LEXX.
        with self.assertRaises(ValueError):
            lexx_params.lexx_for_symbol("Ni", l=1)


if __name__ == "__main__": 
    import unittest
    unittest.main()
