"""Tests for htc.FilesFile."""
from __future__ import print_function, division

import abipy.data as data

from abipy.core.testing import AbipyTest
from abipy.htc.input import *


class LdauLexxTest(AbipyTest):

    def test_nio(self):
        """Test LdauParams and LexxParams."""
        structure = data.structure_from_ucell("NiO")
        pseudos = data.pseudos("28ni.paw", "8o.2.paw")

        u = 8.0
        luj_params = LdauParams(usepawu=1, structure=structure)
        luj_params.luj_for_symbol("Ni", l=2, u=u, j=0.1*u, unit="eV")
        vars = luj_params.to_abivars()

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
