"""Tests for flowtk.abiobjects module."""
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.flowtk.abiobjects import LdauParams, LexxParams


class LdauLexxTest(AbipyTest):

    def test_nio(self):
        """Test LdauParams and LexxParams."""
        aequal, atrue = self.assertEqual, self.assertTrue

        structure = abidata.structure_from_ucell("NiO")
        pseudos = abidata.pseudos("28ni.paw", "8o.2.paw")

        u = 8.0
        luj_params = LdauParams(usepawu=1, structure=structure)
        luj_params.luj_for_symbol("Ni", l=2, u=u, j=0.1*u, unit="eV")
        avars = luj_params.to_abivars()

        self.serialize_with_pickle(luj_params, test_eq=False)

        atrue(avars["usepawu"] == 1)
        aequal(avars["lpawu"], "2 -1"),
        aequal(avars["upawu"], "8.0 0.0 eV")
        aequal(avars["jpawu"], "0.8 0.0 eV")

        # Cannot add UJ for non-existent species.
        with self.assertRaises(ValueError):
            luj_params.luj_for_symbol("Foo", l=2, u=u, j=0.1*u, unit="eV")

        # Cannot overwrite UJ.
        with self.assertRaises(ValueError):
            luj_params.luj_for_symbol("Ni", l=1, u=u, j=0.1*u, unit="eV")

        lexx_params = LexxParams(structure)
        lexx_params.lexx_for_symbol("Ni", l=2)
        avars = lexx_params.to_abivars()

        self.serialize_with_pickle(lexx_params, test_eq=False)

        aequal(avars["useexexch"], 1),
        aequal(avars["lexexch"], "2 -1")

        # Cannot add LEXX for non-existent species.
        with self.assertRaises(ValueError):
            lexx_params.lexx_for_symbol("Foo", l=2)

        # Cannot overwrite LEXX.
        with self.assertRaises(ValueError):
            lexx_params.lexx_for_symbol("Ni", l=1)
