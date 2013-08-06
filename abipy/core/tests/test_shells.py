""""Tests for the kpoints.shells module."""
from __future__ import print_function, division

from abipy.core.shells import *
from abipy.core.testing import *

##########################################################################################

class TestShells(AbipyTest):
    "Unit Tests for Shells"

    def test_base(self):
        """Basic tests for Shells."""
        values = [0, 1, 2, 3, -2, -4, -3, 4, 3, 5, -4]
        func = None
        #func = abs
        shells = Shells(values, func)

        print("len:", len(shells))
        print("shells.values", shells.values)

        sort_values = values[:]
        sort_values.sort()

        sh_values = list()
        for sh in shells:
            sh_values.extend(sh.items)

        self.assertEqual(sort_values, sh_values)

        shells.get_from_value(1.0)
