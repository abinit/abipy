from __future__ import division, print_function

import numpy as np
import abipy.data as data
import abipy.core

from abipy.core.testing import *
from abipy.electrons.gsr import GSR_Reader

class GSR_Reader_TestCase(AbipyTest):

    def setUp(self):
        formulas = ["Si2",]
        self.GSR_paths = d = {}
        for formula in formulas:
            d[formula] = data.ref_file(formula + "_GSR.nc")

    def test_read_Si2(self):
        """Test the reading of GSR file."""
        path = self.GSR_paths["Si2"]

        ref_dims = {
            "number_of_spins": 1
        }

        ref_int_values = {
            "space_group": 227,
            "number_of_states": np.reshape([15, 15], (1,2)),
        }

        ref_float_values = {
            "etotal": -8.85911566912484,
            "primitive_vectors": np.reshape([0, 5.125, 5.125, 5.125, 0, 5.125,
                                             5.125, 5.125, 0], (3,3)),
        }

        with GSR_Reader(path) as r:

            self.assertEqual(r.ngroups, 1)

            print(r.read_varnames())

            # Test dimensions.
            for (dimname, int_ref) in ref_dims.items():
                value = r.read_dimvalue(dimname)
                self.assert_equal(value, int_ref)

            # Test int variables
            for (varname, int_ref) in ref_int_values.items():
                value = r.read_value(varname)
                self.assert_equal(value, int_ref)

            # Test float variables
            for (varname, float_ref) in ref_float_values.items():
                value = r.read_value(varname)
                self.assert_almost_equal(value, float_ref)

            # Reading non-existent variables or dims should raise
            # a subclass of NetcdReder.
            with self.assertRaises(GSR_Reader.Error):
                r.read_value("foobar")

            with self.assertRaises(GSR_Reader.Error):
                r.read_dimvalue("foobar")

            r.print_tree()
            for group in r.walk_tree():
                print("group: " + str(group))

            # Initialize pymatgen structure from GSR.
            structure = r.read_structure()
            self.assertTrue(isinstance(structure, abipy.core.Structure))


if __name__ == "__main__":
    import unittest
    unittest.main()
