from __future__ import division, print_function, unicode_literals

import numpy as np
import abipy.data as data
import abipy.core

from pprint import pprint
from abipy.core.testing import *
from abipy.electrons.gsr import GSR_Reader, GSR_File

class GSRReaderTestCase(AbipyTest):

    def test_read_Si2(self):
        """Test the reading of GSR file."""
        path = data.ref_file("si_scf_GSR.nc")

        ref_dims = {
            "number_of_spins": 1
        }

        ref_int_values = {
            "space_group": 227,
        }

        ref_float_values = {
            "etotal": -8.8652767680604807,
        #    "primitive_vectors": np.reshape([0, 5.125, 5.125, 5.125, 0, 5.125,
        #                                     5.125, 5.125, 0], (3,3)),
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
            with self.assertRaises(GSR_Reader.Error): r.read_value("foobar")
            with self.assertRaises(GSR_Reader.Error): r.read_dimvalue("foobar")

            r.print_tree()
            for group in r.walk_tree():
                print("group: " + str(group))

            # Initialize pymatgen structure from GSR.
            structure = r.read_structure()
            self.assertTrue(isinstance(structure, abipy.core.Structure))


class GSRFileTestCase(AbipyTest):

    def test_methods(self):
        """GSRFile methods"""
        almost_equal = self.assertAlmostEqual

        with GSR_File(data.ref_file("si_scf_GSR.nc")) as gsr:
            print(repr(gsr))
            print(gsr)
            assert gsr.filepath == data.ref_file("si_scf_GSR.nc")
            assert gsr.nsppol == 1 and gsr.nspden == 1 and gsr.nspinor == 1
            assert gsr.mband == 8 and gsr.nband == 8 and gsr.nelect == 8 and len(gsr.kpoints) == 29
            almost_equal(gsr.energy.to("Ha"), -8.86527676798556)

            # Forces and stress
            self.assert_almost_equal(gsr.cartesian_forces.flat, 
                [-5.98330096024095e-30, -5.64111024387213e-30, 1.49693284867669e-29,
                  5.98330096024095e-30, 5.64111024387213e-30, -1.49693284867669e-29])

            #self.assert_almost_equal(gsr.stress_tensor.flat, 
            #    [0.000177139315441305, 0.0001771393154385, 0.000177139315437277, 
            #     0, 0, 5.34581779880965e-15])
            #self.assert_almost_equal(gsr.pressure, 0)

            # Test as_dict
            pprint(gsr.as_dict())
            #import json
            #with open("hello.json", "w") as fp:
            #    json.dump(gsr.as_dict(), fp)
            #assert 0

            # Test pymatgen computed_entries
            for inc_structure in (True, False):
                e = gsr.get_computed_entry(inc_structure=False)
                print(e)
                print(e.as_dict())
                assert gsr.energy == e.energy


if __name__ == "__main__":
    import unittest
    unittest.main()
