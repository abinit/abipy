"""Tests for structure module"""
import numpy as np
import abipy.data as data

from abipy.core.structure import *
from abipy.core.testing import *

class TestStructure(AbipyTest):
    """Unit tests for Structure."""

    def test_structure_from_ncfiles(self):
        """Initialize Structure from Netcdf data files"""
        ncfiles = data.ALL_NCFILES
        assert ncfiles
        for filename in ncfiles:
            print("Reading file %s" % filename)
            structure = Structure.from_file(filename)
            print(structure)

            #structure.visualize(None)
            #xcart = structure.get_xcart()
            #normA = structure.calc_norm(structure.xred, space="r")
            #normB = np.sqrt( [np.dot(x,x) for x in xcart] )
            #self.assert_almost_equal(normA, normB)


if __name__ == "__main__":
    import unittest
    unittest.main()
