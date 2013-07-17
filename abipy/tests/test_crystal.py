#!/usr/bin/env python
"""Tests for core.crystal module."""
from __future__ import print_function, division

from abipy.core import Structure
from abipy.tests import ALL_NCFILES, AbipyTest


class TestStructure(AbipyTest):
    "Unit tests for Structure."

    def test_structure_from_ncfiles(self):
        "Init Structure objects from Netcdf data files"
        ncfiles = ALL_NCFILES
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
