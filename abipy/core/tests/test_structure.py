"""Tests for structure module"""
from __future__ import print_function, division

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
            print("About to read file %s" % filename)
            structure = Structure.from_file(filename)
            print(structure)

            self.assertTrue(structure.has_spacegroup)

            # Call pymatgen machinery to get the high-symmetry stars.
            print(structure.hsym_stars)

            structure.export(".xsf")


if __name__ == "__main__":
    import unittest
    unittest.main()
