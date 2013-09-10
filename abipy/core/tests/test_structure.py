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

        for filename in data.WFK_NCFILES + data.GSR_NCFILES:
            print("About to read file %s" % filename)
            structure = Structure.from_file(filename)
            print(structure)

            # All nc files produced by ABINIT should have info on the spacegroup.
            self.assertTrue(structure.has_spacegroup)

            # Call pymatgen machinery to get the high-symmetry stars.
            print(structure.hsym_stars)

            if self.which("xcrysden") is not None:
                # Export data in Xcrysden format.
                structure.export(".xsf")


if __name__ == "__main__":
    import unittest
    unittest.main()
