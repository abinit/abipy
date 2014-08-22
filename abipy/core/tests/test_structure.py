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

    def test_utils(self):
        """Test utilities for the generation of Abinit inputs."""
        structure = data.structure_from_ucell("MgB2")

        self.serialize_with_pickle(structure)

        pseudos = data.pseudos("12mg.pspnc", "5b.pspnc")
        nval = structure.num_valence_electrons(pseudos)
        self.assertEqual(nval, 8)
        shiftk = structure.calc_shiftk()
        self.assert_equal(shiftk, [[0.0, 0.0, 0.5]])


if __name__ == "__main__":
    import unittest
    unittest.main()
