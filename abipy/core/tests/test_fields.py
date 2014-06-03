#!/usr/bin/env python
"""Tests for core.density module"""
from __future__ import print_function, division

import numpy as np
import abipy.data as data 

from abipy.core.fields import *
from abipy.core.testing import *


class TestScalarField(AbipyTest):
    """Unit tests for ScalarField."""

    def test_base(self):
        """Testing ScalarField."""
        atrue = self.assertTrue
        aequal = self.assertEqual

        structure = data.structure_from_ucell("Si")
        nspinor, nsppol, nspden = 1,1,1
        xyz_shape = (2,3,4)
        datar = np.zeros((nsppol,) + xyz_shape)

        field = ScalarField(nspinor, nsppol, nspden, datar, structure, iorder="c")

        print(field)
        atrue(field.is_collinear)

        #aequal(field.datar.ndim, 2)
        #aequal(field.datar_xyz.ndim, 4)
        #aequal(field.datar_xyz.shape[-3:], xyz_shape)


if __name__ == "__main__":
    import unittest
    unittest.main()
