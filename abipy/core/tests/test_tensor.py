"""Tests for tensor"""
from __future__ import print_function, division

import numpy as np
import abipy.data as data

from pymatgen.core.lattice import Lattice
from abipy.core.tensor import *
from abipy.core.testing import *

class TestTensor(AbipyTest):
    """Unit tests for Tensor."""

    def test_tensor(self):
        """Initialize Tensor"""

        lattice = Lattice.monoclinic(4,5,6,74)

        cartesian_tensor = [[2,3,0],[3,4,0],[0,0,6]]
        tensor = Tensor.from_cartesian_tensor(cartesian_tensor,lattice.reciprocal_lattice)
        red_tensor = tensor.reduced_tensor
        tensor2 = Tensor(red_tensor,lattice.reciprocal_lattice)
        assert(((np.abs(tensor2.cartesian_tensor)-np.abs(cartesian_tensor)) < 1E-8).all())

if __name__ == "__main__":
    import unittest
    unittest.main()
