"""Tests for tensor"""
from __future__ import print_function, division

import numpy as np
import abipy.data as data

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.finder import SymmetryFinder
from abipy.core.tensor import *
from abipy.core.testing import *

class TestTensor(AbipyTest):
    """Unit tests for Tensor."""

    def test_tensor(self):
        """Initialize Tensor"""

        lattice = Lattice.hexagonal(4,6)
        #rprimd = np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
        #rprimd = rprimd*10
        #lattice = Lattice(rprimd)
        structure = Structure(lattice, ["Ga", "As"],
                                      [[0, 0, 0], [0.5, 0.5, 0.5]])

        finder = SymmetryFinder(structure)

        spacegroup = finder.get_spacegroup()
        pointgroup = finder.get_point_group()

        cartesian_tensor = [[2,3,1.2],[3,4,1.0],[1.2,1.0,6]]

        tensor = Tensor.from_cartesian_tensor(cartesian_tensor,lattice.reciprocal_lattice,space="g")
        red_tensor = tensor.reduced_tensor
        tensor2 = Tensor(red_tensor,lattice.reciprocal_lattice,space="g")
        assert(((np.abs(tensor2.cartesian_tensor)-np.abs(cartesian_tensor)) < 1E-8).all())

        #print("non-symmetrized cartesian_tensor = ",tensor2.cartesian_tensor)
        tensor2.symmetrize()

        #print("symmetrized_cartesian_tensor = ",tensor2.cartesian_tensor)


if __name__ == "__main__":
    import unittest
    unittest.main()
