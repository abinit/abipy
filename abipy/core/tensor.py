"""
This module contains classes representing tensors and helper functions to change lattice.
"""
from __future__ import print_function, division

import itertools
import numpy as np

from abipy.core import Structure
from pymatgen.symmetry.finder import SymmetryFinder

__all__ = [
    "Tensor",
    "SymmetricTensor",
]


def from_cart_to_red(cartesian_tensor,lattice):
    mat = lattice.inv_matrix
    red_tensor = np.dot(np.dot(np.transpose(mat), cartesian_tensor), mat)
    return red_tensor

class Tensor(object):
    """Representation of a 3x3 tensor"""
    def __init__(self, red_tensor, lattice, space="r"):
        """
        Args:
            red_tensor:
                array-like object with the 9 cartesian components of the tensor
            lattice:
                Lattice object defining the reference system
        """
        self._reduced_tensor = red_tensor
        self._lattice = lattice
      
        if space == "g":
            self._is_real_space = False
        elif space == "r":
            self._is_real_space = True
        else:
            raise ValueError("space should be either 'g' or 'r'")

    @property
    def reduced_tensor(self):
        return self._reduced_tensor

    @property
    def is_real_space(self):
        return self._is_real_space

    @property
    def cartesian_tensor(self):
        mat = self._lattice.matrix
        return np.dot(np.dot(np.transpose(mat),self._reduced_tensor),mat)

    @classmethod
    def from_cartesian_tensor(cls, cartesian_tensor,lattice,space="r"):
        red_tensor = from_cart_to_red(cartesian_tensor,lattice)
        return cls(red_tensor, lattice,space)

    def symmetrize(self):
        tensor = self._reduced_tensor

        if self._is_real_space:
            real_lattice = self._lattice
        else:
            real_lattice = self._lattice.reciprocal_lattice

        # Create a fake structure to get point group
        real_structure = Structure(real_lattice,["H"],[[0,0,0]])
        real_finder = SymmetryFinder(real_structure)

        real_symmops = real_finder.get_point_group_operations(cartesian=True)

        cartesian_tensor = self.cartesian_tensor

        sym_tensor = np.zeros((3,3))

        my_tensor = cartesian_tensor

        for real_sym in real_symmops:
             mat = real_sym.rotation_matrix
             prod_sym = np.dot(np.transpose(mat),np.dot(cartesian_tensor,mat))
             sym_tensor = sym_tensor + prod_sym

        sym_tensor = sym_tensor/len(real_symmops)

        self._reduced_tensor = from_cart_to_red(sym_tensor,self._lattice)

class SymmetricTensor(Tensor):
    """Representation of a 3x3 symmetric tensor"""
    @classmethod
    def from_directions(cls, qpoints, values, lattice):
        """
        Build a `SymmetricTensor` from the values computed along 6 directions. 

        Args:
            qpoints:
                fractional coordinates of 6 independent q-directions
            values:
                values of (q^T E q)/(q^T q) along the 6 qpoints
            lattice:
                `Lattice` object defining the reference system
        """

        assert len(qpoints) == 6
        assert len(values) == len(qpoints)

        mat = lattice.matrix
        metric = np.dot(np.transpose(mat),mat)

        coeffs_red = np.zeros((6,6))

        for (iqpt,qpt) in enumerate(qpoints):
            metqpt = np.dot(metric,qpt)

            coeffs_red[iqpt,:] = [metqpt[0]**2,metqpt[1]**2,metqpt[2]**2,
                                  2*metqpt[0]*metqpt[1],2*metqpt[0]*metqpt[2],2*metqpt[1]*metqpt[2]]

            normqpt_red = np.dot(np.transpose(qpt),np.dot(metric,qpt))

            coeffs_red[iqpt,:] = coeffs_red[iqpt,:] / (normqpt_red)


        red_symm = np.linalg.solve(coeffs_red,values)


        red_tensor = [[red_symm[0],red_symm[3],red_symm[4]],
                      [red_symm[3],red_symm[1],red_symm[5]],
                      [red_symm[4],red_symm[5],red_symm[2]]]

        return cls(red_tensor,lattice)

