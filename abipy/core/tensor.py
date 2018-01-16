# coding: utf-8
"""
This module contains classes representing tensors and helper functions to change lattice.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import itertools
import numpy as np

from monty.dev import deprecated
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from abipy.core import Structure


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

    # TODO Remove
    #@deprecated(message="abipy.core.Tensor is deprecated and will be replaced by pymatgen tensor in v0.4")
    def __init__(self, red_tensor, lattice, space="r"):
        """
        Args:
            red_tensor: array-like object with the 9 cartesian components of the tensor
            lattice: Lattice object defining the reference system
            space:
                "r" if the lattice is a real space lattice
                "g" if the lattice is a reciprocal space lattice
        """
        self._reduced_tensor = red_tensor
        self._lattice = lattice
        self.space = space

        if space == "g":
            self._is_real_space = False
        elif space == "r":
            self._is_real_space = True
        else:
            raise ValueError("space should be either 'g' or 'r'")

    def __eq__(self, other):
        if other is None:  return False
        return (np.allclose(self.reduced_tensor, other.reduced_tensor) and
                self.lattice == other.lattice and
                self.space == other.space)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return self.to_string()

    def __str__(self):
        return repr(self)

    def to_string(self, verbose=0, with_reduced=False):
        lines = []
        app = lines.append

        app("Tensor in %s space." % self.space)
        app("")
        app("Cartesian coordinates:")
        app(str(self.cartesian_tensor))

        if with_reduced:
            app("")
            app(str(self.lattice))
            app("Reduced coordinates:")
            app(str(self.reduced_tensor))

        return "\n".join(lines)

    @property
    def lattice(self):
        return self._lattice

    @property
    def reduced_tensor(self):
        return self._reduced_tensor

    @property
    def is_real_space(self):
        return self._is_real_space

    @property
    def cartesian_tensor(self):
        mat = self._lattice.matrix
        return np.dot(np.dot(np.transpose(mat), self._reduced_tensor), mat)

    @classmethod
    def from_cartesian_tensor(cls, cartesian_tensor, lattice, space="r"):
        red_tensor = from_cart_to_red(cartesian_tensor, lattice)
        return cls(red_tensor, lattice,space)

    def symmetrize(self, structure):
        tensor = self._reduced_tensor

        if self._is_real_space:
            real_lattice = self._lattice
        else:
            real_lattice = self._lattice.reciprocal_lattice

        # I guess this is the reason why tensor.symmetrize (omega) is so slow!
        real_finder = SpacegroupAnalyzer(structure)

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
    def from_directions(cls, qpoints, values, lattice, space):
        """
        Build a `SymmetricTensor` from the values computed along 6 directions.

        Args:
            qpoints: fractional coordinates of 6 independent q-directions
            values: values of (q^T E q)/(q^T q) along the 6 qpoints
            lattice: `Lattice` object defining the reference system
            space: "r" if the lattice is a real space lattice
                   "g" if the lattice is a reciprocal space lattice
        """
        assert len(qpoints) == 6 and len(values) == len(qpoints)

        mat = lattice.matrix
        metric = np.dot(np.transpose(mat),mat)

        coeffs_red = np.zeros((6,6))

        for (iqpt,qpt) in enumerate(qpoints):
            metqpt = np.dot(metric,qpt)

            coeffs_red[iqpt,:] = [metqpt[0]**2,metqpt[1]**2,metqpt[2]**2,
                                  2*metqpt[0]*metqpt[1],2*metqpt[0]*metqpt[2],2*metqpt[1]*metqpt[2]]

            normqpt_red = np.dot(np.transpose(qpt),np.dot(metric,qpt))

            coeffs_red[iqpt,:] = coeffs_red[iqpt,:] / normqpt_red

        red_symm = np.linalg.solve(coeffs_red,values)

        red_tensor = [[red_symm[0],red_symm[3],red_symm[4]],
                      [red_symm[3],red_symm[1],red_symm[5]],
                      [red_symm[4],red_symm[5],red_symm[2]]]

        return cls(red_tensor, lattice, space)
