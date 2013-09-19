"""
TODO
"""
from __future__ import print_function, division

import itertools
import numpy as np

__all__ = [
    "Tensor",
    "SymmetricTensor",
]


class Tensor(object):
    """Representation of a 3x3 tensor"""
    def __init__(self, red_tensor, lattice):
        """
        Args:
            red_tensor:
                array-like object with the 9 cartesian components of the tensor
            lattice:
                Lattice object defining the reference system
        """
        self._reduced_tensor = red_tensor
        self._lattice = lattice

    @property
    def reduced_tensor(self):
        return self._reduced_tensor

    @property
    def cartesian_tensor(self):
        try:
            return self._cartesian_tensor
        except AttributeError:
            mat = self._lattice.matrix
            self._cartesian_tensor = np.dot(np.dot(np.transpose(mat),self._reduced_tensor),mat)
        
        return self._cartesian_tensor


    @classmethod
    def from_cartesian_tensor(cls, cartesian_tensor,lattice):
        mat = lattice.inv_matrix
        red_tensor = np.dot(np.dot(np.transpose(mat), cartesian_tensor), mat)
        return cls(red_tensor, lattice)

    def symmetrize(self, pointgroup, system="reciprocal_lattice"):
        tensor = self._reduced_tensor

        if system == "reciprocal_lattice":
            real_lattice = self._lattice.reciprocal_lattice
        elif system == "real_lattice":
            real_lattice = self._lattice
        else:
            raise AttributeError("system not understood")

        # TODO symmetrize

class SymmetricTensor(Tensor):
    """Representation of a 3x3 symmetric tensor"""
    @classmethod
    def from_directions(cls, qpoints, values, lattice):
        """
        Args:
            qpoints:
                fractional coordinates of 6 independent q-directions
            values:
                values of (q^T E q)/(q^T q) along the 6 qpoints
            lattice:
                Lattice object defining the reference system
        """

        assert (len(qpoints) == 6)
        assert (len(values) == len(qpoints))

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

