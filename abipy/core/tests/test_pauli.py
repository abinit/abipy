#!/usr/bin/env python
"""Tests for core.field module"""
import numpy as np

from abipy.core.testing import AbipyTest
from abipy.core.pauli import Pauli

class TestPauli(AbipyTest):
    """Unit tests for _Field."""

    def test_pauli(self):
        """Testing Pauli matrices."""

        pauli = Pauli()

        # Operate on a single matrix.
        matrix = np.array([[1,  3+2j],
                           [3-2j, 4]])

        cs = pauli.project_mats(matrix)
        assert len(cs) == 4 and cs.ndim == 1 and cs.shape == (4,)
        same_matrix = cs[0] * pauli.sigma_0 + cs[1] * pauli.sigma_x + cs[2] * pauli.sigma_y + cs[3] * pauli.sigma_z
        self.assert_equal(matrix, same_matrix)

        same_matrix = pauli.mats_from_projections(cs)
        self.assert_equal(matrix, same_matrix)

        # Array of matrices
        matrices = np.array([
            [[1, 2], [2, 4]],
            [[3, 4-2j], [1j, 4+2j]],
        ])

        ## Array of matrices projection
        cs_mat = pauli.project_mats(matrices)
        assert cs_mat.shape == (2, 4)

        for matrix, cs in zip(matrices, cs_mat):
            same_matrix = cs[0] * pauli.sigma_0 + cs[1] * pauli.sigma_x + cs[2] * pauli.sigma_y + cs[3] * pauli.sigma_z
            self.assert_equal(matrix, same_matrix)

        print(cs_mat)
        same_matrices = pauli.mats_from_projections(cs_mat)
        self.assert_equal(matrices, same_matrices)
