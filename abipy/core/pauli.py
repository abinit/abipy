# coding: utf-8
"""Pauli matrices and operations associated to them."""
from __future__ import annotations

import numpy as np


class Pauli:
    """
    Pauli matrices
    """
    def __init__(self):
        self.sigma_0 = np.eye(2)                      # 2x2 identity matrix
        self.sigma_x = np.array([[0, 1], [1, 0]])     # Pauli-X (σₓ)
        self.sigma_y = np.array([[0, -1j], [1j, 0]])  # Pauli-Y (σᵧ)
        self.sigma_z = np.array([[1, 0], [0, -1]])    # Pauli-Z (σ_z)

    def project_mats(self, mats: np.ndarray) -> np.ndarray:
        """
        Project a batch of 2x2 Hermitian matrices onto the Pauli basis.

        This method computes the coefficients of the given matrices in the Pauli basis
        (identity, sigma_x, sigma_y, sigma_z) for each 2x2 matrix in the input array.
        The input array must have shape [..., 2, 2], where the trailing dimensions represent
        individual 2x2 matrices.

        Parameters
        ----------
        mats: np.ndarray
            An array of shape [..., 2, 2], where each trailing 2x2 subarray is a matrix
            to be projected onto the Pauli basis.

        Returns
        -------
        np.ndarray
            An array of shape [..., 4], where the last dimension contains the coefficients
            (a_0, a_x, a_y, a_z) for the identity matrix, sigma_x, sigma_y, and sigma_z, respectively.

        Examples
        --------
        >>> import numpy as np
        >>> mats = np.array([[[1, 0], [0, 1]], [[0, 1], [1, 0]]])  # Shape (2, 2, 2)
        >>> result = project_mats(mats)
        >>> print(result)  # Shape (2, 4)
        [[ 0.5  0.   0.   0.5]
         [ 0.   0.5  0.   0. ]]
        """
        # Ensure input is at least 2D and has shape [..., 2, 2]
        mats = np.asarray(mats)
        if mats.shape[-2:] != (2, 2):
            raise ValueError(f"Input must have shape [..., 2, 2] while it is {mats.shape[-2:]}")

        # Compute traces along the last two axes
        a_0 = 0.5 * np.trace(mats, axis1=-2, axis2=-1)
        a_x = 0.5 * np.trace(np.einsum('ij,...jk->...ik', self.sigma_x, mats), axis1=-2, axis2=-1)
        a_y = 0.5 * np.trace(np.einsum('ij,...jk->...ik', self.sigma_y, mats), axis1=-2, axis2=-1)
        a_z = 0.5 * np.trace(np.einsum('ij,...jk->...ik', self.sigma_z, mats), axis1=-2, axis2=-1)

        # Combines the coefficients  a_0, a_x, a_y, a_z into a single array with shape [..., 4].
        return np.stack((a_0, a_x, a_y, a_z), axis=-1)

    def mats_from_projections(self, coefficients: np.ndarray) -> np.ndarray:
        """
        Reconstruct 2x2 matrices from their Pauli basis coefficients.

        This method takes an array of coefficients (a_0, a_x, a_y, a_z) and reconstructs
        the corresponding 2x2 matrices using the Pauli basis (identity, sigma_x, sigma_y, sigma_z).

        Parameters
        ----------
        coefficients : np.ndarray
            An array of shape [..., 4], where the last dimension contains the coefficients
            (a_0, a_x, a_y, a_z) for the identity matrix, sigma_x, sigma_y, and sigma_z, respectively.

        Returns
        -------
        np.ndarray
            An array of shape [..., 2, 2], where each trailing 2x2 matrix is reconstructed
            from the input coefficients.

        Examples
        --------
        >>> coefficients = np.array([[0.5, 0, 0, 0.5], [0, 0.5, 0, 0]])  # Shape (2, 4)
        >>> result = reconstruct_mats(coefficients)
        >>> print(result)  # Shape (2, 2, 2)
        [[[1. 0.]
          [0. 1.]]
         [[0. 1.]
          [1. 0.]]]
        """
        # Ensure input is at least 1D and has shape [..., 4]
        coefficients = np.asarray(coefficients)
        if coefficients.shape[-1] != 4:
            raise ValueError(f"Input must have shape [..., 4] while it is {coefficients.shape[-1]}")

        # Extract the coefficients
        a_0, a_x, a_y, a_z = coefficients[..., 0], coefficients[..., 1], coefficients[..., 2], coefficients[..., 3]

        # Reconstruct the matrix using the Pauli basis
        reconstructed = (
            a_0[..., np.newaxis, np.newaxis] * self.sigma_0 +
            a_x[..., np.newaxis, np.newaxis] * self.sigma_x +
            a_y[..., np.newaxis, np.newaxis] * self.sigma_y +
            a_z[..., np.newaxis, np.newaxis] * self.sigma_z
        )

        return reconstructed

    #def project_mats_nspden(self, mats_nspden: np.ndarray) -> np.ndarray:
    #    if mats_nspden.shape[0] % 4 != 0:
    #        raise ValueError(f"The first dimension must be divisible by 4 while it is {mats_nspden.shape[0]}")

    #    # Reshape the array to [2, 2, ...]
    #    mats_22 = mats_nspden.reshape(2, 2, *arr.shape[1:])
    #    # Move [2, 2] axes to the end for consistency
    #    mats = np.moveaxis(mats_22, (0, 1), (-2, -1))
    #    coeffs = self.project_mats(mats)
    #    return coeffs
