# coding: utf-8
"""
Functions to perform analytic continuation with Pade'
Some of these routines have been directly translated from the Fortran version
implemented in ABINIT.
"""
from __future__ import annotations

import numpy as np


class SigmaPade:
    """
    High-level interface to perform the analytic continuation of the self-energy with the Pade' method.
    """
    def __init__(self, zs, f_zs):
        """
        Args:
            zs: Complex z-values.
            f_zs: Values of f(zs).
        """
        self.zs, self.f_zs = zs, f_zs
        if len(zs) != len(f_zs):
            raise ValueError(f"{len(zs)=} != {len(f_zs)=}")

    def eval(self, z_evals: np.ndarray) -> tuple[np.ndarray, np.np.ndarray]:
        """
        Compute the Pade fit for a list of zz points.
        Return f(z_evals) and f'(z_evals)
        """
        nn = len(z_evals)
        sws = np.zeros(nn, dtype=complex)
        dsdws = np.zeros(nn, dtype=complex)

        for iz, z_eval in enumerate(z_evals):
            sw, dsw = self._eval_one(z_eval)
            sws[iz] = sw
            dsdws[iz] = dsw

        return sws, dsdws

    def _eval_one(self, z_eval) -> tuple:
        """
        Pade for a single point z_eval.
        """
        # if z_eval is in 2 or 3 quadrant, avoid the branch cut in the complex plane using Sigma(-iw) = Sigma(iw)*.
        # See also sigma_pade_eval in m_dyson_solver.F90
        if z_eval.real > 0.0:
            cval = pade(self.zs, self.f_zs, z_eval)
            dz_cval = dpade(self.zs, self.f_zs, z_eval)
        else:
            cval = pade(-self.zs, self.f_zs.conj(), z_eval)
            dz_cval = dpade(-self.zs, self.f_zs.conj(), z_eval)

        return cval, dz_cval


def pade(zs: np.ndarray, f_zs: np.ndarray, z_eval) -> complex:
    """
    Calculate the Pade approximant of the function f_zs at z_eval.

    Args
      zs: Input array of complex numbers.
      f_zs: Input array of complex numbers.
      z_eval: Point at which to evaluate the Pade approximant.

    Returns: The Pade approximant value at z_eval.
    """
    if (n := len(zs)) != len(f_zs):
        raise ValueError(f"{len(zs)=} != {len(f_zs)=}")

    if not np.iscomplexobj(zs):
        raise TypeError(f"zs should be complex, but got {zs.dtype=}")

    # Calculate Pade coefficients
    a = calculate_pade_a(zs, f_zs)

    # Initialize Az and Bz arrays
    Az = np.zeros(n + 1, dtype=complex)
    Bz = np.zeros(n + 1, dtype=complex)
    Az[0] = 0.0 + 0.0j  # czero
    Az[1] = a[0]
    Bz[0] = 1.0 + 0.0j  # cone
    Bz[1] = 1.0 + 0.0j  # cone

    # Calculate Az and Bz recursively
    for i in range(1, n):
        Az[i + 1] = Az[i] + (z_eval - zs[i - 1]) * a[i] * Az[i - 1]
        Bz[i + 1] = Bz[i] + (z_eval - zs[i - 1]) * a[i] * Bz[i - 1]

    # Compute the Pade approximant
    pade_value = Az[n] / Bz[n]

    # Optional: print debugging information
    # print("Pade approximant:", pade_value)
    # print("Bz(n):", Bz[n])
    # if np.isclose(Bz[n].real, 0.0) and np.isclose(Bz[n].imag, 0.0):
    #     print("Warning: Bz(n) is close to zero:", Bz[n])

    return pade_value


def dpade(zs: np.ndarray, f_zs: np.ndarray, z_eval: complex) -> complex:
    """
    Calculate the derivative of the Pade approximant of the function f_zs at z_eval.

    Args:
        zs: Input array of complex numbers.
        f_zs: Input array of complex numbers.
        z_eval: Point at which to evaluate the derivative of the Pade approximant.

    Returns: The derivative of the Pade approximant value at z_eval.
    """
    if (n := len(zs)) != len(f_zs):
        raise ValueError(f"{len(zs)=} != {len(f_zs)=}")

    if not np.iscomplexobj(zs):
        raise TypeError(f"zs should be complex, but got {zs.dtype=}")

    # Calculate Pade coefficients
    a = calculate_pade_a(zs, f_zs)

    # Initialize Az, Bz, dAz, dBz arrays
    Az = np.zeros(n + 1, dtype=complex)
    Bz = np.zeros(n + 1, dtype=complex)
    dAz = np.zeros(n + 1, dtype=complex)
    dBz = np.zeros(n + 1, dtype=complex)

    Az[0] = 0.0 + 0.0j  # czero
    Az[1] = a[0]
    Bz[0] = 1.0 + 0.0j  # cone
    Bz[1] = 1.0 + 0.0j  # cone
    dAz[0] = 0.0 + 0.0j  # czero
    dAz[1] = 0.0 + 0.0j  # czero
    dBz[0] = 0.0 + 0.0j  # czero
    dBz[1] = 0.0 + 0.0j  # czero

    # Calculate Az, Bz, dAz, and dBz recursively
    for i in range(1, n):
        Az[i + 1] = Az[i] + (z_eval - zs[i - 1]) * a[i] * Az[i - 1]
        Bz[i + 1] = Bz[i] + (z_eval - zs[i - 1]) * a[i] * Bz[i - 1]
        dAz[i + 1] = dAz[i] + a[i] * Az[i - 1] + (z_eval - zs[i - 1]) * a[i] * dAz[i - 1]
        dBz[i + 1] = dBz[i] + a[i] * Bz[i - 1] + (z_eval - zs[i - 1]) * a[i] * dBz[i - 1]

    # Compute the derivative of the Pade approximant
    dpade_value = dAz[n] / Bz[n] - Az[n] * dBz[n] / (Bz[n] * Bz[n])

    # Optional: print debugging information
    # print("Bz(n):", Bz[n])
    # if np.isclose(Bz[n].real, 0.0) and np.isclose(Bz[n].imag, 0.0):
    #     print("Warning: Bz(n) is close to zero:",
    return dpade_value


def calculate_pade_a(zs: np.ndarray, f_zs: np.ndarray) -> np.ndarray:
    """
    Calculate Pade coefficients.

    Args:
        zs: Input array of complex numbers.
        f_zs: Input array of complex numbers.

    Returns: Array of complex Pade coefficients
    """
    if (n := len(zs)) != len(f_zs):
        raise ValueError(f"{len(zs)=} != {len(f_zs)=}")

    # Initialize the g array
    g = np.zeros((n, n), dtype=complex)
    g[0, :n] = f_zs[:n]

    # Compute the divided differences
    for i in range(1, n):
        for j in range(i, n):
            #if np.real(g[i-1, j]) == 0.0 and np.imag(g[i-1, j]) == 0.0:
            #    print(f"g_i(z_j): i={i+1}, j={j+1}, g={g[i, j]}")
            g[i, j] = (g[i-1, i-1] - g[i-1, j]) / ((zs[j] - zs[i-1]) * g[i-1, j])

    # Extract the coefficients a(i)
    a = np.diag(g[:n, :n])
    # Optional: print coefficients for debugging
    # print('a:', a)
    return a
