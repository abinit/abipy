# coding: utf-8
"""
Functions to perform analytic continuation with Pade'
Some of these routines have been directly translated from the Fortran versions
implemented in ABINIT.
"""
from __future__ import annotations

import numpy as np


def pade(z, f, zz) -> complex:
    """
    Calculate the Pade approximant of the function f at zz.

    Args
      z (numpy.ndarray): Input array of complex numbers (shape: n).
      f (numpy.ndarray): Input array of complex numbers (shape: n).
      zz (complex): Point at which to evaluate the Pade approximant.

    Returns:
        complex: The Pade approximant value at zz.
    """
    if (n := len(z)) != len(f):
        raise ValueError(f"{len(z)=} != {len(f)=}")

    # Calculate Pade coefficients
    a = calculate_pade_a(z, f)

    # Initialize Az and Bz arrays
    Az = np.zeros(n + 1, dtype=complex)
    Bz = np.zeros(n + 1, dtype=complex)
    Az[0] = 0.0 + 0.0j  # czero
    Az[1] = a[0]
    Bz[0] = 1.0 + 0.0j  # cone
    Bz[1] = 1.0 + 0.0j  # cone

    # Calculate Az and Bz recursively
    for i in range(1, n):
        Az[i + 1] = Az[i] + (zz - z[i - 1]) * a[i] * Az[i - 1]
        Bz[i + 1] = Bz[i] + (zz - z[i - 1]) * a[i] * Bz[i - 1]

    # Compute the Pade approximant
    pade_value = Az[n] / Bz[n]

    # Optional: print debugging information
    # print("Pade approximant:", pade_value)
    # print("Bz(n):", Bz[n])
    # if np.isclose(Bz[n].real, 0.0) and np.isclose(Bz[n].imag, 0.0):
    #     print("Warning: Bz(n) is close to zero:", Bz[n])

    return pade_value


def dpade(z, f, zz) -> complex
    """
    Calculate the derivative of the Pade approximant of the function f at zz.

    Parameters:
    z (numpy.ndarray): Input array of complex numbers (shape: n).
    f (numpy.ndarray): Input array of complex numbers (shape: n).
    zz (complex): Point at which to evaluate the derivative of the Pade approximant.

    Returns:
        complex: The derivative of the Pade approximant value at zz.
    """
    n = len(z)
    # Calculate Pade coefficients
    a = calculate_pade_a(z, f)

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
        Az[i + 1] = Az[i] + (zz - z[i - 1]) * a[i] * Az[i - 1]
        Bz[i + 1] = Bz[i] + (zz - z[i - 1]) * a[i] * Bz[i - 1]
        dAz[i + 1] = dAz[i] + a[i] * Az[i - 1] + (zz - z[i - 1]) * a[i] * dAz[i - 1]
        dBz[i + 1] = dBz[i] + a[i] * Bz[i - 1] + (zz - z[i - 1]) * a[i] * dBz[i - 1]

    # Compute the derivative of the Pade approximant
    dpade_value = dAz[n] / Bz[n] - Az[n] * dBz[n] / (Bz[n] * Bz[n])

    # Optional: print debugging information
    # print("Bz(n):", Bz[n])
    # if np.isclose(Bz[n].real, 0.0) and np.isclose(Bz[n].imag, 0.0):
    #     print("Warning: Bz(n) is close to zero:",


def calculate_pade_a(z, f) -> np.ndarray:
    """
    Calculate Pade coefficients.

    Parameters:
    n (int): Number of points.
    z (numpy.ndarray): Input array of complex numbers (shape: n).
    f (numpy.ndarray): Input array of complex numbers (shape: n).

    Returns:
        numpy.ndarray: Array of complex Pade coefficients (shape: n).
    """
    if (n := len(z)) != len(f):
        raise ValueError(f"{len(z)=} != {len(f)=}")

    # Initialize the g array
    g = np.zeros((n, n), dtype=complex)
    g[0, :n] = f[:n]

    # Compute the divided differences
    for i in range(1, n):
        for j in range(i, n):
            #if np.real(g[i-1, j]) == 0.0 and np.imag(g[i-1, j]) == 0.0:
            #    print(f"g_i(z_j): i={i+1}, j={j+1}, g={g[i, j]}")
            g[i, j] = (g[i-1, i-1] - g[i-1, j]) / ((z[j] - z[i-1]) * g[i-1, j])

    # Extract the coefficients a(i)
    a = np.diag(g[:n, :n])
    # Optional: print coefficients for debugging
    # print('a:', a)
    return a
