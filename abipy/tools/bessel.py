# coding: utf-8
"""This module provides functions to compute integrals of Bessel functions."""
import numpy as np

from scipy.special import spherical_jn
from scipy.interpolate import UnivariateSpline
from scipy.integrate import simps # cumtrapz, quad


_DEFAULTS = {"numq": 3001, "numr": 3001}


def spline_int_jlqr(l, qmax, rcut, numq=None, numr=None):
    r"""
    Compute :math:`j_n(z) = \int_0^{rcut} r^2 j_l(qr) dr`
    where :math:`j_l` is the Spherical Bessel function.

    Args:
        l: Angular momentum
        qmax: Max :math:`|q|` in integral in Ang-1
        rcut: Sphere radius in Angstrom.
        numq: Number of q-points in qmesh.
        numr: Number of r-points for integration.

    Return:
        Spline object.
    """
    numq = numq if numq is not None else _DEFAULTS["numq"]
    numr = numr if numr is not None else _DEFAULTS["numr"]

    rs = np.linspace(0, rcut, num=numr)
    r2 = rs ** 2
    qmesh = np.linspace(0, qmax, num=numq)
    values = []
    for q in qmesh:
        qrs = q * rs
        #qrs = 2 * np.pi * q * rs
        ys = spherical_jn(l, qrs) * r2
        intg = simps(ys, x=rs)
        values.append(intg)

    return UnivariateSpline(qmesh, values, s=0)
