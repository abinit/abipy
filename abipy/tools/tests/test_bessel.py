# coding: utf-8
"""Tests for bessel module."""
import numpy as np

from abipy.tools import bessel
from abipy.core.testing import AbipyTest


class TestBessels(AbipyTest):

    def test_spline_int_jlqr(self):
        """Testing spline_int_jlqr."""
        rcut, qmax = 1.3, 3
        l = 0
        qvals = np.linspace(0, qmax, num=200)
        spline = bessel.spline_int_jlqr(l, qmax, rcut, numq=1024, numr=1024)
        fq = spline(qvals)
        assert len(fq) == len(qvals)
        assert fq.shape == (len(qvals),)
        self.assert_almost_equal(fq[0], rcut ** 3 / 3)
        def primitive(x):
            return -x * np.cos(x) + np.sin(x)
        self.assert_almost_equal(fq[-1], (1 / qmax**3) * (primitive(qmax*rcut) - primitive(0)))
