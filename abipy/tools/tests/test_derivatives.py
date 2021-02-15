# coding: utf-8
"""Tests for derivatives module."""
import numpy as np

from abipy.tools.derivatives import finite_diff
from abipy.core.testing import AbipyTest


class FiniteDiffTest(AbipyTest):

    def test_complex(self):
        """complex functions are not supported"""
        x, h = np.linspace(0, 1,  800, retstep=True)
        cf = 1j*x
        with self.assertRaises(ValueError):
            return finite_diff(cf, h)

    def test_poly(self):
        """Test derivatives of polynomials."""
        x, h = np.linspace(0, 1,  800, retstep=True)
        orders = [1,2,3,4]
        accuracies = [2,4,6]

        f = x**4
        dpolys = {
            1: 4 * x**3,
            2: 12 * x**2,
            3: 24 * x,
            4: 24 * np.ones(len(x)),
        }

        decs = {
            1: 5,
            2: 4,
            3: 4,
            4: 1,
        }

        for order in orders:
            for acc in accuracies:
                if order == 4 and acc == 6: continue
                #print("order %s, acc %s" % (order, acc))
                yder = finite_diff(f, h, order=order, acc=acc)
                #print(np.max(np.abs(yder - dpolys[order])))
                self.assert_almost_equal(yder, dpolys[order], decimal=decs[order])

    def test_exp(self):
        """Test derivatives of exp(x)."""
        x, h = np.linspace(0, 2, 800, retstep=True)
        orders = [1,2,3,4]
        accuracies = [2,4,6]
        exp = np.exp(x)

        decs = {
            1: 4,
            2: 4,
            3: 4,
            4: 2,
        }

        for order in orders:
            for acc in accuracies:
                if order == 4 and acc == 6: continue
                #print("order %s, acc %s" % (order, acc))
                yder = finite_diff(exp, h, order=order, acc=acc)
                #print(np.max(np.abs(yder - exp)))
                self.assert_almost_equal(yder, exp, decs[order])

                d = finite_diff(exp, h, order=order, acc=acc, index=100)
                assert yder[100] == d.value
