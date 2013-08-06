"""Tests for func1d module"""
from __future__ import print_function, division

import numpy as np
import collections
import cStringIO as StringIO

from abipy.core.func1d import *
from abipy.core.testing import *

class TestFunction1D(AbipyTest):
    """Test Function1d."""

    def setUp(self):
        self.xmin, self.xmax = 0, 2*np.pi
        mesh, self.h = np.linspace(self.xmin, self.xmax, 500, retstep=True)
        self.sinf = Function1D.from_func(np.sin, mesh)
        self.cosf = Function1D.from_func(np.cos, mesh)
        self.eix = self.cosf + 1j*self.sinf

    def test_base(self):
        """Basic methods of TestFunction1D."""
        sinf, cosf, eix = self.sinf, self.cosf, self.eix

        print(sinf)
        self.assertTrue(isinstance(sinf, collections.Iterable))
        self.assertTrue(len(sinf) == len(sinf.mesh))

        self.assertTrue(self.h == sinf.h)

        cumint = cosf.integral()
        self.assert_almost_equal(cumint.values, sinf.values, decimal=5)

        self.assert_almost_equal(cosf.spline_integral(), 0.0)
        self.assert_almost_equal(cosf.spline.roots(), [np.pi/2, 3*np.pi/2])

        stream = StringIO.StringIO()
        cosf.to_file(stream)
        stream.seek(0)
        newcosf = Function1D.from_file(stream)
        self.assertTrue(cosf == newcosf)

        self.assertTrue(cosf == cosf)
        self.assertFalse(cosf == sinf)
        self.assertTrue(cosf.has_same_mesh(sinf))
        self.assertTrue(cosf-cosf, Function1D(sinf.mesh, np.zeros(len(cosf))))

        self.assertTrue(2*cosf == Function1D(cosf.mesh, 2*cosf.values))
        self.assertTrue(cosf/2. == Function1D(cosf.mesh, cosf.values/2.))

        twof = Function1D(cosf.mesh, 2 * np.ones(len(cosf)))
        self.assertTrue(cosf/2. == cosf / twof)
        self.assertTrue(cosf*2. == cosf * twof)
        self.assertTrue(1 - cosf**2 == sinf**2)
        self.assert_almost_equal((1-cosf).integral()[-1], 2*np.pi)
        self.assert_almost_equal(sinf.l1_norm, 4., decimal=4)

        def primitive(x):
            return 0.5 * x - 0.25 * np.sin(2*x)
        self.assert_almost_equal(sinf.l2_norm, primitive(self.xmax) - primitive(self.xmin))

        # Derivatives
        self.assertTrue(sinf.finite_diff(order=1, acc=4) == cosf)
        self.assertTrue(sinf.finite_diff(order=2, acc=4) == -sinf)
        self.assertTrue(sinf.finite_diff(order=3, acc=4) == -cosf)
        self.assertTrue(np.abs((sinf.finite_diff(order=4, acc=4) - sinf).max) < 1.e-5)

        self.assertFalse(sinf.iscomplexobj)
        self.assertTrue(eix.iscomplexobj)
        self.assert_almost_equal(cosf.max, 1.)
        self.assert_almost_equal(abs(sinf).min, 0)
        self.assert_almost_equal(eix.max, 1.)
        self.assert_almost_equal(eix.min, 1.)

    def test_fft(self):
        """Test FFT transforms."""
        sinf, cosf, eix = self.sinf, self.cosf, self.eix
        # Test linearity
        self.assertTrue(2*sinf.fft() == (2*sinf).fft())
        self.assertTrue(sinf.fft() + cosf.fft() == (sinf + cosf).fft())
        # Test F F^{-1} = I
        same_sinf = sinf.fft().ifft()
        self.assert_almost_equal(same_sinf.values, sinf.values)
        self.assert_almost_equal(same_sinf.mesh, sinf.mesh)

if __name__ == "__main__": 
    import unittest
    unittest.main()
