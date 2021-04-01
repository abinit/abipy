"""Tests for func1d module"""
import numpy as np
import collections
import tempfile

from abipy.core.testing import AbipyTest
from abipy.core.func1d import *


class TestFunction1D(AbipyTest):
    """Test Function1d."""

    def setUp(self):
        self.xmin, self.xmax = 0, 2 * np.pi
        mesh, self.h = np.linspace(self.xmin, self.xmax, 500, retstep=True)
        self.sinf = Function1D.from_func(np.sin, mesh)
        self.cosf = Function1D.from_func(np.cos, mesh)
        self.eix = self.cosf + 1j * self.sinf

    def test_base(self):
        """Basic methods of TestFunction1D."""
        sinf, cosf, eix = self.sinf, self.cosf, self.eix

        repr(sinf); str(sinf)
        assert isinstance(sinf, collections.abc.Iterable)
        assert len(sinf) == len(sinf.mesh)
        assert self.h == sinf.h

        cumint = cosf.integral()
        self.assert_almost_equal(cumint.values, sinf.values, decimal=5)

        self.assert_almost_equal(cosf.spline_integral(), 0.0)
        self.assert_almost_equal(cosf.spline.roots(), [np.pi/2, 3*np.pi/2])
        assert cosf.spline_on_mesh(cosf.mesh) == cosf

        _, path = tempfile.mkstemp(text=True)
        cosf.to_file(path)
        newcosf = Function1D.from_file(path)
        assert cosf == newcosf

        assert cosf == cosf
        assert not cosf == sinf
        assert cosf.has_same_mesh(sinf)
        assert cosf - cosf == Function1D(sinf.mesh, np.zeros(len(cosf)))

        assert 2 * cosf == Function1D(cosf.mesh, 2*cosf.values)
        assert cosf / 2. == Function1D(cosf.mesh, cosf.values / 2.)

        twof = Function1D(cosf.mesh, 2 * np.ones(len(cosf)))
        assert cosf/2. == cosf / twof
        assert cosf*2. == cosf * twof
        assert 1 - cosf ** 2 == sinf ** 2
        self.assert_almost_equal((1-cosf).integral()[-1], 2*np.pi)
        self.assert_almost_equal(sinf.l1_norm, 4., decimal=4)
        self.assert_almost_equal(sinf.integral_value, 0.0, decimal=4)

        one = Function1D.from_constant(cosf.mesh, 1.0)
        assert one == sinf ** 2 + cosf ** 2

        def primitive(x):
            return 0.5 * x - 0.25 * np.sin(2*x)

        self.assert_almost_equal(sinf.l2_norm**2, primitive(self.xmax) - primitive(self.xmin))

        # Derivatives
        assert sinf.finite_diff(order=1, acc=4) == cosf
        assert sinf.finite_diff(order=2, acc=4) == -sinf
        assert sinf.finite_diff(order=3, acc=4) == -cosf
        assert np.abs((sinf.finite_diff(order=4, acc=4) - sinf).max) < 1.e-5

        assert not sinf.iscomplexobj
        assert eix.iscomplexobj
        self.assert_almost_equal(cosf.max, 1.)
        self.assert_almost_equal(abs(sinf).min, 0)
        self.assert_almost_equal(eix.max, 1.)
        self.assert_almost_equal(eix.min, 1.)

        # Convolution
        #sg = sinf.gaussian_convolution(0.01)
        #self.assert_almost_equal(sg.integral_value, sinf.integral_value, decimal=3)
        #sl = sinf.lorentzian_convolution(0.01)
        #self.assert_almost_equal(sl.integral_value, sinf.integral_value, decimal=3)

        # Test Kramers-Kronig methods.
        # TODO: This is not a response function. Should use realistic values.
        #for with_div in (True, False):
        #    real_part = eix.real_from_kk(with_div=with_div)
        #    imag_part = real_part.imag_from_kk(with_div=with_div)

        if self.has_matplotlib():
            cosf.plot(show=False)
            eix.plot(show=False)
            eix.plot(cplx_mode="re", exchange_xy=True, xfactor=2, yfactor=3, show=False)

        if self.has_plotly():
            from abipy.tools.plotting import get_fig_plotly
            fig, go = get_fig_plotly()
            eix.plotly_traces(fig, cplx_mode="re", exchange_xy=True, xfactor=2, yfactor=3)

    def test_fft(self):
        """Test FFT transforms."""
        sinf, cosf, eix = self.sinf, self.cosf, self.eix
        # Test linearity
        assert 2 * sinf.fft() == (2 * sinf).fft()
        assert sinf.fft() + cosf.fft() == (sinf + cosf).fft()
        # Test F F^{-1} = I
        same_sinf = sinf.fft().ifft()
        self.assert_almost_equal(same_sinf.values, sinf.values)
        self.assert_almost_equal(same_sinf.mesh, sinf.mesh)
