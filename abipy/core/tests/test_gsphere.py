"""Tests for gsphere"""
import numpy as np

from abipy.core import Mesh3D
from abipy.core.gsphere import *
from abipy.core.testing import AbipyTest


class TestGSphere(AbipyTest):
    """Unit tests for GSphere"""

    def test_base(self):
        """Basic G-sphere methods"""
        ecut = 2
        lattice = np.array([1.,0,0, 0,1,0, 0,0,1])
        lattice.shape = (3, 3)
        kpoint = [0, 0, 0]
        gvecs = np.array([[0, 0, 0], [1, 0, 0]])

        gsphere = GSphere(ecut, lattice, kpoint, gvecs, istwfk=1)
        repr(gsphere)
        str(gsphere)

        assert len(gsphere) == 2
        assert [1, 0, 0] in gsphere
        assert gsphere.index([1, 0, 0]) == 1
        assert gsphere.count([1, 0, 0]) == 1

        self.serialize_with_pickle(gsphere, protocols=[-1])

        same_gsphere = gsphere.copy()
        assert gsphere == same_gsphere
        same_gsphere.kpt = [0.5, 0.1, 0.3]
        assert np.all(gsphere.kpoint.frac_coords == kpoint)

        assert np.all(gsphere.zeros() == 0)
        values = gsphere.czeros()
        assert values.dtype == complex
        assert np.all(values == 0)

        assert len(gsphere.empty()) == len(gsphere)
        assert len(gsphere.cempty()) == len(gsphere)

    def test_fft(self):
        """FFT transforms"""
        rprimd = np.array([1.,0,0, 0,1,0, 0,0,1])
        rprimd.shape = (3,3)

        mesh = Mesh3D( (12,3,5), rprimd)

        extra_dims = [(), 1, (2,), (3,4)]
        types = [float, complex]

        for exdim in extra_dims:
            for typ in types:
                fg = mesh.random(dtype=typ, extra_dims=exdim)

                fr = mesh.fft_g2r(fg)
                same_fg = mesh.fft_r2g(fr)
                self.assert_almost_equal(fg, same_fg)

                int_r = mesh.integrate(fr)
                int_g = fg[...,0,0,0]
                self.assert_almost_equal(int_r, int_g)
