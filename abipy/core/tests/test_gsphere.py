"""Tests for gsphere"""
from __future__ import print_function, division

import numpy as np

from abipy.core import Mesh3D
from abipy.core.gsphere import *
from abipy.tests import AbipyTest


class TestGSphere(AbipyTest):
    """Unit tests for GSphere"""

    def test_base(self):
        """Basic G-sphere methods"""

        ecut = 2
        gprimd = np.array([1.,0,0, 0,1,0, 0,0,1])
        gprimd.shape = (3,3)
        kpoint = [0,0,0]
        gvecs = np.array([0,0,0])

        gsphere = GSphere(ecut, gprimd, kpoint, gvecs, istwfk=1)
        print(gsphere)

        same_gsphere = gsphere.copy()
        self.assertTrue(gsphere == same_gsphere)
        same_gsphere.kpt = [0.5, 0.1, 0.3]
        self.assertTrue(np.all(gsphere.kpoint.frac_coords == kpoint))

        gsphere.zeros()
        gsphere.czeros()

        gsphere.empty()
        gsphere.cempty()

    def test_fft(self):
        """FFT transforms"""
        rprimd = np.array([1.,0,0, 0,1,0, 0,0,1])
        rprimd.shape = (3,3)

        mesh = Mesh3D( (12,3,5), rprimd)

        extra_dims = [(), 1, (2,), (3,4)]
        types = [np.float, np.complex]

        for exdim in extra_dims:
            for typ in types:
                fg = mesh.random(dtype=typ, extra_dims=exdim)

                fr = mesh.fft_g2r(fg)
                same_fg = mesh.fft_r2g(fr)
                self.assert_almost_equal(fg, same_fg)

                int_r = mesh.integrate(fr)
                int_g = fg[...,0,0,0]
                self.assert_almost_equal(int_r, int_g)
