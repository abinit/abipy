"""Tests for pwwave module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from abipy.core import Mesh3D
from abipy.core.testing import AbipyTest
from abipy.waves.pwwave import *


class TestPWWave(AbipyTest):
    """Test PWWave"""

    def test_base(self):
        """Basic tests for PWWave"""
        vectors = np.array([1.,0,0, 0,1,0, 0,0,1])
        vectors.shape = (3, 3)

        mesh_443 = Mesh3D((4, 4, 3), vectors)
        mesh_444 = Mesh3D((4, 4, 4), vectors)
        repr(mesh_444); str(mesh_444)
        assert not mesh_443 == mesh_444
        #mesh_444.get_gvec()
        #mesh_444.get_rpoints()

    def test_fft(self):
        """FFT transforms"""
        vectors = np.array([1.,0,0, 0,1,0, 0,0,1])
        vectors.shape = (3, 3)

        mesh = Mesh3D((12, 3, 5), vectors)
        extra_dims = [(), 1, (2,), (3,4)]
        types = [np.float, np.complex]

        for exdim in extra_dims:
            for typ in types:
                fg = mesh.random(dtype=typ, extra_dims=exdim)

                fr = mesh.fft_g2r(fg)
                same_fg = mesh.fft_r2g(fr)
                self.assert_almost_equal(fg, same_fg)

                int_r = mesh.integrate(fr)
                int_g = fg[..., 0, 0, 0]
                self.assert_almost_equal(int_r, int_g)
