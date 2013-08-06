"""Tests for pwwave"""
from __future__ import print_function, division
import numpy as np

from abipy.core import Mesh3D
from abipy.waves.pwwave import *
from abipy.core.testing import *


class TestPWWave(AbipyTest):
    """Test PWWave"""

    def test_base(self):
        """basic tests"""
        rprimd = np.array([1.,0,0, 0,1,0, 0,0,1])
        rprimd.shape = (3,3)

        mesh_443 = Mesh3D( (4,4,3), rprimd)
        mesh_444 = Mesh3D( (4,4,4), rprimd)

        print(mesh_444)
        self.assertNotEqual(mesh_443, mesh_444)

        #mesh_444.get_gvec()
        #mesh_444.get_rpoints()

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
