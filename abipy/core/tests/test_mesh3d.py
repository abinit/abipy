"""Tests for mesh3d"""
from __future__ import print_function, division

import numpy as np

from abipy.core.mesh3d import *
from abipy.tests import AbipyTest

class TestMesh3D(AbipyTest):
    "Test Mesh3d"

    def test_base(self):
        """basic tools"""
        rprimd = np.array([1.,0,0, 0,1,0, 0,0,1])
        rprimd.shape = (3,3)

        mesh_443 = Mesh3D((4,4,3), rprimd)
        mesh_444 = Mesh3D((4,4,4), rprimd)
        #
        # Test __eq__
        self.assertNotEqual(mesh_443, mesh_444)

        print(mesh_444)

        mesh_444.get_gvecs()
        mesh_444.get_rpoints()

    def test_fft(self):
        """FFT transforms"""
        rprimd = np.array([1.,0,0, 0,1,0, 0,0,1])
        rprimd.shape = (3,3)

        mesh = Mesh3D( (12,3,5), rprimd)

        extra_dims = [(), 1, (2,), (3,1)]
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


    #def test_trilinear_interp(self):
    #    return
    #    rprimd = np.array([1.,0,0, 0,1,0, 0,0,1])
    #    rprimd.shape = (3,3)

    #    mesh = Mesh3D( (12,3,5), rprimd)

    #    mesh_points = mesh.get_rpoints()

    #    # Linear function.
    #    def lfunc(rr, c = (1,1,1,1)): # xyz, xy, xz, yz
    #        x = rr[:,0]
    #        y = rr[:,1]
    #        z = rr[:,2]
    #        return c[0] * (x*y*z) + c[1] * (x*y) + c[2] * (x*z) + c[3] * (y*z)

    #    fr = lfunc(mesh_points)

    #    xx = mesh_points
    #    same_fr = mesh.trilinear_interp(fr, xx)
    #    self.assert_almost_equal(fr, same_fr)

    #    #shift = 0.5 *
    #    xx = mesh_points + shift
    #    same_fr = mesh.trilinear_interp(fr, xx)
    #    self.assert_almost_equal(fr, same_fr)
