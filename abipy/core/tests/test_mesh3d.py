"""Tests for mesh3d module"""
import numpy as np

from abipy.core.mesh3d import *
from abipy.core.testing import AbipyTest


class TestMesh3D(AbipyTest):
    """Test Mesh3d"""

    def test_base(self):
        """Testing mesh3d methods"""
        rprimd = np.reshape([1., 0, 0, 0, 1, 0, 0, 0, 1], (3, 3))

        mesh_443 = Mesh3D((4, 4, 3), rprimd)
        assert len(mesh_443) == 4 * 4 * 3
        assert mesh_443.nx == 4 and mesh_443.ny == 4 and mesh_443.nz == 3
        self.serialize_with_pickle(mesh_443)

        mesh_444 = Mesh3D((4, 4, 4), rprimd)
        repr(mesh_444); str(mesh_444)

        # Test __eq__
        assert mesh_443 == mesh_443
        assert mesh_443 != mesh_444

        self.assert_almost_equal(np.dot(mesh_444.vectors, mesh_444.inv_vectors), np.eye(3))
        gvecs = mesh_444.gvecs
        assert gvecs is mesh_444.gvecs
        assert gvecs.shape == (64, 3)
        self.assert_equal(gvecs[1], [0, 0, 1])

        gmods = mesh_444.gmods
        assert gmods.shape == mesh_444.size

        assert gmods[0] == 0
        self.assert_almost_equal(gmods[1], 2 * np.pi)

        rpoints = mesh_444.rpoints
        assert rpoints.shape == (64, 3)
        assert rpoints is mesh_444.rpoints

        empty = mesh_444.empty()
        assert empty.shape == mesh_444.shape and empty.dtype == float
        cempty = mesh_444.cempty()
        assert cempty.shape == mesh_444.shape and cempty.dtype == complex
        rand_vas = mesh_443.random()
        assert rand_vas.shape == mesh_443.shape and rand_vas.dtype == float
        crand_vas = mesh_443.crandom()
        assert crand_vas.shape == mesh_443.shape and crand_vas.dtype == complex

        zeros = mesh_444.zeros()
        assert zeros.shape == mesh_444.shape and np.all(zeros == 0) and zeros.dtype == float

        czeros = mesh_444.czeros()
        assert czeros.shape == mesh_444.shape and np.all(czeros == 0) and czeros.dtype == complex

        # Iteration
        # i = iz + iy * nz + ix * ny * nz
        #ny, nz = mesh_443.ny, mesh_443.nz
        #nyz = mesh_443.ny * nz
        #for i, r in enumerate(mesh_443):
        #    #iyz, ix = divmod(i, nyz)
        #    ix, iyz = divmod(i, nyz)
        #    iy, iz = divmod(iyz, ny)
        #    assert np.all(mesh_443.rpoint(ix, iy, iz) == r)

        for ixyz, r in mesh_443.iter_ixyz_r():
            self.assert_equal(mesh_443.rpoint(*ixyz), r)

        shift = 0.2 * mesh_443.dvx + 0.1 * mesh_443.dvy + 0.3 * mesh_443.dvz
        for ix in range(mesh_443.nx):
            for iy in range(mesh_443.ny):
                for iz in range(mesh_443.nz):
                    r = mesh_443.rpoint(ix, iy, iz)
                    self.assert_equal(mesh_443.i_closest_gridpoints(r), [[ix, iy, iz]])
                    r += shift
                    self.assert_equal(mesh_443.i_closest_gridpoints(r), [[ix, iy, iz]])

    def test_fft(self):
        """Test FFT transforms with mesh3d"""
        rprimd = np.array([1.,0,0, 0,1,0, 0,0,1])
        rprimd.shape = (3, 3)
        mesh = Mesh3D( (12, 3, 5), rprimd)

        extra_dims = [(), 1, (2,), (3, 1)]
        types = [float, complex]

        for exdim in extra_dims:
            for typ in types:
                fg = mesh.random(dtype=typ, extra_dims=exdim)

                fr = mesh.fft_g2r(fg)
                same_fg = mesh.fft_r2g(fr)
                self.assert_almost_equal(fg, same_fg)

                int_r = mesh.integrate(fr)
                int_g = fg[..., 0, 0, 0]
                self.assert_almost_equal(int_r, int_g)

    #def test_trilinear_interp(self):
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
