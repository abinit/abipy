# coding: utf-8
"""This module contains the class defining Uniform 3D meshes."""

import numpy as np

#from itertools import product as iproduct
from monty.functools import lazy_property
#from numpy.random import random
from numpy.fft import fftn, ifftn, fftshift, ifftshift, fftfreq
from abipy.tools import duck


__all__ = [
    "Mesh3D",
]


class Mesh3D(object):
    r"""
    Descriptor-class for uniform 3D meshes.

    A Mesh3D object holds information on how functions, such
    as wave functions and electron densities, are discretized in a
    certain domain in space. The main information here is how many
    grid points are used in each direction of the unit cell.

    There are methods for tasks such as allocating arrays, performing
    symmetry operations and integrating functions over space.

    This is how a 2x2x2 3D array arr[x,y,z] is stored in memory::

        3-----7
        |\    |
        | \   |
        |  1-----5      z
        2--|--6  |   y  |
         \ |   \ |    \ |
          \|    \|     \|
           0-----4      +-----x

    """
    def __init__(self, shape, vectors):
        """
        Construct ``Mesh3D`` object.

        Args:
            shape: 3 int's Number of grid points along axes.
            vectors: unit cell vectors in real space.

        Attributes:

        ==========  ========================================================
        ``shape``   Array of the number of grid points along the three axes.
        ``dv``      Volume per grid point.
        ==========  ========================================================
        """
        self.shape = tuple(np.asarray(shape, int))
        self.size = np.prod(self.shape)
        self.vectors = np.reshape(vectors, (3, 3))

        cross12 = np.cross(self.vectors[1], self.vectors[2])
        self.dv = abs(np.sum(self.vectors[0] * cross12.T)) / self.size
        self.dvx = self.vectors[0] / self.nx
        self.dvy = self.vectors[1] / self.ny
        self.dvz = self.vectors[2] / self.nz

    def __len__(self):
        return self.size

    def __eq__(self, other):
        return (np.all(self.shape == other.shape) and
                np.all(self.vectors == other.vectors))

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return self.to_string()

    def __iter__(self):
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    yield ix * self.dvx + iy * self.dvy + iz * self.dvz

    def iter_ixyz_r(self):
        """
        Iterator returning (ixyz, rr) where ixyz gives the index of the point and r is the point on the grid.
        """
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    rr = ix * self.dvx + iy * self.dvy + iz * self.dvz
                    yield np.array((ix, iy, iz), dtype=int), rr

    def rpoint(self, ix, iy, iz):
        """The vector corresponding to the (ix, iy, iz) indices"""
        return ix * self.dvx + iy * self.dvy + iz * self.dvz

    def to_string(self, verbose=0):
        """String representation."""
        return self.__class__.__name__ + ": nx=%d, ny=%d, nz=%d" % self.shape

    @property
    def nx(self):
        """Number of points along x."""
        return self.shape[0]

    @property
    def ny(self):
        """Number of points along y."""
        return self.shape[1]

    @property
    def nz(self):
        """Number of points along z."""
        return self.shape[2]

    @lazy_property
    def inv_vectors(self):
        return np.linalg.inv(self.vectors)

    def _new_array(self, dtype=float, zero=True, extra_dims=()):
        shape = self.shape

        if duck.is_intlike(extra_dims):
            extra_dims = (extra_dims,)

        shape = extra_dims + tuple(shape)

        if zero:
            return np.zeros(shape, dtype)
        else:
            return np.empty(shape, dtype)

    def zeros(self, dtype=float, extra_dims=()):
        """
        Returns new zeroed 3D array for this domain.

        The type can be set with the ``dtype`` keyword.
        Extra dimensions can be added with ``extra_dims``.
        """
        return self._new_array(dtype=dtype, zero=True, extra_dims=extra_dims)

    def czeros(self, extra_dims=()):
        """Returns new zeroed 3D complex array for this domain."""
        return self._new_array(dtype=complex, zero=True, extra_dims=extra_dims)

    def empty(self, dtype=float, extra_dims=()):
        """
        Returns new uninitialized 3D |numpy-array| for this domain.

        The type can be set with the ``dtype`` keyword.
        Extra dimensions can be added with ``extra_dims``.
        """
        return self._new_array(dtype=dtype, zero=False, extra_dims=extra_dims)

    def cempty(self, extra_dims=()):
        """Returns new uninitialized 3D complex |numpy-array| for this domain."""
        return self._new_array(dtype=complex, zero=False, extra_dims=extra_dims)

    def random(self, dtype=float, extra_dims=()):
        """Returns random real |numpy-array| for this domain with val in [0.0, 1.0)."""
        shape = self.shape
        if duck.is_intlike(extra_dims):
            extra_dims = (extra_dims,)

        shape = extra_dims + tuple(shape)

        re = np.random.random(shape)
        if dtype == float:
            return re
        elif dtype == complex:
            im = self.random(extra_dims=extra_dims)
            return re + 1j*im
        else:
            raise ValueError("Wrong dtype: %s" % str(dtype))

    def crandom(self, extra_dims=()):
        """Returns random complex |numpy-array| for this domain with val in [0.0, 1.0)."""
        return self.random(dtype=complex, extra_dims=extra_dims)

    def reshape(self, arr):
        """
        Reshape the array arr defined on the FFT box.

        Returns |numpy-array| with 4 dimensions (?, nx, ny, nz) where ?*nx*ny*nz == arr.size
        """
        #if duck.is_intlike(extra_dims): extra_dims = (extra_dims,)
        #shape = extra_dims + self.shape)
        return np.reshape(arr, (-1,) + self.shape)

    def fft_r2g(self, fr, shift_fg=False):
        """
        FFT of array ``fr`` given in real space.
        """
        ndim, shape = fr.ndim, fr.shape

        if ndim == 1:
            fr = np.reshape(fr, self.shape)
            return self.fft_r2g(fr, shift_fg=shift_fg).flatten()

        elif ndim == 3:
            assert self.size == np.prod(shape[-3:])
            fg = fftn(fr)
            if shift_fg: fg = fftshift(fg)

        elif ndim > 3:
            assert self.size == np.prod(shape[-3:])
            axes = np.arange(ndim)[-3:]
            fg = fftn(fr, axes=axes)
            if shift_fg: fg = fftshift(fg, axes=axes)

        else:
            raise NotImplementedError("ndim < 3 are not supported")

        return fg / self.size

    def fft_g2r(self, fg, fg_ishifted=False):
        """
        FFT of array ``fg`` given in G-space.
        """
        ndim, shape = fg.ndim, fg.shape

        if ndim == 1:
            fg = np.reshape(fg, self.shape)
            return self.fft_g2r(fg, fg_ishifted=fg_ishifted).flatten()

        if ndim == 3:
            assert self.size == np.prod(shape[-3:])
            if fg_ishifted: fg = ifftshift(fg)
            fr = ifftn(fg)

        elif ndim > 3:
            assert self.size == np.prod(shape[-3:])
            axes = np.arange(ndim)[-3:]
            if fg_ishifted: fg = ifftshift(fg, axes=axes)
            fr = ifftn(fg, axes=axes)

        else:
            raise NotImplementedError("ndim < 3 are not supported")

        return fr * self.size

    #def fourier_interp(self, data, new_mesh, inspace="r"):
    #    """
    #    Fourier interpolation of data.

    #    Args:

    #        data: Input array defined on this mesh
    #        new_mesh: Mesh where data is interpolated
    #        inspace: string specifying if data is given in real space "r" or in reciprocal space "g".

    #    Return:
    #        Numpy array in real space on the new_mesh
    #    """
    #    raise NotImplementedError("gtransfer is missing")
    #    assert inspace in ("r", "g")

    #    # Insert data in the FFT box of new mesh.
    #    if inspace == "r": data = self.fft_r2g(data)
    #    intp_datag = new_mesh.gtransfer_from(self, data)

    #    # FFT transform G --> R.
    #    return new_mesh.fft_g2r(intp_datag)

    def integrate(self, fr):
        """
        Integrate array(s) fr.
        """
        shape, ndim = fr.shape, fr.ndim
        assert self.size == np.prod(shape[-3:])

        if ndim == 3:
            return fr.sum() * self.dv

        elif ndim > 3:
            sums = np.sum(np.reshape(fr, shape[:-3] + (-1,)), axis=-1)
            return sums * self.dv

        else:
            raise NotImplementedError("ndim < 3 are not supported")

    @lazy_property
    def gvecs(self):
        """
        Array with the reduced coordinates of the G-vectors.

        .. note::

            These are the G-vectors of the FFT box and should be used
            when we compute quantities on the FFT mesh.
            These vectors differ from the gvecs stored in |GSphere| that
            are k-centered and enclosed by a sphere whose radius is defined by ecut.
        """
        gx_list = np.rint(fftfreq(self.nx) * self.nx)
        gy_list = np.rint(fftfreq(self.ny) * self.ny)
        gz_list = np.rint(fftfreq(self.nz) * self.nz)
        #print(gz_list, gy_list, gx_list)

        gvecs = np.empty((self.size, 3), dtype=int)

        idx = -1
        for gx in gx_list:
            for gy in gy_list:
                for gz in gz_list:
                    idx += 1
                    gvecs[idx, :] = gx, gy, gz

        return gvecs

    @lazy_property
    def gmods(self):
        """[ng] |numpy-array| with :math:`|G|`"""
        gmet = np.dot(self.inv_vectors.T, self.inv_vectors)
        gmods = np.empty(self.size)
        for i, g in enumerate(self.gvecs):
            gmods[i] = np.dot(g, np.dot(gmet, g))

        return 2 * np.pi * np.sqrt(gmods)

    #@lazy_property
    #def gmax(self)
    #    return self.gmods.max()

    @lazy_property
    def rpoints(self):
        """|numpy-array| with the points in real space in reduced coordinates."""
        nx, ny, nz = self.nx, self.ny, self.nz
        rpoints = np.empty((self.size, 3))

        grid_points = np.meshgrid(np.linspace(0, 1, nx, endpoint=False),
                                  np.linspace(0, 1, ny, endpoint=False),
                                  np.linspace(0, 1, nz, endpoint=False))

        rpoints[:, 0] = grid_points[0].ravel()
        rpoints[:, 1] = grid_points[1].ravel()
        rpoints[:, 2] = grid_points[2].ravel()

        return rpoints

    #def ogrid_rfft(self):
    #    return np.ogrid[0:1:1/self.nx,
    #                    0:1:1/self.ny,
    #                    0:1:1/self.nz]

    #def line_inds(self, line):
    #    """
    #    Returns an ogrid with the indices associated to the specified line.

    #    Args:
    #        line:
    #            String specifying the type of line (in reduced coordinates),
    #            e.g. "x", "y", "z".
    #    """
    #    line = line.lower()

    #    if line == "x":
    #        return np.ogrid[0:self.nx, 0:1, 0:1]
    #    elif line == "y":
    #        return np.ogrid[0:1, 0:self.ny, 0:1]
    #    elif line == "z":
    #        return np.ogrid[0:1, 0:1, 0:self.nz]
    #    else:
    #        raise ValueError("Wrong line %s" % line)

    #def plane_inds(self, plane, h):
    #    """
    #    Returns an ogrid with the indices associated to the specified plane.

    #    Args:
    #        plane:
    #            String specifying the type of plane (in reduced coordinates),
    #            e.g. "xy" "xz" ...
    #        h:
    #            Index giving the position of the plane along the perpendicular.
    #    """
    #    plane = plane.lower()
    #    nx, ny, nz = self.nx, self.ny, self.nz

    #    if plane in ("xy", "yx"):
    #        return np.ogrid[0:nx, 0:ny, h:h+1]
    #    elif plane in ("xz", "zx"):
    #        return np.ogrid[0:nx, h:h+1, 0:nz]
    #    elif plane in ("yz", "zy"):
    #        return np.ogrid[h:h+1, 0:ny, 0:nz]
    #    else:
    #        raise ValueError("Wrong plane %s" % plane)

    #def irottable(self, symmops):
    #    nsym = len(symmops)
    #    nx, ny, nz = self.nx, self.ny, self.nz

    #    red2fft = np.diag([nx, ny, nz])
    #    fft2red = np.diag([1/nx, 1/ny, 1/nz])

    #    # For a fully compatible mesh, each mat in rotsm1_fft should be integer
    #    rotsm1_fft, tnons_fft = np.empty((nsym, 3, 3)), np.empty((nsym, 3))

    #    for isym, symmop in enumerate(symmops):
    #        rotm1_r, tau = symmop.rotm1_r, symmop.tau
    #        rotsm1_fft[isym] = np.dot(np.dot(red2fft, rotm1_r), fft2red)
    #        tnons_fft[isym] = np.dot(red2fft, tau)

    #    # Indeces of $R^{-1}(r-\tau)$ in the FFT box.
    #    irottable = np.empty((nsym, nx*ny*nz), dtype=int)

    #    #max_err = 0.0
    #    nxyz = np.array((nx, ny, nz), np.int)
    #    for isym in range(nsym):
    #        rm1_fft = rotsm1_fft[isym]
    #        tau_fft = tnons_fft[isym]
    #        for ifft, p1_fft in enumerate(iproduct(range(nx), range(ny), range(nz))):
    #            # Form R^-1 (r-\tau) in the FFT basis.
    #            p1_fft = np.array(p1_fft)
    #            prot_fft = np.dot(rm1_fft, p1_fft - tau_fft)
    #            #err = ABS(prot_fft - (ix, iy, iz)) / (nx, ny, nz)
    #            prot_fft = np.round(prot_fft)
    #            jx, jy, jz = prot_fft % nxyz
    #            irottable[isym, ifft] = jz + (jy * nz) + (jx * nx * ny)

    #    # Test
    #    #for isym in range(nsym):
    #    #    irottable[isym, ifft]
    #    #    rm1_fft = rotsm1_fft[isym]
    #    #    tau_fft = tnons_fft[isym]
    #    #    for ifft, p1_fft in enumerate(itertools.product(range(nx), range(ny), range(nz))):
    #    #        irot_fft == irottable[isym, ifft]

    #    return irottable

    def i_closest_gridpoints(self, points):
        """
        Given a list of points, this function return a |numpy-array| with the indices of the closest gridpoint.
        """
        points = np.reshape(points, (-1, 3))
        inv_vectors = self.inv_vectors
        fcoords = [np.dot(point, inv_vectors) for point in points]
        inds = [[np.mod(int(np.rint(pc[ii]*self.shape[ii])), self.nx) for ii in range(3)] for pc in fcoords]
        return np.array(inds)

        # return [(int(np.rint(pc[ii]*self.nx)), int(np.rint(pc[0]*self.nx)),int(np.rint(pc[0]*self.nx))) for pc in coords]
        # ix = int(np.rint(coords[0]*self.nx))
        # iy = int(np.rint(coords[1]*self.ny))
        # iz = int(np.rint(coords[2]*self.nz))
        # return (ix, iy, iz)

    # @DW TODO Add test.
    def dist_gridpoints_in_spheres(self, points, radius):
        # c_ab = np.cross(self.vectors[0], self.vectors[1])
        # c_bc = np.cross(self.vectors[1], self.vectors[2])
        # c_ca = np.cross(self.vectors[2], self.vectors[0])
        # h_ab = np.abs(np.dot(c_ab, self.vectors[2]) / np.linalg.norm(c_ab))
        # h_bc = np.abs(np.dot(c_bc, self.vectors[0]) / np.linalg.norm(c_bc))
        # h_ca = np.abs(np.dot(c_ca, self.vectors[1]) / np.linalg.norm(c_ca))
        maxdiag = max([np.linalg.norm(self.dvx+self.dvy+self.dvz),
                       np.linalg.norm(self.dvx+self.dvy-self.dvz),
                       np.linalg.norm(self.dvx-self.dvy+self.dvz),
                       np.linalg.norm(self.dvx-self.dvy-self.dvz)])
        c_ab = np.cross(self.dvx, self.dvy)
        c_bc = np.cross(self.dvy, self.dvz)
        c_ca = np.cross(self.dvz, self.dvx)
        h_ab = np.abs(np.dot(c_ab, self.dvz) / np.linalg.norm(c_ab))
        h_bc = np.abs(np.dot(c_bc, self.dvx) / np.linalg.norm(c_bc))
        h_ca = np.abs(np.dot(c_ca, self.dvy) / np.linalg.norm(c_ca))
        a_factor = 1.01 * (radius+0.5*maxdiag) / h_bc
        b_factor = 1.01 * (radius+0.5*maxdiag) / h_ca
        c_factor = 1.01 * (radius+0.5*maxdiag) / h_ab
        # print('HEIGHTS')
        # print(h_ab, h_bc, h_ca)
        # print(c_factor*h_ab, a_factor* h_bc, b_factor*h_ca)
        mins = np.array(np.floor([-a_factor, -b_factor, -c_factor]), dtype=int)
        maxes = np.array(np.ceil([a_factor, b_factor, c_factor]), dtype=int)
        # print('AAAMINS AND MAXES', mins, maxes)
        # maxes = np.ceil(nmax)
        i_closest_gridpoint_points = self.i_closest_gridpoints(points=points)
        dist_gridpoints_points = []
        r2 = radius**2
        for ipoint, point_i_closest_gridpoint in enumerate(i_closest_gridpoint_points):
            pp = points[ipoint]
            dist_gridpoints = []
            for ix in range(int(mins[0]), int(maxes[0])):
                ipx = point_i_closest_gridpoint[0] + ix
                for iy in range(int(mins[1]), int(maxes[1])):
                    ipy = point_i_closest_gridpoint[1] + iy
                    for iz in range(int(mins[2]), int(maxes[2])):
                        ipz = point_i_closest_gridpoint[2] + iz
                        gp = ipx*self.dvx + ipy * self.dvy + ipz * self.dvz
                        dist2_gp_pp = np.dot(pp-gp, pp-gp)
                        if dist2_gp_pp <= r2:
                            dist_gridpoints.append(((np.mod(ipx, self.nx), np.mod(ipy, self.ny), np.mod(ipz, self.nz)),
                                                    np.sqrt(dist2_gp_pp), (ipx, ipy, ipz)))
            dist_gridpoints_points.append(dist_gridpoints)
        return dist_gridpoints_points

    # def dist2_gridpoints_in_spheres(self, points, radius):
    #     # c_ab = np.cross(self.vectors[0], self.vectors[1])
    #     # c_bc = np.cross(self.vectors[1], self.vectors[2])
    #     # c_ca = np.cross(self.vectors[2], self.vectors[0])
    #     # h_ab = np.abs(np.dot(c_ab, self.vectors[2]) / np.linalg.norm(c_ab))
    #     # h_bc = np.abs(np.dot(c_bc, self.vectors[0]) / np.linalg.norm(c_bc))
    #     # h_ca = np.abs(np.dot(c_ca, self.vectors[1]) / np.linalg.norm(c_ca))
    #     maxdiag = max([np.linalg.norm(self.dvx+self.dvy+self.dvz),
    #                    np.linalg.norm(self.dvx+self.dvy-self.dvz),
    #                    np.linalg.norm(self.dvx-self.dvy+self.dvz),
    #                    np.linalg.norm(self.dvx-self.dvy-self.dvz)])
    #     c_ab = np.cross(self.dvx, self.dvy)
    #     c_bc = np.cross(self.dvy, self.dvz)
    #     c_ca = np.cross(self.dvz, self.dvx)
    #     h_ab = np.abs(np.dot(c_ab, self.dvz) / np.linalg.norm(c_ab))
    #     h_bc = np.abs(np.dot(c_bc, self.dvx) / np.linalg.norm(c_bc))
    #     h_ca = np.abs(np.dot(c_ca, self.dvy) / np.linalg.norm(c_ca))
    #     a_factor = 1.01 * (radius+0.5*maxdiag) / h_bc
    #     b_factor = 1.01 * (radius+0.5*maxdiag) / h_ca
    #     c_factor = 1.01 * (radius+0.5*maxdiag) / h_ab
    #     # print('HEIGHTS')
    #     # print(h_ab, h_bc, h_ca)
    #     # print(c_factor*h_ab, a_factor* h_bc, b_factor*h_ca)
    #     mins = np.array(np.floor([-a_factor, -b_factor, -c_factor]), dtype=int)
    #     maxes = np.array(np.ceil([a_factor, b_factor, c_factor]), dtype=int)
    #     # print('AAAMINS AND MAXES', mins, maxes)
    #     # maxes = np.ceil(nmax)
    #     i_closest_gridpoint_points = self.i_closest_gridpoints(points=points)
    #     dist_gridpoints_points = []
    #     r2 = radius**2
    #     for ipoint, point_i_closest_gridpoint in enumerate(i_closest_gridpoint_points):
    #         pp = points[ipoint]
    #         dist_gridpoints = []
    #         for ix in range(int(mins[0]), int(maxes[0])):
    #             ipx = point_i_closest_gridpoint[0] + ix
    #             for iy in range(int(mins[1]), int(maxes[1])):
    #                 ipy = point_i_closest_gridpoint[1] + iy
    #                 for iz in range(int(mins[2]), int(maxes[2])):
    #                     ipz = point_i_closest_gridpoint[2] + iz
    #                     gp = ipx*self.dvx + ipy * self.dvy + ipz * self.dvz
    #                     dist2_gp_pp = np.linalg.norm(pp-gp)
    #                     if dist_gp_pp <= r2:
    #                         dist_gridpoints.append(((np.mod(ipx, self.nx), np.mod(ipy, self.ny), np.mod(ipz, self.nz)),
    #                                                 dist_gp_pp, (ipx, ipy, ipz)))
    #         dist_gridpoints_points.append(dist_gridpoints)
    #     return dist_gridpoints_points
