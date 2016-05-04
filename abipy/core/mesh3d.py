# coding: utf-8
"""This module contains the class defining Uniform 3D meshes."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
from itertools import product as iproduct

from monty.functools import lazy_property
from numpy.random import random
from numpy.fft import fftn, ifftn, fftshift, ifftshift, fftfreq

__all__ = [
    "Mesh3D",
]


class Mesh3D(object):
    """
    Descriptor-class for uniform 3D meshes.

    A Mesh3D object holds information on how functions, such
    as wave functions and electron densities, are discretized in a
    certain domain in space.  The main information here is how many
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
        Construct `Mesh3D` object.

        Args:
            shape:
                3 int's Number of grid points along axes.
            vectors:
                unit cell vectors

        Attributes:

        ==========  ========================================================
        ``shape``   Array of the number of grid points along the three axes.
        ``dv``      Volume per grid point.
        ==========  ========================================================
        """
        self.shape = tuple( np.asarray(shape, np.int) )
        self.size = np.prod(self.shape)
        self.vectors = np.reshape(vectors, (3,3))

        cross12 = np.cross(self.vectors[1], self.vectors[2])
        self.dv = abs(np.sum(self.vectors[0] * cross12.T)) / self.size

    def __len__(self):
        return self.size

    def __eq__(self, other):
        return (np.all(self.shape == other.shape) and
                np.all(self.vectors == other.vectors))

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return self.to_string()

    def to_string(self, prtvol=0):
        s = self.__class__.__name__ + ": nx=%d, ny=%d, nz=%d" % self.shape
        return s

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

    def _new_array(self, dtype=np.float, zero=True, extra_dims=()):

        shape = self.shape

        if isinstance(extra_dims, int): extra_dims = (extra_dims,)
        shape = extra_dims + tuple(shape)

        if zero:
            return np.zeros(shape, dtype)
        else:
            return np.empty(shape, dtype)

    def zeros(self, dtype=np.float, extra_dims=()):
        """
        Returns new zeroed 3D array for this domain.

        The type can be set with the `dtype` keyword.
        Extra dimensions can be added with `extra_dims`.
        """
        return self._new_array(dtype=dtype, zero=True, extra_dims=extra_dims)

    def czeros(self, extra_dims=()):
        """Returns new zeroed 3D complex array for this domain."""
        return self._new_array(dtype=np.complex, zero=True, extra_dims=extra_dims)

    def empty(self, dtype=np.float, extra_dims=()):
        """
        Returns new uninitialized 3D array for this domain.

        The type can be set with the `dtype` keyword.
        Extra dimensions can be added with `extra_dims`.
        """
        return self._new_array(dtype=dtype, zero=False, extra_dims=extra_dims)

    def cempty(self, extra_dims=()):
        """Returns new uninitialized 3D complex array for this domain."""
        return self._new_array(dtype=np.complex, zero=False, extra_dims=extra_dims)

    def random(self, dtype=np.float, extra_dims=()):
        """Returns random real array for this domain with val in [0.0, 1.0)."""
        shape = self.shape
        if isinstance(extra_dims, int): extra_dims = (extra_dims,)
        shape = extra_dims + tuple(shape)

        re = np.random.random(shape)
        if dtype == np.float:
            return re
        elif dtype == np.complex:
            im = self.random(extra_dims=extra_dims)
            return re + 1j*im
        else:
            raise ValueError("Wrong dtype = " + str(dtype))

    def crandom(self, extra_dims=()):
        """Returns random complex array for this domain with val in [0.0, 1.0)."""
        return random(self, dtype=np.complex, extra_dims=extra_dims)

    def reshape(self, arr):
        """
        Reshape the array arr defined on the FFT box.

        Returns ndarray with 4 dimensions (?,nx,ny,nz) where ?*nx*ny*nz == arr.size
        """
        #if isinstance(extra_dims, int): extra_dims = (extra_dims,)
        #shape = extra_dims + self.shape)
        return np.reshape(arr, (-1,) + self.shape)

    def fft_r2g(self, fr, shift_fg=False):
        """FFT of array fr given in real space."""
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
        """FFT of array fg given in G-space."""
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

    def fourier_interp(self, data, new_mesh, inspace="r"):
        """
        Fourier interpolation of data.

        :param data: Input array defined on this mesh
        :param new_mesh: Mesh where data is interpolated
        :param inspace: string specifying if data is given in real space "r" or in reciprocal space "g".
        :return: Numpy array in real space on the new_mesh
        """
        raise NotImplementedError("gtransfer is missing")
        assert inspace in ("r", "g")

        # Insert data in the FFT box of new mesh.
        if inspace == "r": data = self.fft_r2g(data)
        intp_datag = new_mesh.gtransfer_from(self, data)

        # FFT transform G --> R.
        return new_mesh.fft_g2r(intp_datag)

    def integrate(self, fr):
        """Integrate array(s) fr."""
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
        """Array with the reduced coordinates of the G-vectors."""
        gx_list = np.rint(fftfreq(self.nx) * self.nx)
        gy_list = np.rint(fftfreq(self.ny) * self.ny)
        gz_list = np.rint(fftfreq(self.nz) * self.nz)
        #print(gz_list, gy_list, gx_list)
                                                      
        gvecs = np.empty((self.size,3), dtype=np.int)
                                                      
        idx = -1
        for gx in gx_list:
            for gy in gy_list:
                for gz in gz_list:
                    idx += 1
                    gvecs[idx,:] = gx, gy, gz

        return gvecs

    @lazy_property
    def rpoints(self):
        """Array with the points in real space in reduced coordinates."""
        nx, ny, nz = self.nx, self.ny, self.nz
        rpoints = np.empty((self.size,3))

        for ifft, p1_fft in enumerate(iproduct(range(nx), range(ny), range(nz))):
            rpoints[ifft,:] = p1_fft[0]/nx, p1_fft[1]/ny, p1_fft[2]/nz

        return rpoints

    #def ogrid_rfft(self):
    #    return np.ogrid[0:1:1/self.nx, 
    #                    0:1:1/self.ny,
    #                    0:1:1/self.nz]

    def line_inds(self, line):
        """
        Returns an ogrid with the indices associated to the specified line.

        Args:
            line: 
                String specifying the type of line (in reduced coordinates),
                e.g. "x", "y", "z".
        """
        line = line.lower()

        if line == "x":
            return np.ogrid[0:self.nx, 0:1, 0:1]
        elif line == "y":
            return np.ogrid[0:1, 0:self.ny, 0:1]
        elif line == "z":
            return np.ogrid[0:1, 0:1, 0:self.nz]
        else:
            raise ValueError("Wrong line %s" % line)

    def plane_inds(self, plane, h):
        """
        Returns an ogrid with the indices associated to the specified plane.

        Args:
            plane: 
                String specifying the type of plane (in reduced coordinates),
                e.g. "xy" "xz" ...
            h:
                Index giving the position of the plane along the perpendicular.
        """
        plane = plane.lower()
        nx, ny, nz = self.nx, self.ny, self.nz

        if plane in ("xy", "yx"):
            return np.ogrid[0:nx, 0:ny, h:h+1]
        elif plane in ("xz", "zx"):
            return np.ogrid[0:nx, h:h+1, 0:nz]
        elif plane in ("yz", "zy"):
            return np.ogrid[h:h+1, 0:ny, 0:nz]
        else:
            raise ValueError("Wrong plane %s" % plane)

    def irottable(self, symmops):
        nsym = len(symmops)
        nx, ny, nz = self.nx, self.ny, self.nz

        red2fft = np.diag([nx, ny, nz])
        fft2red = np.diag([1/nx, 1/ny, 1/nz])

        # For a fully compatible mesh, each mat in rotsm1_fft should be integer 
        rotsm1_fft, tnons_fft = np.empty((nsym,3,3)), np.empty((nsym,3)) 

        for isym, symmop in enumerate(symmops):
            rotm1_r, tau = symmop.rotm1_r, symmop.tau
            rotsm1_fft[isym] = np.dot(np.dot(red2fft, rotm1_r), fft2red)
            tnons_fft[isym] = np.dot(red2fft, tau)

        # Indeces of $R^{-1}(r-\tau)$ in the FFT box.
        irottable = np.empty((nsym, nx*ny*nz), dtype=np.int) 

        #max_err = 0.0
        nxyz = np.array((nx, ny, nz), np.int)
        for isym in range(nsym):
            rm1_fft = rotsm1_fft[isym]
            tau_fft = tnons_fft[isym]
            for ifft, p1_fft in enumerate(iproduct(range(nx), range(ny), range(nz))):
                # Form R^-1 (r-\tau) in the FFT basis.
                p1_fft = np.array(p1_fft)
                prot_fft = np.dot(rm1_fft, p1_fft - tau_fft)
                #err = ABS(prot_fft - (ix, iy, iz)) / (nx, ny, nz)
                prot_fft = np.round(prot_fft)
                jx, jy, jz = prot_fft % nxyz
                irottable[isym, ifft] = jz + (jy * nz) + (jx * nx * ny)

        # Test
        #for isym in range(nsym):
        #    irottable[isym, ifft]
        #    rm1_fft = rotsm1_fft[isym]
        #    tau_fft = tnons_fft[isym]
        #    for ifft, p1_fft in enumerate(itertools.product(range(nx), range(ny), range(nz))):
        #        irot_fft == irottable[isym, ifft]

        return irottable
