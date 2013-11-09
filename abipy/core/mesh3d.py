"""This module contains the class defining Uniform 3D meshes."""
from __future__ import division, print_function

import numpy as np

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
        self.dv = abs(np.sum( self.vectors[0] * cross12.T)) / self.size

    def __len__(self):
        return self.size

    def __eq__(self, other):
        return (np.all(self.shape == other.shape) and
                np.all(self.vectors == other.vectors))

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        s = self.__class__.__name__
        s += ": nx=%d, ny=%d, nz=%d" % self.shape
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
            im = self.random(extra_dims = extra_dims)
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
        return np.reshape(arr, (-1,) + self.shape)

    #def zero_pad(self, arr):
    #  """
    #  Pad array with zeros as first element along non-periodic directions.
    #  """
    #  assert np.all(arr.shape[-3:] == (self.N_c + self.pbc - 1))
    #  if self.pbc.all(): return arr
    #  npbx, npby, npbz = 1 - self.pbc
    #  b_xg = np.zeros(arr.shape[:-3] + tuple(self.N_c), dtype=arr.dtype)
    #  b_xg[..., npbx:, npby:, npbz:] = arr
    #  return b_xg

    def fft_r2g(self, fr, shift_fg=False):
        """FFT of array fr given in real space."""
        ndim, shape = fr.ndim, fr.shape
        assert self.size == np.prod(shape[-3:])

        if ndim == 3:
            fg = fftn(fr)
            if shift_fg: fg = fftshift(fg)

        elif ndim > 3:
            axes = np.arange(ndim)[-3:]
            fg = fftn(fr, axes=axes)
            if shift_fg: fg = fftshift(fg, axes=axes)

        else:
            raise NotImplementedError("ndim < 3 are not supported")

        return fg / self.size

    def fft_g2r(self, fg, fg_ishifted=False):
        """FFT of array fg given in G-space."""
        ndim, shape  = fg.ndim, fg.shape
        assert self.size == np.prod(shape[-3:])

        if ndim == 3:
            if fg_ishifted: fg = ifftshift(fg)
            fr = ifftn(fg)

        elif ndim > 3:
            axes = np.arange(ndim)[-3:]
            if fg_ishifted: fg = ifftshift(fg, axes=axes)
            fr = ifftn(fg, axes=axes)

        else:
            raise NotImplementedError("ndim < 3 are not supported")

        return fr * self.size

    def integrate(self, fr):
        """Integrate array(s) fr."""
        shape, ndim = fr.shape, fr.ndim
        assert self.size == np.prod(shape[-3:])

        if ndim == 3:
            return fr.sum() * self.dv

        elif ndim > 3 :
            sums = np.sum(np.reshape(fr, shape[:-3] + (-1,)), axis=-1)
            return sums * self.dv

        else:
            raise NotImplementedError("ndim < 3 are not supported")
    
    def get_rpoints(self):
        rpoints = np.empty((self.size,3))
        idx = -1
        for x in range(self.nx):
            for y in range(self.ny):
                for z in range(self.nz):
                    idx += 1
                    rpoints[idx,0] = x / self.nx
                    rpoints[idx,1] = y / self.ny
                    rpoints[idx,2] = z / self.nz

        return rpoints

    #def get_rgrid(self):
    #  grid_x, grid_y, grid_z = np.mgrid[0:1:100j, 0:1:200j]

    def get_gvecs(self):
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
                    gvecs[idx,0] = gx
                    gvecs[idx,1] = gy
                    gvecs[idx,2] = gz

        return gvecs

    def trilinear_interp(self, fr, xx):
        """Interpolate fr on points."""
        raise NotImplementedError()
        fr = self.reshape(fr)

        xx = np.atleast_2d(xx)
        nx = len(xx)

        oarr = np.empty(xx, fr.dtype)

        #So we assume your grid has the corners (0.0, 0.0, 0.0) and (max_x, max_y, max_z)
        #and is aligned with the coordinate system.
        #We denote the number of cells along each axis by (n_x, n_y, n_z) respectively
        #and the point you wish to evaluate at by (x, y, z) (all of type float).
        #Then your logic might be something similar to

        #a_x = x * n_x / max_x
        #a_y = y * n_y / max_y
        #a_z = z * n_z / max_z
        #i_x = math.floor(a_x)
        #i_y = math.floor(a_y)
        #i_z = math.floor(a_z)
        #l_x = a_x - i_x
        #l_y = a_y - i_y
        #l_z = a_z - i_z

        #for idx, pt in enumerate(xx):
            #pass
            # Find the indices of nodes enclosing the point.

            # Linear interpolation
            #new[idx] = f..

        #if npts > 1:
        #    return new
        #else:
        #    return new[0]
