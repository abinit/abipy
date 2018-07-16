# coding: utf-8
"""This module contains the class defining the G-sphere for wavefunctions, densities and potentials"""
from __future__ import print_function, division, unicode_literals, absolute_import

import collections
import numpy as np

from .kpoints import Kpoint
from abipy.tools import duck


__all__ = [
    "GSphere",
]


class GSphere(collections.Sequence):
    """Descriptor-class for the G-sphere."""

    def __init__(self, ecut, lattice, kpoint, gvecs, istwfk=1):
        """
        Args:
            ecut: Cutoff energy in Hartree.
            lattice: Reciprocal lattice.
            kpoint: Reduced coordinates of the k-point.
            gvecs: Array with the reduced coordinates of the G-vectors.
            istwfk: Storage option (time-reversal symmetry, see abinit variable)
        """
        self.ecut = ecut
        self.lattice = lattice
        self.kpoint = Kpoint.as_kpoint(kpoint, lattice)

        self._gvecs = np.reshape(np.array(gvecs), (-1, 3))
        self.npw = self.gvecs.shape[0]

        self.istwfk = istwfk
        if istwfk != 1:
            raise NotImplementedError("istwfk %d is not implemented" % self.istwfk)

    @property
    def gvecs(self):
        """|numpy-array| with the G-vectors in reduced coordinates."""
        return self._gvecs

    #@property
    #def get_kpg2(self):
    #    """ndarray with |k+G|**2. in atomic unit"""
    #    # note that now we use pymatgen lattice hence we have to convert to a.u.
    #    self.kpg2 =
    #    return self._kpg2

    # Sequence protocol
    def __len__(self):
        return self.gvecs.shape[0]

    def __getitem__(self, slice):
        return self.gvecs[slice]

    def __iter__(self):
        return self.gvecs.__iter__()

    def __contains__(self, gvec):
        return np.asarray(gvec) in self.gvecs

    def index(self, gvec):
        """
        return the index of the G-vector ``gvec`` in self.
        Raises: `ValueError` if the value is not present.
        """
        gvec = np.asarray(gvec)
        for i, g in enumerate(self):
            if np.all(g == gvec):
                return i
        else:
            raise ValueError("Cannot find %s in Gsphere" % str(gvec))

    def count(self, gvec):
        """Return number of occurrences of gvec."""
        return np.count_nonzero(np.all(g == gvec) for g in self)

    def __str__(self):
        return self.to_string()

    def __eq__(self, other):
        if other is None: return False
        return (self.ecut == other.ecut and
                np.all(self.lattice == other.lattice) and
                self.kpoint == other.kpoint and
                np.all(self.gvecs == other.gvecs) and
                self.istwfk == other.istwfk)

    def __ne__(self, other):
        return not (self == other)

    def copy(self):
        """shallow copy."""
        return self.__class__(self.ecut, self.lattice.copy(), self.kpoint.copy(), self.gvecs.copy(), istwfk=self.istwfk)

    def to_string(self, verbose=0):
        """String representation."""
        name = str(self.__class__)
        s = name + ": kpoint: %(kpoint)s, ecut: %(ecut)f, npw: %(npw)d, istwfk: %(istwfk)d" % self.__dict__
        return s

    def _new_array(self, dtype=np.float, zero=True, extra_dims=()):
        """Returns a numpy array defined on the sphere."""
        shape = (self.npw,)

        if duck.is_intlike(extra_dims):
            extra_dims = (extra_dims,)

        shape = extra_dims + tuple(shape)

        if zero:
            return np.zeros(shape, dtype)
        else:
            return np.empty(shape, dtype)

    def zeros(self, dtype=np.float, extra_dims=()):
        """
        Returns new zeroed 1D |numpy-array|.

        The type can be set with the ``dtype`` keyword.
        Extra dimensions can be added with ``extra_dims``.
        """
        return self._new_array(dtype=dtype, zero=True, extra_dims=extra_dims)

    def czeros(self, extra_dims=()):
        """New zeroed 1D complex |numpy-array|."""
        return self._new_array(dtype=np.complex, zero=True, extra_dims=extra_dims)

    def empty(self, dtype=np.float, extra_dims=()):
        """
        Returns new uninitialized 1D |numpy-array|.

        The type can be set with the ``dtype`` keyword.
        Extra dimensions can be added with ``extra_dims``.
        """
        return self._new_array(dtype=dtype, zero=False, extra_dims=extra_dims)

    def cempty(self, extra_dims=()):
        """Returns new uninitialized 1D complex |numpy-array|."""
        return self._new_array(dtype=np.complex, zero=False, extra_dims=extra_dims)

    #def build_fftbox(self, boxsph_ratio=1.05):
    #  """Returns the number of divisions of the FFT box enclosing the sphere."""
    #  #return ndivs

    def tofftmesh(self, mesh, arr_on_sphere):
        """
        Insert the array ``arr_on_sphere`` given on the sphere inside the FFT mesh.

        Args:
            mesh:
            arr_on_sphere:
        """
        arr_on_sphere = np.atleast_2d(arr_on_sphere)
        ishape = arr_on_sphere.shape
        s0 = ishape[0]
        assert self.npw == ishape[-1]

        arr_on_mesh = np.zeros((s0,) + mesh.shape, dtype=arr_on_sphere.dtype)

        if self.istwfk == 1:
            #do ipw=1,npw
            #  i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
            #  i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
            #  i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
            #end do

            n1, n2, n3 = mesh.shape
            for sph_idx, gvec in enumerate(self.gvecs):
                i1 = gvec[0]
                if i1 < 0: i1 = i1 + n1
                i2 = gvec[1]
                if i2 < 0: i2 = i2 + n2
                i3 = gvec[2]
                if i3 < 0: i3 = i3 + n3
                arr_on_mesh[..., i1, i2, i3] = arr_on_sphere[..., sph_idx]

        else:
            raise NotImplementedError("istwfk = %s not implemented" % self.istwfk)

        if s0 == 1:
            # Reinstate input shape
            arr_on_mesh.shape = mesh.shape

        return arr_on_mesh

    def fromfftmesh(self, mesh, arr_on_mesh):
        """
        Transfer ``arr_on_mesh`` given on the FFT mesh to the G-sphere.
        """
        indim =  arr_on_mesh.ndim
        arr_on_mesh = mesh.reshape(arr_on_mesh)
        ishape = arr_on_mesh.shape
        s0 = ishape[0]

        arr_on_sphere = np.empty((s0,) + (self.npw,), dtype=arr_on_mesh.dtype)

        if self.istwfk == 1:
            #do ig=1,npwout
            #  i1=kg_kout(1,ig); if (i1<0) i1=i1+n1; i1=i1+1
            #  i2=kg_kout(2,ig); if (i2<0) i2=i2+n2; i2=i2+1
            #  i3=kg_kout(3,ig); if (i3<0) i3=i3+n3; i3=i3+1
            #end do
            n1, n2, n3 = mesh.shape

            for sph_idx, gvec in enumerate(self.gvecs):
                i1 = gvec[0]
                if i1 < 0: i1 = i1 + n1
                i2 = gvec[1]
                if i2 < 0: i2 = i2 + n2
                i3 = gvec[2]
                if i3 < 0: i3 = i3 + n3
                arr_on_sphere[..., sph_idx] = arr_on_mesh[..., i1, i2, i3]

        else:
            raise NotImplementedError("istwfk %s is not implemented" % self.istwfk)

        if s0 == 1 and indim == 1:
            # Reinstate input shape
            arr_on_sphere.shape = self.npw

        return arr_on_sphere

    #def rotate(self, symmop):
    #    """
    #    Returns a new `GSphere` centered on Sk.

    #    Args:
    #        symmop: Symmetry operation object.
    #    """
    #    # The problem in this approach is that G-spheres centered on the
    #    # same k-point might have G-vectors ordered in a different way
    #    # and therefore one cannot operate on two wavefunctions in reciprocal space
    #    # on the G-sphere without checking first that gvecs1 == gvecs2.
    #    # The best solution is to compute the list of g-vectors with a deterministic
    #    # algorithm, similar to the one used in Abinit and then create tables
    #    # defining the mapping btw the two sets
    #    if self.istwfk != 1:
    #        raise ValueError("istwfk %d not coded" % self.istwfk)

    #    # Rotate the k-point and the G-vectors
    #    rot_kpt = symmop.rotate_k(self.kpoint.frac_coords, wrap_tows=False)
    #    rot_gvecs = symmop.rotate_gvecs(self.gvecs)

    #    #rot_istwfk = istwfk(rot_kpt)
    #    rot_istwfk = self.istwfk

    #    new = self.__class__(self.ecut, self.lattice, rot_kpt, rot_gvecs, istwfk=rot_istwfk)
    #    return new


#def kpg_sphere(lattice, kcoords, ecut):
#    """
#    Set up the list of G vectors inside a sphere out to $ (1/2)*(2*\pi*(k+G))^2=ecut $
#    """
#    # Set up standard search sequence for grid points, in standard storage mode i.e.
#    # 0 1 2 3 ... g_max g_min ... -1
#    from pymatgen.core.units import Energy
#    ecut = Energy(ecut, "Ha").to("eV")
#    two_ecut = 2 * ecut
#
#    g1d_list = 3 * [None]
#
#    def norm2(vec):
#        return lattice.dot(vec, vec)
#
#    for dim in range(3):
#        rec_vec = lattice.matrix[dim,:]
#
#        # Compute g_max and g_min for this direction.
#        import itertools
#        for ig in itertools.count(start=0, step=1):
#            kpg = kcoords + (ig * rec_vec)
#            kpg2 = norm2(kpg)
#            if kpg2 > two_ecut:
#                g_max = ig
#                break
#
#        for ig in itertools.count(start=-1, step=-1):
#            kpg = kcoords + (ig * rec_vec)
#            kpg2 = norm2(kpg)
#            if kpg2 > two_ecut:
#                g_min = ig
#                break
#
#        g1d_list[dim] = list(range(g_max)) + list(range(-1, g_min, -1))
#
#    gx_list = g1d_list[0]
#    gy_list = g1d_list[1]
#    gz_list = g1d_list[2]
#
#    # Compute the list of G-vectors. Note that the Gs are ordered
#    # according to the Fortran convention.
#    gvecs = []
#    app = gvecs.append
#
#    for gvec in itertools.product(gz_list, gy_list, gx_list):
#        gvec = np.array(gvec)
#        kpg2 = norm2(kcoords + gvec)
#        if kpg2 <= two_ecut:
#            app(gvec)
#
#    return np.array(gvecs, dtype=np.int)
