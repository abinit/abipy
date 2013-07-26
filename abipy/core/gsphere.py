"""This module contains the class defining the G-sphere for wavefunctions, densities and potentials"""
from __future__ import print_function, division

import collections
import numpy as np

from abipy.core.kpoints import Kpoint

__all__ = [
    "GSphere",
]


class GSphere(collections.Sequence):
    """Descriptor-class for the G-sphere."""

    def __init__(self, ecut, lattice, kpoint, gvecs, istwfk=1):
        """
        Args:
            ecut:
                Cutoff energy in Hartree.
            lattice:
                Reciprocal lattice.
            kpoint:
                Reduced coordinates of the k-point.
            gvecs:
                Array with the reduced coordinates of the G-vectors.
            istwfk:
                Storage option (time-reversal symmetry, see abinit variable)
        """
        self.ecut = ecut
        self.lattice = lattice
        self.kpoint = Kpoint.askpoint(kpoint, lattice)
        self._gvecs = np.reshape(gvecs, (-1, 3))
        self.npw = self.gvecs.shape[0]

        self.istwfk = istwfk

        if istwfk != 1:
            raise NotImplementedError("istwfk %d is not implemented" % self.istwfk)

        # min and Max reduced component for each direction.
        #self.gmin = gvecs.min(axis=0)
        #self.gmax = gvecs.max(axis=0)

    @property
    def gvecs(self):
        """G-vectors in reduced coordinates."""
        return self._gvecs

    # Sequence protocol
    def __len__(self):
        return self.gvecs.shape[0]

    def __getitem__(self, slice):
        return self.gvecs[slice]

    def __iter__(self):
        return self.gvecs.__iter__()

    def __contains__(self, gvec):
        return np.asarray(gvec) in self

    def index(self, gvec):
        """
        return the index of the G-vector gvec in self.
        Raises ValueError if the value is not present.
        """
        for i, g in enumerate(self):
            if np.all(g == gvec):
                return i
        else:
            raise ValueError("Cannot find %s in Gsphere" % gvec)

    def count(self, gvec):
        """Return number of occurrences of gvec."""
        return np.count_nonzero(np.all(g == gvec) for g in self)

    def __str__(self):
        return self.tostring()

    def __eq__(self, other):
        if other is None: return False
        return (self.ecut == other.ecut and
                np.all(self.lattice == other.lattice) and
                self.kpoint == other.kpoint and
                np.all(self.gvecs == other.gvecs) and
                self.istwfk == other.istwfk
                )

    def __ne__(self, other):
        return not self == other

    def copy(self):
        """Deep copy."""
        return GSphere(self.ecut, self.lattice.copy(), self.kpoint.copy(), self.gvecs.copy(), istwfk=self.istwfk)

    def tostring(self, prtvol=0):
        """String representation."""
        s = "GSphere: kpoint = %(kpoint)s, ecut = %(ecut)f, npw = %(npw)d, istwfk = %(istwfk)d" % (
            self.__dict__)
        return s

    def _new_array(self, dtype=np.float, zero=True, extra_dims=()):
        """Returns a numpy array defined on the sphere."""
        shape = (self.npw,)

        if isinstance(extra_dims, int): extra_dims = (extra_dims,)
        shape = extra_dims + tuple(shape)

        if zero:
            return np.zeros(shape, dtype)
        else:
            return np.empty(shape, dtype)

    def zeros(self, dtype=np.float, extra_dims=()):
        """
        Returns new zeroed 1D array.

        The type can be set with the ``dtype`` keyword.
        Extra dimensions can be added with ``extra_dims``.
        """
        return self._new_array(dtype=dtype, zero=True, extra_dims=extra_dims)

    def czeros(self, extra_dims=()):
        """New zeroed 1D complex array."""
        return self._new_array(dtype=np.complex, zero=True, extra_dims=extra_dims)

    def empty(self, dtype=np.float, extra_dims=()):
        """
        Returns new uninitialized 1D array.

        The type can be set with the ``dtype`` keyword.
        Extra dimensions can be added with ``extra_dims``.
        """
        return self._new_array(dtype=dtype, zero=False, extra_dims=extra_dims)

    def cempty(self, extra_dims=()):
        """Returns New uninitialized 1D complex array."""
        return self._new_array(dtype=np.complex, zero=False, extra_dims=extra_dims)

    #def build_fftbox(self, boxsph_ratio=1.05):
    #  """Returns the number of divisions of the FFT box enclosing the sphere."""
    #  #return ndivs

    def tofftmesh(self, mesh, arr_on_sphere):
        """Insert arr_on_sphere in the FFT mesh."""
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
             #  do idat=1,ndat
             #    cfft(1,i1,i2,i3+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
             #    cfft(2,i1,i2,i3+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
             #  end do
             #end do

            (n1, n2, n3) = mesh.shape
            for sph_idx, gvec in enumerate(self.gvecs):
                i1 = gvec[0]
                if i1 < 0: i1 = i1 + n1
                i2 = gvec[1]
                if i2 < 0: i2 = i2 + n2
                i3 = gvec[2]
                if i3 < 0: i3 = i3 + n3
                arr_on_mesh[...,i1,i2,i3] = arr_on_sphere[...,sph_idx]
        else:
            raise NotImplementedError("istwfk = %s not implemented" % self.istwfk)

        if s0 == 1:  # Reinstate input shape
            arr_on_mesh.shape = mesh.shape

        return arr_on_mesh

    def fromfftmesh(self, mesh, arr_on_mesh):
        """Transfer arr_on_mesh given on the FFT mesh to the G-sphere."""
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
           #  do idat=1,ndat
           #    fofgout(1,ig+npwout*(idat-1))=fofr(1,i1,i2,i3+n3*(idat-1))
           #    fofgout(2,ig+npwout*(idat-1))=fofr(2,i1,i2,i3+n3*(idat-1))
           #  end do
           #end do
            (n1, n2, n3) = mesh.shape

            for sph_idx, gvec in enumerate(self.gvecs):
                i1 = gvec[0]
                if i1 < 0: i1 = i1 + n1
                i2 = gvec[1]
                if i2 < 0: i2 = i2 + n2
                i3 = gvec[2]
                if i3 < 0: i3 = i3 + n3
                arr_on_sphere[...,sph_idx] = arr_on_mesh[...,i1,i2,i3]
        else:
            raise NotImplementedError("istwfk = " + str(self.istwfk) + " is not implemented")

        if s0 == 1 and indim == 1:  # Reinstate input shape
            arr_on_sphere.shape = self.npw

        return arr_on_sphere

#########################################################################################

#def kpgsphere(kpoint, ecut, gmet, istwfk)
# """Set up the list of G vectors inside a sphere out to $ (1/2)*(2*\pi*(k+G))^2=ecut $"""
# return gvec
