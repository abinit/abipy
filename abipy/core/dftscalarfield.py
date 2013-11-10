"""This module contains the class describing densities in real space on uniform 3D meshes."""
from __future__ import print_function, division

import numpy as np

from abipy.tools import transpose_last3dims
from abipy.iotools import Visualizer, xsf
from abipy.core.mesh3d import Mesh3D

__all__ = [
    "DFTScalarField",
]


class DFTScalarField(object):

    def __init__(self, nspinor, nsppol, nspden, datar, structure, iorder="c"):
        """
        Args:
            nspinor:
                Number of spinorial components.
            nsppol:
                Number of spins.
            nspden:
                Number of spin density components.
            datar:
                numpy array with the scalar field in real space. shape [..., nx, ny, nz]
            structure:
                `Structure` object describing the crystalline structure.
            iorder:
                Order of the array. "c" for C ordering, "f" for Fortran ordering.
        """
        self.nspinor, self.nsppol, self.nspden = nspinor, nsppol, nspden
        self.structure = structure

        iorder = iorder.lower()
        assert iorder in ["f", "c"]

        if iorder == "f": # (z,x,y) --> (x,y,z)
            datar = transpose_last3dims(datar)

        # Init Mesh3D
        mesh_shape = datar.shape[-3:]
        self._mesh = Mesh3D(mesh_shape, structure.lattice_vectors())

        # Make sure we have the correct shape.
        self._datar = np.reshape(datar, (nspden,) + self.mesh.shape)

        # FFT R --> G.
        self._datag = self.mesh.fft_r2g(self.datar)
        #print self.datag[...,0,0,0] * structure.volume / np.product(datar.shape[-3:])

    def __len__(self):
        return len(self.datar)

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        """String representation"""

        s  = "%s: nspinor = %i, nsppol = %i, nspden = %i" % (
            self.__class__.__name__, self.nspinor, self.nsppol, self.nspden)
        s += "  " + self.mesh.tostring(prtvol)
        if prtvol > 0:
            s += "  " + str(self.structure)

        return s

    @property
    def datar(self):
        """`ndarrray` with data in real space."""
        return self._datar

    @property
    def datag(self):
        """`ndarrray` with data in reciprocal space."""
        return self._datag

    @property
    def mesh(self):
        """`Mesh3D`"""
        return self._mesh

    @property
    def shape(self):
        """Shape of the array."""
        shape_r, shape_g = self.datar.shape, self.datag.shape
        assert np.all(shape_r == shape_g)
        return shape_r

    @property
    def nx(self):
        """Number of points along x."""
        return self.mesh.nx

    @property
    def ny(self):
        """Number of points along y."""
        return self.mesh.ny

    @property
    def nz(self):
        """Number of points along z."""
        return self.mesh.nz

    @property
    def is_collinear(self):
        """True if collinear i.e. nspinor==1."""
        return self.nspinor == 1

    #@property
    #def datar_xyz(self):
    #    """
    #    Returns a copy with datar[nspden, nx, ny, nz]. 
    #    Mainly used for post-processing.
    #    """
    #    return self.mesh.reshape(self.datar).copy()

    #def braket_waves(self, bra_wave, ket_wave):
    #    """
    #    Compute the matrix element of the datar in real space
    #    """
    #    if bra_wave.mesh != self.mesh:
    #       bra_ur = bra_wave.fft_ug(self.mesh)
    #    else:
    #       bra_ur = bra_wave.ur

    #    if ket_wave.mesh != self.mesh:
    #       ket_ur = ket_wave.fft_ug(self.mesh)
    #    else:
    #       ket_ur = ket_wave.ur

    #    assert self.nspinor == 1
    #    assert bra_wave.spin == ket_wave.spin

    #    spin = bra_wave.spin
    #    datar = self.datar[spin]

    #    return self.mesh.integrate(bra_ur.conj() * datar * ket_ur)

    #def interpolate(self, points, method="linear", space="r")

    #def fourier_interp(self, new_mesh):
    #  intp_datar =
    #  return DFTScalarField(self.nspinor, self.nsppol, self.nspden, self.structure, intp_datar)

    def export(self, filename):
        """
        Export the real space data on file filename. 
        Format is defined by the extension in filename.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        if "." not in filename:
            raise ValueError(" Cannot detect file extension in filename: %s " % filename)

        tokens = filename.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]: # filename == ".ext" ==> Create temporary file.
            import tempfile
            filename = tempfile.mkstemp(suffix="."+ext, text=True)[1]

        with open(filename, mode="w") as fh:
            if ext == "xsf":
                # xcrysden
                xsf.xsf_write_structure(fh, self.structure)
                xsf.xsf_write_data(fh, self.structure, self.datar, add_replicas=True)
            else:
                raise NotImplementedError("extension %s is not supported." % ext)

        return Visualizer.from_file(filename)

    def visualize(self, visualizer):
        """
        Visualize data with visualizer.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        extensions = Visualizer.exts_from_appname(visualizer)
                                                                                                 
        for ext in extensions:
            ext = "." + ext
            try:
                return self.export(ext)
            except Visualizer.Error:
                pass
        else:
            raise Visualizer.Error("Don't know how to export data for visualizer %s" % visualizer)

    #def get_plane(self, plane, h):
    #    x, y, z = self.mesh.plane_inds(plane, h=h)
    #    plane = self.datar_xyz[:, x, y, z]
    #    new_shape = (plane.shape[0],) + tuple([s for s in plane.shape[-3:-1] if s > 1])
    #    return np.reshape(plane, new_shape)
