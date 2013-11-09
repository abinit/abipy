"""This module contains the class describing densities in real space on uniform 3D meshes."""
from __future__ import print_function, division

import numpy as np
import tempfile

from abipy.tools import transpose_last3dims
from abipy.iotools import Visualizer, xsf
from .mesh3d import Mesh3D

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
                numpy array with the scalar field in real space.
            structure:
                `Structure` object describing the crystalline structure.
            iorder:
                Order of the array. "c" for C ordering, "f" for Fortran ordering.
        """
        self.nspinor, self.nsppol, self.nspden = nspinor, nsppol, nspden
        self.datar = datar
        self.structure = structure

        iorder = iorder.lower()
        assert iorder in ["f", "c"]

        if iorder == "f": # (z,x,y) --> (x,y,z)
            self.datar = transpose_last3dims(self.datar)

        # Init Mesh3D
        mesh_shape = self.datar.shape[-3:]
        self.mesh = Mesh3D(mesh_shape, structure.lattice_vectors())

        # Make sure we have the correct shape.
        self.datar = np.reshape(self.datar, (nspden,) + self.mesh.shape)

        # FFT R --> G.
        self.datag = self.mesh.fft_r2g(self.datar)
        #print self.datag[...,0,0,0] * structure.volume / np.product(datar.shape[-3:])

    def __len__(self):
        return len(self.datar)

    def __getitem__(self, slice):
        return self.datar[slice], self.datag[slice]

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        """String representation"""

        s  = "ScalarField: nspinor = %i, nsppol = %i, nspden = %i" % (
            self.nspinor, self.nsppol, self.nspden)
        s += "  " + self.mesh.tostring(prtvol)
        if prtvol > 0:
            s += "  " + str(self.structure)

        return s

    @property
    def shape(self):
        shape_r, shape_g = self.datar.shape, self.datag.shape
        assert np.all(shape_r == shape_g)

        return shape_r

    @property
    def iscollinear(self):
        """True if collinear i.e. nspinor==1."""
        return self.nspinor == 1

    def _set_datar(self, datar):
        """Set the value of datar (datag is updated accordingly)."""
        self.datar = datar
        self.datag = self.mesh.fft_r2g(datar)

    def _set_datag(self, datag):
        """Set the value of datag (datar is updated accordingly)."""
        self.datag = datag
        self.datar = self.mesh.fft_g2r(datag)

    #def bracket(self, bra_wave, ket_wave):
    #    """
    #    Compute the matrix element of the local operator v(r) in real space
    #    """
    #    if bra_wave.mesh != self.mesh:
    #       bra_ur = bra_wave.fft_ug(self.mesh)
    #    else:
    #       bra_ur = bra_wave.ur

    #    if ket_wave.mesh != self.mesh:
    #       ket_ur = ket_wave.fft_ug(self.mesh)
    #    else:
    #       ket_ur = ket_wave.ur

    #    return self.mesh.integrate(bra_ur.conj() * self.datar * ur_right)

    #def interpolate(self, points, method="linear", space="r")

    #def fourier_interp(self, new_mesh):
    #  intp_datar =
    #  return DFTScalarField(self.nspinor, self.nsppol, self.nspden, self.structure, intp_datar)

    def export(self, filename):
        """
        Export the DftScalarField on file filename. Format is defined by the extension in filename.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        if "." not in filename:
            raise ValueError(" Cannot detect file extension in filename: %s " % filename)

        tokens = filename.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]: # filename == ".ext" ==> Create temporary file.
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
