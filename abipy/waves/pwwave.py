# coding: utf-8
"""This module contains the class describing a planewave wavefunction."""
from __future__ import print_function, division, unicode_literals, absolute_import

import tempfile
import copy
import six
import itertools
import numpy as np

from abipy.iotools import Visualizer
from abipy.iotools.xsf import xsf_write_structure, xsf_write_data
from abipy.core import Mesh3D
from abipy.core.kpoints import Kpoint

__all__ = [
    "PWWaveFunction",
]


class WaveFunction(object):
    """
    Abstract class defining base and abstract methods for wavefunction objects.
    """
    def __eq__(self, other):
        if other is None: return False
        if self.gsphere != other.gsphere: return False
        return np.allclose(self.ug, other.ug)

    def __ne__(self, other):
        return not (self == other)

    def __iter__(self):
        """Yields G, ug[0:nspinor, G]"""
        if six.PY2:
            return itertools.izip(self.gvecs, self.ug.T)
        else:
            return zip(self.gvecs, self.ug.T)

    def __getitem__(self, slice):
        return self.gvecs[slice], self.ug[:, slice]

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.to_string()

    @property
    def shape(self):
        """Shape of ug i.e. (nspinor, npw)"""
        return self.nspinor, self.npw

    @property
    def gsphere(self):
        """:class:`GSphere` object"""
        return self._gsphere

    @property
    def kpoint(self):
        """:class:`Kpoint` object"""
        return self.gsphere.kpoint

    @property
    def gvecs(self):
        """G-vectors in reduced coordinates."""
        return self.gsphere.gvecs

    @property
    def npw(self):
        """Number of G-vectors."""
        return len(self.gsphere)

    @property
    def ecut(self):
        """Cutoff energy in Hartree."""
        return self.gsphere.ecut

    @property
    def isnc(self):
        """True if norm-conserving wavefunction."""
        return isinstance(self, PWWaveFunction)

    @property
    def ispaw(self):
        """True if PAW wavefunction."""
        return isinstance(self, PAW_WaveFunction)

    @property
    def ug(self):
        """Periodic part of the wavefunctions in G-space."""
        return self._ug

    def set_ug(self, ug):
        """Set the value of the u(nspinor, G) array."""
        assert ug.shape == self.shape
        self._ug = ug
        self.delete_ur()

    @property
    def ur(self):
        """Periodic part of the wavefunction in real space."""
        try:
            return self._ur
        except AttributeError:
            self._ur = self.fft_ug()
            return self._ur

    def delete_ur(self):
        """Delete _u(r) (if it has been computed)."""
        try:
            del self._ur
        except AttributeError:
            pass

    @property
    def ur_xyz(self):
        """
        Returns a copy with ur[nspinor, nx, ny, nz]. Mainly used for post-processing.
        """
        return self.mesh.reshape(self.ur).copy()

    @property
    def ur2_xyz(self):
        """
        Returns ur2[nx, ny, nz]. Mainly used for post-processing.
        """
        return self.mesh.reshape(self.ur2)

    @property
    def mesh(self):
        """The mesh used for the FFT."""
        return self._mesh

    def set_mesh(self, mesh):
        """Set the FFT mesh. :math:`u(r)` is computed on this box."""
        assert isinstance(mesh, Mesh3D)
        self._mesh = mesh
        self.delete_ur()

    def deepcopy(self):
        """Deep copy of self."""
        return copy.deepcopy(self)

    def ug_mesh(self, mesh=None):
        """
        Returns u(G) on the FFT mesh,

        Args:
            mesh: :class:`Mesh3d` object. If mesh is None, self.mesh is used.
        """
        mesh = self.mesh if mesh is None else mesh
        ug_mesh = self.gsphere.tofftmesh(mesh, self.ug)
        return ug_mesh

    def fft_ug(self, mesh=None):
        """
        Performs the FFT transform of :math:`u(g)` on mesh.

        Args:
            mesh: :class:`Mesh3d` object. If mesh is None, self.mesh is used.

        Returns:
            :math:`u(r)` on the real space FFT box.
        """
        mesh = self.mesh if mesh is None else mesh
        ug_mesh = self.ug_mesh(mesh)
        return mesh.fft_g2r(ug_mesh, fg_ishifted=False)

    def to_string(self, prtvol=0):
        """String representation."""
        lines = []
        app = lines.append
        app("%s: nspinor = %d, spin = %d, band = %d " % (
            self.__class__.__name__, self.nspinor, self.spin, self.band))

        if hasattr(self, "gsphere"):
            app(self.gsphere.tostring(prtvol))

        if hasattr(self, "mesh"):
            app(self.mesh.to_string(prtvol))

        return "\n".join(lines)

    @property
    def ur2(self):
        """Return :math:`||u(r)||^2` in real space."""
        ur2 = self.ur.conj() * self.ur
        #ur2 = self.ur[0].conj() * self.ur[0]
        #if self.nspinor == 2:
        #    ur2 += self.ur[1].conj() * self.ur[1]

        # copy to have contiguous data.
        return ur2.real.copy()


class PWWaveFunction(WaveFunction):
    """
    This object describes a wavefunction expressed in a plane-wave basis set.
    """
    def __init__(self, nspinor, spin, band, gsphere, ug):
        """
        Creation method.

        Args:
            nspinor: number of spinorial components.
            spin: spin index.
            band: band index (>=0)
            gsphere :class:`GSphere` instance.
            ug: 2D array containing u[nspinor,G] for G in gsphere.
        """
        self.nspinor, self.spin, self.band = nspinor, spin, band
        # Sanity check.
        assert ug.ndim == 2
        assert ug.shape[0] == nspinor
        assert ug.shape[1] == gsphere.npw

        self._gsphere = gsphere
        self._ug = np.array(ug)

    def norm2(self, space="g"):
        """Return :math:`||\psi||^2` computed in G- or r-space."""
        space = space.lower()

        if space == "g":
            return np.real(np.vdot(self.ug, self.ug))

        elif space == "r":
            return np.real(self.mesh.integrate(self.ur2))

        else:
            raise ValueError("Wrong space: %s" % space)

    def export_ur2(self, filename, structure, visu=None):
        """
        Export u(r)**2 on file filename.

        Args:
            filename: String specifying the file path and the file format.
                The format is defined by the file extension. filename="prefix.xsf", for example,
                will produce a file in XSF format. An *empty* prefix, e.g. ".xsf" makes the code use a temporary file.
            structure: :class:`Structure` object.
            visu: :class:`Visualizer` subclass. By default, this method returns the first available
                visualizer that supports the given file format. If visu is not None, an
                instance of visu is returned. See :class:`Visualizer` for the list of
                applications and formats supported.

        Returns:
            Instance of :class:`Visualizer`
        """
        if "." not in filename:
            raise ValueError("Cannot detect file extension in: %s" % filename)

        tokens = filename.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]: # fname == ".ext" ==> Create temporary file.
            filename = tempfile.mkstemp(suffix="." + ext, text=True)[1]
            print("Creating temporary file: %s" % filename)

        # Compute |u(r)|2 and write data according to ext.
        ur2 = np.reshape(self.ur2, (1,) + self.ur2.shape)

        with open(filename, mode="w") as fh:
            if ext == "xsf":
                # xcrysden
                xsf_write_structure(fh, structures=[structure])
                xsf_write_data(fh, structure, ur2, add_replicas=True)
            else:
                raise NotImplementedError("extension %s is not supported." % ext)

        if visu is None:
            return Visualizer.from_file(filename)
        else:
            return visu(filename)

    def visualize_ur2(self, structure, visu_name):
        """
        Visualize u(r)**2 visualizer.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        # Get the Visualizer subclass from the string.
        visu = Visualizer.from_name(visu_name)

        # Try to export data to one of the formats supported by the visualizer
        # Use a temporary file (note "." + ext)
        for ext in visu.supported_extensions():
            ext = "." + ext
            try:
                return self.export_ur2(ext, structure, visu=visu)
            except visu.Error:
                pass
        else:
            raise visu.Error("Don't know how to export data for %s" % visu_name)

    #def tkin(self):
    #    """Computes the matrix element of the kinetic operator in reciprocal space."""
    #    tug = -0.5 * self.gsphere.kpg2 * self.ug
    #    return np.vdot(self.ug, tug).sum()

    def braket(self, other, space="g"):
        """
        Returns the scalar product <u1|u2> of the periodic part of two wavefunctions
        computed in G-space or r-space, depending on the value of space.

        Args:
            other: Other wave (right-hand side)
            space:  Integration space. Possible values ["g", "gsphere", "r"]
                if "g" or "r" the scalar product is computed in G- or R-space on the FFT box.
                if space="gsphere" the integration is done on the G-sphere. Note that
                this option assumes that self and other have the same list of G-vectors.
        """
        space = space.lower()

        if space == "g":
            ug1_mesh = self.gsphere.tofftmesh(self.mesh, self.ug)
            ug2_mesh = other.gsphere.tofftmesh(self.mesh, other.ug)
            return np.vdot(ug1_mesh, ug2_mesh)

        elif space == "gsphere":
            return np.vdot(self.ug, other.ug)

        elif space == "r":
            return np.vdot(self.ur, other.ur) * self.mesh.dv

        else:
            raise ValueError("Wrong space: %s" % space)

    #def pww_translation(self, gvector, rprimd):
    #    """Returns the pwwave of the kpoint translated by one gvector."""
    #    gsph = self.gsphere.copy()
    #    wpww = PWWaveFunction(self.nspinor, self.spin, self.band, gsph, self.ug.copy())
    #    wpww.mesh = self.mesh
    #    wpww.pww_translation_inplace(gvector, rprimd)
    #    return wpww

    #def pww_translation_inplace(self, gvector, rprimd):
    #    """Translates the pwwave from 1 kpoint by one gvector."""
    #    self.gsphere.kpoint = self.gsphere.kpoint + gvector
    #    self.gsphere.gvecs = self.gsphere.gvecs + gvector
    #    fft_ndivs = (self.mesh.shape[0] + 2, self.mesh.shape[1] + 2, self.mesh.shape[2] + 2)
    #    newmesh = Mesh3D(fft_ndivs, rprimd, pbc=True)
    #    self.mesh = newmesh

    #def pwwtows_inplace(self):
    #    """Wrap the kpoint to the interval ]-1/2,1/2] and update pwwave accordingly."""
    #    kpoint = Kpoint(self.gsphere.kpoint, self.gsphere.gprimd)
    #    wkpt = kpoint.wrap_to_ws()

    #    if np.allclose(wkpt.rcoord, kpoint.rcoord):
    #        return

    #    #@David FIXME this is wrong
    #    gvector = np.array(kpoint.rcoord - wkpt.rcoord, np.int)
    #    self.gsphere.gvecs = self.gsphere.gvecs + gvector
    #    self.gsphere.kpoint = wkpt.rcoord

    #def pwwtows(self):
    #    """Return the pwwave of the kpoint wrapped to the interval ]-1/2,1/2]."""
    #    gsph = self.gsphere.copy()
    #    wpww = PWWaveFunction(self.nspinor, self.spin, self.band, gsph, self.ug.copy())
    #    wpww.pwwtows_inplace()
    #    wpww.mesh = self.mesh
    #    return wpww

    def rotate(self, symmop, mesh=None):
        """
        Rotate the pwwave by the symmetry operation symmop.

        Args:
            symmop: :class:`Symmetry` operation
            mesh: mesh for the FFT, if None the mesh of self is used.

        Returns:
            New wavefunction object.
        """
        if self.nspinor != 1:
            raise ValueError("Spinor rotation not available yet.")

        rot_gsphere = self.gsphere.rotate(symmop)
        #rot_istwfk = istwfk(rot_kpt)

        if not np.allclose(symmop.tau, np.zeros(3)):
            rot_ug = np.empty_like(self.ug)
            rot_gvecs = rot_gsphere.gvecs
            rot_kpt = rot_gsphere.kpoint.frac_coords

            ug = self._ug
            #phase = np.exp(-2j * np.pi * (np.dot(rot_gvecs + rot_kpt, symmop.tau)))
            for ig in range(self.npw):
                rot_ug[:, ig] = ug[:, ig] * np.exp(-2j * np.pi * (np.dot(rot_gvecs[ig] + rot_kpt, symmop.tau)))
        else:
            rot_ug = self.ug.copy()

        # Invert the collinear spin if we have an AFM operation
        rot_spin = self.spin
        if self.nspinor == 1:
            rot_spin = self.spin if symmop.is_fm else (self.spin + 1) % 2

        # Build new wave and set the mesh.
        new = self.__class__(self.nspinor, rot_spin, self.band, rot_gsphere, rot_ug)
        new.set_mesh(mesh if mesh is not None else self.mesh)
        return new


class PAW_WaveFunction(WaveFunction):
    """
    All the methods that are related to the all-electron representation should start with ae.
    """
