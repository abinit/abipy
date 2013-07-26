"""This module contains the class describing a planewave wavefunction."""
from __future__ import print_function, division

import numpy as np
import tempfile

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

    def __mul__(self, other):
        raise NotImplementedError("Dot product not yet implemented for PWWave objects")

    #def __iter__(self):
    #def __getitem__(self, slice):

    @property
    def shape(self):
        return self.nspinor, self.npw

    @property
    def gsphere(self):
        """`GSphere` object"""
        return self._gsphere

    @property
    def kpoint(self):
        """`Kpoint` object"""
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
    def ug(self):
        """Periodic part of the wavefunctions in G-space."""
        return self._ug

    def set_ug(self, ug):
        assert ug.shape == self.shape
        set._ug = ug
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

    def set_mesh(self, mesh):
        """Set the FFT mesh. :math:`u(r)` is computed on this box."""
        assert isinstance(mesh, Mesh3D)
        self.mesh = mesh
        self.delete_ur()

    def fft_ug(self):
        """
        Performs the FFT transform of :math:`u(g)`.
        Returns :math:`u(r)` on the real space FFT box.
        """
        ug_mesh = self.gsphere.tofftmesh(self.mesh, self.ug)
        return self.mesh.fft_g2r(ug_mesh, fg_ishifted=False)

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        """String representation."""
        lines = []
        app = lines.append
        app("%s: nspinor = %d, spin = %d, band = %d " % (
            self.__class__.__name__, self.nspinor, self.spin, self.band))

        if hasattr(self, "gsphere"):
            app(self.gsphere.tostring(prtvol))
        if hasattr(self, "mesh"):
            app(self.mesh.tostring(prtvol))

        return "\n".join(lines)

    @property
    def ur2(self):
        """Return :math:`||u(r)||^2` in real space."""
        ur2 = self.ur.conj() * self.ur
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
            nspinor:
                number of spinorial components.
            spin:
                spin index.
            band:
                band index (>=0)
            gsphere
                Gsphere instance.
            ug:
                2D array containing u(nspinor,G) for G in gsphere.
        """
        self.nspinor, self.spin, self.band = nspinor, spin, band
        # Sanity check.
        assert ug.ndim == 2
        assert ug.shape[0] == nspinor
        assert ug.shape[1] == gsphere.npw

        self._gsphere = gsphere
        self._ug = ug

    def norm2(self, space="g"):
        """Return :math:`||\psi||^2` computed in G- or r-space."""
        if space.lower() == "g":
            return np.real(np.vdot(self.ug, self.ug))
        elif space.lower() == "r":
            return np.real(self.mesh.integrate(self.ur2))
        else:
            raise ValueError("Wrong space: %s" % space)

    def export_ur2(self, filename, structure):
        """
        Export the wavefunction on file filename.
        Format is defined by the extension in filename.
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

        return Visualizer.from_file(filename)

    def braket(self, wav2, space="g"):
        """Return <u|u2> computed in G- or r-space."""
        if space.lower() == "g":
            ug_mesh = self.gsphere.tofftmesh(self.mesh, self.ug)
            wav2_ug_mesh = wav2.gsphere.tofftmesh(self.mesh, wav2.ug)
            return np.vdot(ug_mesh, wav2_ug_mesh)
        elif space.lower() == "r":
            return np.vdot(self.ur, wav2.ur) / self.mesh.size
        else:
            raise ValueError("Wrong space: %s" % space)

    def pww_translation(self, gvector, rprimd):
        """Returns the pwwave of the kpoint translated by one gvector."""
        gsph = self.gsphere.copy()
        wpww = PWWaveFunction(self.nspinor, self.spin, self.band, gsph, self.ug.copy())
        wpww.mesh = self.mesh
        wpww.pww_translation_inplace(gvector, rprimd)
        return wpww

    def pww_translation_inplace(self, gvector, rprimd):
        """Translates the pwwave from 1 kpoint by one gvector."""
        self.gsphere.kpoint = self.gsphere.kpoint + gvector
        self.gsphere.gvecs = self.gsphere.gvecs + gvector
        fft_ndivs = (self.mesh.shape[0] + 2, self.mesh.shape[1] + 2, self.mesh.shape[2] + 2)
        newmesh = Mesh3D(fft_ndivs, rprimd, pbc=True)
        self.mesh = newmesh

    def pwwtows_inplace(self):
        """Wrap the kpoint to the interval ]-1/2,1/2] and update pwwave accordingly."""
        kpoint = Kpoint(self.gsphere.kpoint, self.gsphere.gprimd)
        wkpt = kpoint.wrap_to_ws()

        if np.allclose(wkpt.rcoord, kpoint.rcoord):
            return

        #@David FIXME this is wrong
        gvector = np.array(kpoint.rcoord - wkpt.rcoord, np.int)
        self.gsphere.gvecs = self.gsphere.gvecs + gvector
        self.gsphere.kpoint = wkpt.rcoord

    def pwwtows(self):
        """Return the pwwave of the kpoint wrapped to the interval ]-1/2,1/2]."""
        gsph = self.gsphere.copy()
        wpww = PWWaveFunction(self.nspinor, self.spin, self.band, gsph, self.ug.copy())
        wpww.pwwtows_inplace()
        wpww.mesh = self.mesh
        return wpww

    def rotate_inplace(self, symmop):
        """Rotate the pwwave by the symmetry operation symmop."""
        # CHANGE THIS !! Have a method that rotates the Ug's then calls the PWWaveFunction constructor!
        rkpt = symmop.rotate_k(self.kpoint, wrap_tows=False)
        rgvecs = symmop.rotate_g(self.gsphere.gvecs)

        if not np.allclose(symmop.tau, np.zeros(3)):
            rug = np.zeros_like(self.ug)
            for ii in range(np.shape(self.ug)[1]):
                rug[:, ii] = self.ug[:, ii] * np.exp(-2j * np.pi * (np.dot(rgvecs + rkpt, symmop.tau)))
            self._ug = rug
        self.gsphere.kpoint = rkpt
        self.gsphere.gvecs = rgvecs

        # Delete _u(r) if it has been computed.
        self.delete_ur()

    def rotate(self, symmop):
        """Return rotated pwwave by the symmetry operation symmop."""
        gsph = self.gsphere.copy()
        rpww = PWWaveFunction(self.nspinor, self.spin, self.band, gsph, self.ug.copy())
        rpww.rotate_inplace(symmop)
        rpww.mesh = self.mesh
        return rpww


class PAW_Wavefunction(PWWaveFunction):
    """
    A PAW wavefunction extends PWWavefunction adding new methods
    """
