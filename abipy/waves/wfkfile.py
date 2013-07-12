"""Wavefunction file."""
from __future__ import print_function, division

import os
import numpy as np

from abipy.core import Mesh3D, GSphere, Structure
from abipy.iotools import ETSF_Reader, Visualizer
from abipy.electrons import ElectronBands
from abipy.kpoints import kpoints_factory
from abipy.waves.pwwave import PWWaveFunction

__all__ = [
    "WFK_File",
]


class WFK_File(object):
    """
    This object provides simple interfaces to access and analyze
    the data stored in the WFK file produced by ABINIT.
    """
    def __init__(self, path):
        """
        Initialized the object from a Netcdf file.
        """
        self.path = os.path.abspath(path)

        # Initialize the  structure from file.
        self.structure = Structure.from_ncfile(path)

        # Initialize the band energies.
        self.bands = ElectronBands.from_ncfile(path)

        self.kpoints = kpoints_factory(path)
        self.nkpt = len(self.kpoints)

        with WFK_Reader(self.path) as reader:
            assert reader.has_pwbasis_set
            assert reader.cplex_ug == 2
            print(reader)
            self.npwarr = reader.npwarr
            self.nband_sk = reader.nband_sk

            self.nspinor = reader.nspinor
            self.nsppol = reader.nsppol
            self.nspden = reader.nspden

        # FFT mesh (augmented divisions reported in the WFK file)
        self.fft_mesh = Mesh3D(reader.fft_divs, self.structure.lattice_vectors())

        # Build G-spheres for each k-point
        gspheres = len(self.kpoints) * [None]
        ecut = reader.ecut
        for k, kpoint in enumerate(self.kpoints):
            gvec_k, istwfk = reader.read_gvec_istwfk(k)
            gspheres[k] = GSphere(ecut, self.structure.reciprocal_lattice, kpoint, gvec_k, istwfk=istwfk)

        self._gspheres = tuple(gspheres)

        # Save reference to the reader.
        self.reader = reader

    @classmethod
    def from_ncfile(cls, path):
        """Initialize the object from a Netcdf file."""
        return cls(path)

    @property
    def basename(self):
        """Basename of the WFK file"""
        return os.path.basename(self.path)

    @property
    def gspheres(self):
        """List of `GSphere` ordered by k-points."""
        return self._gspheres

    @property
    def mband(self):
        """Maximum band index"""
        return np.max(self.nband_sk)

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        """
        String representation

        Args:
            prtvol:
                verbosity level.
        """
        keys = ["nspinor", "nspden"]
        lines = []
        app = lines.append
        for k in keys:
            try:
                value = self.__dict__[k]
                if prtvol == 0 and isinstance(value, np.ndarray):
                    continue
                app("%s = %s" % (k, value))
            except KeyError:
                pass

        return "\n".join(lines)

    def kindex(self, kpoint):
        """The index of the k-point in the file. Accepts: `Kpoint` object or int."""
        return self.reader.kindex(kpoint)

    def get_structure(self):
        """Returns the `Structure`"""
        return self.structure

    def get_bands(self):
        """Return an instance of `ElectronBands`"""
        return self.bands

    def get_wave(self, spin, kpoint, band):
        """
        Read and return the wavefunction with the given spin, band and kpoint.

        Args:
            spin:
                spin index (0,1)
            kpoint:
                Either `Kpoint` instance or integer giving the sequential index in the IBZ (C-convention).
            band:
                band index.

            returns:
                `WaveFunction` instance.
        """
        k = self.kindex(kpoint)

        if (spin not in range(self.nsppol) or k not in range(self.nkpt) or
            band not in range(self.nband_sk[spin,k])):
            raise ValueError("Wrong (spin, band, kpt) indices")

        ug_skb = self.reader.read_ug(spin, kpoint, band)

        # Istanciate the wavefunction object and set the FFT mesh
        # using the divisions reported in the WFK file.
        wave = PWWaveFunction(self.nspinor, spin, band, self.gspheres[k], ug_skb)
        wave.set_mesh(self.fft_mesh)

        return wave

    def plot_bands(self, klabels=None, *args, **kwargs):
        """
        Plot the energy bands. See the :func:`ElectronBands.plot` for the meaning of the variables""
        """
        return self.bands.plot(klabels=klabels, *args, **kwargs)

    def export_structure(self, path):
        """
        Export the structure on file.

        returns:
            Instance of :class:`Visualizer`
        """
        return self.structure.export(path)

    def visualize_structure_with(self, visualizer):
        """
        Visualize the crystalline structure with visualizer.

        See :class:`Visualizer` for the list of applications and formats supported.
        """
        extensions = Visualizer.exts_from_appname(visualizer)
                                                                                                 
        for ext in extensions:
            ext = "." + ext
            try:
                return self.export_structure(ext)
            except Visualizer.Error as exc:
                #print(exc)
                pass
        else:
            raise Visualizer.Error(
                "Don't know how to export data for visualizer %s" % visualizer)

    def export_ur2(self, path, spin, kpoint, band):
        """
        Export :math:`|u(r)|^2` on file filename.

        returns:
            Instance of :class:`Visualizer`
        """
        # Read the wavefunction from file.
        wave = self.get_wave(spin, kpoint, band)

        # Export data uding the format specified by filename.
        return wave.export_ur2(path, self.structure)

    #def visualize_ur2_with(self, spin, kpoint, bands, visualizer):
    #    """
    #    Visualize :math:`|u(r)|^2`  with visualizer.

    #    See :class:`Visualizer` for the list of applications and formats supported.
    #    """
    #    extensions = Visualizer.exts_from_appname(visualizer)
    #
    #    for ext in extensions:
    #        ext = "." + ext
    #        try:
    #            return self.export_ur2(ext)
    #        except Visualizer.Error:
    #            pass
    #    else:
    #        raise Visualizer.Error(
    #            "Don't know how to export data for visualizer %s" % visualizer)

#########################################################################################


class WFK_Reader(ETSF_Reader):
    """This object reads data from the WFK file."""
    def __init__(self, path):
        """Initialize the object from a filename."""
        super(WFK_Reader, self).__init__(path)

        self.structure = self.read_structure()
        self.kpoints = kpoints_factory(path)

        self.nfft1 = self.read_dimvalue("number_of_grid_points_vector1")
        self.nfft2 = self.read_dimvalue("number_of_grid_points_vector2")
        self.nfft3 = self.read_dimvalue("number_of_grid_points_vector3")

        self.cplex_ug = self.read_dimvalue("real_or_complex_coefficients")

        self.nspinor = self.read_dimvalue("number_of_spinor_components")
        self.nsppol = self.read_dimvalue("number_of_spins")
        self.nspden = self.read_dimvalue("number_of_components")

        self.ecut = self.read_value("kinetic_energy_cutoff")
        self.nband_sk = self.read_value("number_of_states")
        self.istwfk = self.read_value("istwfk")
        self.npwarr = self.read_value("number_of_coefficients")

        # Gvectors
        self._kg = self.read_value("reduced_coordinates_of_plane_waves")

        # Wavefunctions (complex array)
        # TODO use variables to avoid storing the full block.
        if self.cplex_ug == 2:
            self.set_of_ug = self.read_value("coefficients_of_wavefunctions", cmode="c")
        else:
            raise NotImplementedError("")

    @property
    def basis_set(self):
        """String defining the basis set."""
        try:
            return self._basis_set
        except AttributeError:
            basis_set = self.read_value("basis_set")
            self._basis_set = "".join([c for c in basis_set]).strip()
            return self._basis_set

    @property
    def has_pwbasis_set(self):
        """True if the plane-wave basis set is used."""
        return self.basis_set == "plane_waves"

    @property
    def fft_divs(self):
        """FFT divisions used to compute the data in the WFK file."""
        return self.nfft1, self.nfft2, self.nfft3

    def kindex(self, kpoint):
        """
        Index of the k-point in the internal tables.

        Accepts: `Kpoint` instance of integer.
        """
        if isinstance(kpoint, int):
            return kpoint
        else:
            try:
                return self.kpoints.index(kpoint)
            except:
                raise

    def read_kpoints(self):
        return self.kpoints

    def read_gvec_istwfk(self, kpoint):
        """
        Read the set of G-vectors and the value of istwfk for the given k-point.
        Accepts `Kpoint` object or integer.
        """
        k = self.kindex(kpoint)
        npw_k, istwfk = self.npwarr[k], self.istwfk[k]
        return self._kg[k,:npw_k,:], istwfk

    def read_ug(self, spin, kpoint, band):
        """Read the Fourier components of the wavefunction."""
        k = self.kindex(kpoint)
        npw_k, istwfk = self.npwarr[k], self.istwfk[k]
        # TODO use variables to avoid storing the full block.
        return self.set_of_ug[spin,k,band,:,:npw_k]

