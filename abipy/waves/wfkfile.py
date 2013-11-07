"""Wavefunction file."""
from __future__ import print_function, division

import numpy as np

from abipy.core import Mesh3D, GSphere, Structure
from abipy.iotools import ETSF_Reader, Visualizer, AbinitNcFile, Has_Structure, Has_ElectronBands
from abipy.electrons import ElectronsReader
from abipy.waves.pwwave import PWWaveFunction

__all__ = [
    "WFK_File",
]

def iuptri(items, diago=True, with_inds=False):
    """
    Iterate over the upper triangle of the matrix (items x items)

    Args:
        items:
            Iterable object with elements [e0, e1, ...]
        diago:
            False if diagonal matrix elements should be excluded
        with_inds:
            If True, (i,j) (e_i, e_j) is returned else (e_i, e_j)

    >>> for (ij, mate) in iuptri([0,1], with_inds=True): 
    ...     print("ij:", ij, "mate:", mate)
    ij: (0, 0) mate: (0, 0)
    ij: (0, 1) mate: (0, 1)
    ij: (1, 1) mate: (1, 1)
    """
    for (ii, item1) in enumerate(items):
        for (jj, item2) in enumerate(items):
            do_yield = (jj >= ii) if diago else (jj > ii)
            if do_yield:
                if with_inds:
                    yield (ii, jj), (item1, item2)
                else:
                    yield item1, item2


def ilotri(items, diago=True, with_inds=False):
    """
    Iterate over the lower triangle of the matrix (items x items)

    Args:
        items:
            Iterable object with elements [e0, e1, ...]
        diago:
            False if diagonal matrix elements should be excluded
        with_inds:
            If True, (i,j) (e_i, e_j) is returned else (e_i, e_j)

    >>> for (ij, mate) in ilotri([0,1], with_inds=True): 
    ...     print("ij:", ij, "mate:", mate)
    ij: (0, 0) mate: (0, 0)
    ij: (1, 0) mate: (1, 0)
    ij: (1, 1) mate: (1, 1)
    """
    for (ii, item1) in enumerate(items):
        for (jj, item2) in enumerate(items):
            do_yield = (jj <= ii) if diago else (jj < ii)
            if do_yield:
                if with_inds:
                    yield (ii, jj), (item1, item2)
                else:
                    yield item1, item2


class WFK_File(AbinitNcFile, Has_Structure, Has_ElectronBands):
    """
    This object provides a simple interface to access and analyze
    the data stored in the WFK file produced by ABINIT.
    """
    def __init__(self, filepath):
        """
        Initialize the object from a Netcdf file.
        """
        super(WFK_File, self).__init__(filepath)

        with WFK_Reader(filepath) as reader:
            # Read the electron bands 
            self._ebands = reader.read_ebands()

            assert reader.has_pwbasis_set
            assert reader.cplex_ug == 2
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
            gvec_k, istwfk = reader.read_gvecs_istwfk(k)
            gspheres[k] = GSphere(ecut, self.structure.reciprocal_lattice, kpoint, gvec_k, istwfk=istwfk)

        self._gspheres = tuple(gspheres)

        # Save reference to the reader.
        self.reader = reader

    @property
    def structure(self):
        """`Structure` object"""
        return self.ebands.structure

    @property
    def ebands(self):
        """`ElectronBands` object"""
        return self._ebands

    @property
    def kpoints(self):
        return self.ebands.kpoints

    @property
    def nkpt(self):
        return len(self.kpoints)

    @property
    def gspheres(self):
        """List of `GSphere` objects ordered by k-points."""
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

        if (spin not in range(self.nsppol) or 
            k not in range(self.nkpt) or
            band not in range(self.nband_sk[spin, k])):
            raise ValueError("Wrong (spin, band, kpt) indices")

        ug_skb = self.reader.read_ug(spin, kpoint, band)

        # Istantiate the wavefunction object and set the FFT mesh
        # using the divisions reported in the WFK file.
        wave = PWWaveFunction(self.nspinor, spin, band, self.gspheres[k], ug_skb)
        wave.set_mesh(self.fft_mesh)

        return wave

    def export_ur2(self, filepath, spin, kpoint, band):
        """
        Export :math:`|u(r)|^2` on file filename.

        returns:
            Instance of :class:`Visualizer`
        """
        # Read the wavefunction from file.
        wave = self.get_wave(spin, kpoint, band)

        # Export data uding the format specified by filename.
        return wave.export_ur2(filepath, self.structure)

    def classify_ebands(self, spin, kpoint, bands_range, tol_ediff=1e-3):
        """
        Analyze the caracter of the bands at the given k-point and spin.

        Args:
            spin:
                Spin index.
            kpoint:
                K-point index or `Kpoint` object 
            bands_range:
                List of band indices to analyze.
            tol_ediff:
                Tolerance on the energy difference (in eV)
        """
        # Extract the k-point index to speed up the calls belows
        k = self.kindex(kpoint)
        kpoint = self.kpoints[k]

        # Find the set of degenerate states at the given spin and k-point.
        deg_ebands = self.ebands.degeneracies(spin, k, bands_range, tol_ediff=tol_ediff)

        # Create list of tuples (energy, waves) for each degenerate set.
        deg_ewaves = []
        for e, bands in deg_ebands:
            deg_ewaves.append((e, [self.get_wave(spin, k, band) for band in bands])) 
        #print(deg_ewaves)

        # Find the little group of the k-point
        #ltgk = spgrp.get_little_group(kpoint)

        # Compute the D(S) matrices for each degenerate subset.
        #dmats = compute_dmats(ltgk, deg_ewaves)

        #for idg, (e, waves) in enumerate(deg_ewaves):
        #    for symmops in ltgk.groupby_class()
        #        op = symmops[0]
        #        row = -1
        #        for (ii,jj), (wave1, wave2) in iuptri(waves):
        #            if ii != row:
        #               ii = row
        #               rot_wave1 = wave1.rotate(op)
        #            prod = wave2.product(rot_wave1)
        #            dmats[idg][symmop][ii,jj] = prod

        # Locate the D(S) in the lookup table.
        #for trace_atol in [0.1, 0.01, 0.001]:
        #try
        #   dmats.classify(trace_atol)
        #   return dmats
        #except dmats.ClassificationError:
        #   pass
        #
        # Case with accidental degeneraties. 
        # Try to decompose the reducible representations
        #return dmats.decompose()

    #def visualize_ur2(self, spin, kpoint, band, visualizer):
    #    """
    #    Visualize :math:`|u(r)|^2`  with visualizer.
    #    See :class:`Visualizer` for the list of applications and formats supported.
    #    """
    #    extensions = Visualizer.exts_from_appname(visualizer)
    #
    #    for ext in extensions:
    #        ext = "." + ext
    #        try:
    #            return self.export_ur2(ext, spin, kpoint, band)
    #        except Visualizer.Error:
    #            pass
    #    else:
    #        msg = "Don't know how to export data for visualizer %s" % visualizer
    #        raise Visualizer.Error(msg)



class WFK_Reader(ElectronsReader):
    """This object reads data from the WFK file."""

    def __init__(self, filepath):
        """Initialize the object from a filename."""
        super(WFK_Reader, self).__init__(filepath)

        self.kpoints = self.read_kpoints()

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
            return self.kpoints.index(kpoint)

    def read_gvecs_istwfk(self, kpoint):
        """
        Read the set of G-vectors and the value of istwfk for the given k-point.
        Accepts `Kpoint` object or integer.
        """
        k = self.kindex(kpoint)
        npw_k, istwfk = self.npwarr[k], self.istwfk[k]
        return self._kg[k, :npw_k, :], istwfk

    def read_ug(self, spin, kpoint, band):
        """Read the Fourier components of the wavefunction."""
        k = self.kindex(kpoint)
        npw_k, istwfk = self.npwarr[k], self.istwfk[k]
        # TODO use variables to avoid storing the full block.
        return self.set_of_ug[spin, k, band, :, :npw_k]

