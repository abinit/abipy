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


# Tools and helper functions.
def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


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
        """List of k-points in the WFK file."""
        return self.ebands.kpoints

    @property
    def nkpt(self):
        """Number of k-points."""
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

    def export_ur2(self, filepath, spin, kpoint, band, visu=None):
        """
        Export :math:`|u(r)|^2` on file filename.

        returns:
            Instance of :class:`Visualizer`
        """
        # Read the wavefunction from file.
        wave = self.get_wave(spin, kpoint, band)

        # Export data uding the format specified by filename.
        if visu is None:
            return wave.export_ur2(filepath, self.structure)
        else:
            return wave.export_ur2(filepath, self.structure, visu=visu)

    def visualize_ur2(self, spin, kpoint, band, visu_name):
        """
        Visualize :math:`|u(r)|^2`  with visualizer.
        See :class:`Visualizer` for the list of applications and formats supported.
        """
        visu = Visualizer.from_name(visu_name)
    
        for ext in visu.supported_extensions():
           ext = "." + ext
           try:
               return self.export_ur2(ext, spin, kpoint, band, visu=visu)
           except visu.Error:
               pass
        else:
            raise visu.Error("Don't know how to export data for visualizer %s" % visu_name)

    def classify_ebands(self, spin, kpoint, bands_range, tol_ediff=1e-3):
        """
        Analyze the caracter of the bands at the given k-point and spin.

        Args:
            spin:
                Spin index.
            kpoint:
                K-point index or `Kpoint` object 
            bands_range:
                Range of band indices to analyze.
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

        print("degeneracies detected:", deg_ebands)
        #print(deg_ewaves)

        # Find the little group of the k-point
        ltk = self.structure.spacegroup.find_little_group(kpoint)
        #assert ltk.is_group()

        # Compute the D(R) matrices for each degenerate subset.
        dmats = DMatrices(ltk, deg_ewaves)
        return dmats


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
            self.ug_block = self.read_value("coefficients_of_wavefunctions", cmode="c")
        else:
            raise NotImplementedError("")

    @property
    def basis_set(self):
        """String defining the basis set."""
        try:
            return self._basis_set
        except AttributeError:
            basis_set = self.read_value("basis_set")
            self._basis_set = "".join(basis_set).strip()
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
        return self.ug_block[spin, k, band, :, :npw_k]



class DmatsError(Exception):
    """Base error class."""


class DmatsClassificationError(DmatsError):
    """Bands cannot be classified."""


class DmatsDecompositionError(DmatsError):
    """Accidental degeneracies detected, decompostion failed."""


class DMatrices(object):
    """
    * Let M(R_t) the irreducible representation associated to the space group symmetry (R_t).
    * By convention M(R_t) multiplies wave functions as a row vector:

       $ R_t \psi_a(r) = \psi_a (R^{-1}(r-\tau)) = \sum_b M(R_t)_{ba} \psi_b $

      Therefore, if R_t belongs to the little group of k (i.e. Sk=k+G0), one obtains:

       $ M_ab(R_t) = e^{-i(k+G0).\tau} \int e^{iG0.r} u_{ak}(r)^* u_{bk}(R^{-1}(r-\tau)) \,dr $.

    * The irreducible representation of the small _point_ group of k, M_ab(R), suffices to
      classify the degenerate eigenstates provided that particular conditions are fulfilled 
      (see limitations below). The matrix is indeed given by:

       $ M_ab(R) = e^{+ik.\tau} M_ab(R_t) = e^{-iG0.\tau} \int e^{iG0.r} u_{ak}(r)^* u_{bk}(R^{-1}(r-\tau))\,dr $
     
      The phase factor outside the integral should be zero since symmetry analysis at border zone in non-symmorphic
      space groups is not available. Anyway it is included in our expressions for the sake of consistency.

    * For PAW there is an additional onsite terms involving <phi_i|phi_j(R^{-1}(r-\tau)> and 
      the pseudized version that can be  evaluated using the rotation matrix for 
       real spherical harmonis, zarot(mp,m,l,R). $ Y_{lm}(Rr)= \sum_{m'} zarot(m',m,ll,R) Y_{lm'}(r) $

       $ M^{onsite}_ab(R_t) = sum_{c ij} <\tpsi_a| p_i^c>  <p_j^{c'}|\tpsi_b\> \times 
          [ <\phi_i^c|\phi_j^{c'}> - <\tphi_i^c|\tphi_j^{c'}> ]. $

       $ [ <\phi_i^c|\phi_j^{c'}> - <\tphi_i^c|\tphi_j^{c'}> ] = s_{ij} D_{\mi\mj}^\lj(R^{-1}) $

      where c' is the rotated atom i.e c' = R^{-1}( c-\tau) and D is the rotation matrix for 
      real spherical harmonics.

      Remember that zarot(m',m,l,R)=zarot(m,m',l,R^{-1})
      and $ Y^l_m(ISG) = sum_{m'} D_{m'm}(S) Y_{m'}^l(G) (-i)^l $
          $ D_{m'm}^l (R) = D_{m,m'}^l (R^{-1}) $

    * LIMITATIONS: The method does not work if k is at zone border and the little group of k 
                   contains a non-symmorphic fractional translation. 
    """
    Error = DmatsError
    ClassificationError = DmatsClassificationError
    DecompositionError = DmatsDecompositionError

    def __init__(self, ltk, deg_ewaves):
        """
        Compute the D(R) matrices for each degenerate subset.
        """
        kgroup = ltk.kgroup
        #print("g0vecs", ltk.g0vecs)

        #The main problem here is represented by the fact 
        #that the classes in ltk might not have the 
        #same order as the classes reported in the Bilbao database.
        #Hence we have to shuffle the last dimension of my_character
        #so that we can compare the two array correctly
        #The most robust approach consists in matching class invariants 
        #such as the trace, the determinant of the rotation matrices.
        #as well as the order of the rotation and inv_root!
        #Perhaps I can make LatticeRotation hashable with
        #__hash__ = return det + 10 * trace + 100 * order + 1000 * inv_root
        #The pseudo code below computes the table to_bilbao_classes:

        # Get the Bilbao entry for this point group
        from abipy.core.symmetries import bilbao_ptgroup
        try:
            self.bilbao_ptg = bilbao_ptgroup(kgroup.sch_symbol)
            #to_bilbao_classes = self.bilbao_ptg.map_rotclasses(kgroup.rotclasses)
        except:
            raise self.Error(straceback())

        # Number of degeneracies, 
        # number of spatial rotation in the group of k, number of classes in kgroup.
        self.num_degs = num_degs = len(deg_ewaves)
        num_rotk, num_classes = len(kgroup), kgroup.num_classes

        # Allocate D(R) for each set of degenerated bands.
        # Init them with np.inf. If everything goes well 
        # only the diagonal element for the first operation 
        # in each class has to be computed.
        dmats = self.num_degs * [None]
        for idg, (e, waves) in enumerate(deg_ewaves):
            nb = len(waves)
            dr_bbp = np.empty((num_rotk, nb, nb), dtype=np.complex)
            dr_bbp[:,:,:] = np.inf
            dmats[idg] = dr_bbp

        # FIXME
        # Here there's a problem since I have to exclude time_rev and afm.
        #ltk_symmops = ltk.symmops[:48]
        ltk_symmops = ltk.symmops[:8]

        for idg, (e, waves) in enumerate(deg_ewaves):
            nb = len(waves)
            for isym, symmop in enumerate(ltk_symmops):
                for (ii, wave) in enumerate(waves):
                    rot_wave = wave.rotate(symmop)
                    prod = wave.braket(rot_wave)
                    # Update the entry
                    dmats[idg][isym,ii,ii] = prod
            print("idg", idg, "shape", dmats[idg].shape)

        self.dmats = dmats

        for idg in range(num_degs):
            if idg != 2: continue
            print(self.all_traces(idg))
        
        #my_character[idg] = my_character[idg][to_bilbao_classes]
        #from pymatgen.util.num_utils import iuptri
        # Loop over the set of degenerate states.
        # For each degenerate set compute the full D_(R)
        # mats = [None] * num_degs 

        #for idg, (e, waves) in enumerate(deg_ewaves):
        #    nb = len(waves)
        #    for icl, symmops in enumerate(ltk.groupby_class())
        #        op = symmops[0]
        #        row = -1
        #        for (ii,jj), (wave1, wave2) in iuptri(waves):
        #            if ii != row:
        #               ii = row
        #               rot_wave1 = wave1.rotate(op)
        #            prod = wave2.product(rot_wave1)
        #            # Update the entry
        #            mats[idg][icl][ii,jj] = prod
        #self.mats = mats

        # Locate the D(R) in the lookup table.
        #deg_labels = None
        #for trace_atol in [0.001, 0.01, 0.01]:
        #    try:
        #       deg_labels = dmats.classify(trace_atol)
        #    except dmats.ClassificationError:
        #       pass

        #if deg_labels is not None:
        #    self.deg_labels = deg_labels
        #    return

        # Case with accidental degeneraties. 
        # Try to decompose the reducible representations
        #try:
        #    dmats.decompose()
        #    return dmats
        #except dmats.DecompositionError:
        #    raise 

    #def __str__(self):
    #    lines = []
    #    app = lines.append
    #    return "\n".join(lines)

    #def _compute_diagonal_per_class(self):
    #def _compute_all_dmats(self):

    def all_traces(self, idg):
        """Return the calculated character given the degeneracy index."""
        ds = self.dmats[idg]
        return [d.trace() for d in ds]

    #def my_character(self, idg):
    #    """Return the calculated character given the degeneracy index."""
    #    ds = self.dmats[idg]
    #    return ds[isym,:,:].trace() for isym in class_ids

    def classify(self, trace_rtol=1e-05, trace_atol=1e-08):
        """
        Raises:
            ClassificationError
        """
        # Now my_character has the same ordered as in the Bilbao database.
        # The remaining part is trivial and we only have to deal with possible
        # failures due to numerical errors and accidental degeneracies.
        num_degs = self.num_degs
        deg_labels = [None] * num_degs

        for idg in range(self.num_degs):
            mychar = self.my_character(idg)
            for irrep in self.bilbao_ptg.irreps:
                if np.allclose(irrep.character, mychar, rtol=trace_rtol, atol=trace_atol):
                    deg_labels[idg] = irrep.name
                    break
                                                                                  
        if any(label is None for label in deg_labels):
            raise ClassificationError()

        return deg_labels

    def decompose(self):
        """
        Raises:
            DecompositionError
        """
        raise NotImplementedError()
