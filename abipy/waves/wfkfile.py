# coding: utf-8
"""Wavefunction file."""
from __future__ import print_function, division, unicode_literals, absolute_import

import six
import numpy as np

from monty.functools import lazy_property
from monty.string import marquee # is_string, list_strings,
from abipy.core import Mesh3D, GSphere, Structure
from abipy.core.mixins import AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.iotools import ETSF_Reader, Visualizer
from abipy.electrons.ebands import ElectronsReader
from abipy.waves.pwwave import PWWaveFunction
from abipy.tools import duck

__all__ = [
    "WfkFile",
]


class WfkFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This object provides a simple interface to access and analyze
    the data stored in the WFK file produced by ABINIT.

    Usage example:

    .. code-block:: python

        wfk = WfkFile("foo_WFK.nc")

        # Plot band energies.
        wfk.ebands.plot_ebands()

        # Visualize crystalline structure with vesta.
        wfk.visualize_structure_with("vesta")

        # Visualize u(r)**2 with vesta.
        wfk.visualize_ur2(spin=0, kpoint=0, band=0, appname="vesta")

        # Get a wavefunction.
        wave = wfk.get_wave(spin=0, kpoint=[0, 0, 0], band=0)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: WfkFile
    """
    def __init__(self, filepath):
        """
        Initialize the object from a Netcdf file.
        """
        super(WfkFile, self).__init__(filepath)
        self.reader = reader = WFK_Reader(filepath)
        assert reader.has_pwbasis_set

        # Read the electron bands
        self._ebands = reader.read_ebands()

        self.npwarr = reader.npwarr
        self.nband_sk = reader.nband_sk

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

    def close(self):
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @property
    def ebands(self):
        """|ElectronBands| object"""
        return self._ebands

    @property
    def nkpt(self):
        """Number of k-points."""
        return len(self.kpoints)

    @property
    def gspheres(self):
        """List of :class:`GSphere` objects ordered by k-points."""
        return self._gspheres

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """
        String representation

        Args:
            verbose: verbosity level.
        """
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))

        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)

    def kindex(self, kpoint):
        """The index of the k-point in the file. Accepts :class:`Kpoint` object or int."""
        return self.reader.kindex(kpoint)

    def get_wave(self, spin, kpoint, band):
        """
        Read and return the wavefunction with the given spin, band and kpoint.

        Args:
            spin: spin index. Must be in (0, 1)
            kpoint: Either :class:`Kpoint` instance or integer giving the sequential index in the IBZ (C-convention).
            band: band index.

            returns:
                :class:`WaveFunction` instance.
        """
        ik = self.kindex(kpoint)

        if (spin not in range(self.nsppol) or ik not in range(self.nkpt) or
            band not in range(self.nband_sk[spin, ik])):
            raise ValueError("Wrong (spin, band, kpt) indices")

        ug_skb = self.reader.read_ug(spin, kpoint, band)

        # Istantiate the wavefunction object and set the FFT mesh
        # using the divisions reported in the WFK file.
        wave = PWWaveFunction(self.structure, self.nspinor, spin, band, self.gspheres[ik], ug_skb)
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
            return wave.export_ur2(filepath)
        else:
            return wave.export_ur2(filepath, visu=visu)

    def visualize_ur2(self, spin, kpoint, band, appname="vesta"):
        """
        Visualize :math:`|u(r)|^2`  with visualizer.
        See :class:`Visualizer` for the list of applications and formats supported.
        """
        visu = Visualizer.from_name(appname)

        for ext in visu.supported_extensions():
            ext = "." + ext
            try:
                v = self.export_ur2(ext, spin, kpoint, band, visu=visu)
                return v()
            except visu.Error:
                pass
        else:
            raise visu.Error("Don't know how to export data for visualizer %s" % appname)

    #def classify_states(self, spin, kpoint, band_range=None, energy_range=None, atol=1e-3):
    #    """
    #    Classify electronic eigenstates

    #    Args:
    #        spin: spin index. Must be in (0, 1)
    #        kpoint: Either :class:`Kpoint` instance or integer giving the sequential index in the IBZ (C-convention).
    #        band_range: Define the set of bands included in the classification. See also `energy_range`.
    #            `None` means all bands available in the WFK file. Accepts also range object, use e.g.
    #            `band_range=(0, 5)` to include [0, 1, 2, 3, 4].
    #        energy_range: Define the set of bands included in the classification. Mutually exclusive with `band_range`.
    #            Accepts: "gap"
    #        atol: Absolute tolerance on the energy in eV. Two states are considered degerate if
    #            their energy differ by less than `atol`.

    #    Return:
    #    """
    #    if band_range is not None and energy_range is not None:
    #        raise ValueError("band_range and energy_range are mutually exclusive.")

    #    if energy_range is None:
    #        bids = ... if band_range is None else list(range(band_range))
    #    else:
    #        raise NotImplementedError("energy_range")

    #    # Find little group of the k-point.
    #    ik = self.kindex(kpoint)
    #    kpoint = self.kpoints[ik]
    #    lgk = self.structure.abi_spacegroup.find_little_group(kpoint)
    #    #wclass = WavefuntionsClassifier(kpoint)

    #    # Select bands to include, find degenerate states and group them.
    #    from abipy.core.skw import find_degs_sk
    #    enes = self.ebands.eigens[spin, ik, bids]
    #    deg_list = find_degs_sk(enes, atol)

    #    for deg in deg_list:
    #        print("deg", deg)
    #        nb = len(deg)
    #        waves = [self.get_wave(spin, ik, band) for band in deg]
    #        mats = []
    #        for iclass, op_class in enumerate(lgk.groupby_class()):
    #            op = op_class[0]
    #            #print("op_class[0]\n", op)
    #            #op_waves = [wave.rotate(op) for wave in waves]
    #            op_waves = waves

    #            # Compute <u1|Op|u2>
    #            cmat = np.empty((nb, nb), dtype=np.complex)
    #            for i in range(nb):
    #                #print("int", waves[i].norm2())
    #                for j in range(nb):
    #                    cmat[i, j] = waves[i].braket(op_waves[j], space="r")
    #                    #cmat[i, j] = 1.0

    #            print(cmat)
    #            mats.append(cmat)

    #        #if wclass.classify_characters(deg, mats) /= 0:
    #        #    cprint("Warning")
    #        #    wclass.classify_characters_accidental(deg, full_mats)

    #    #return wclass

    def ipw_visualize_widget(self): # pragma: no cover
        """
        Return an ipython widget with controllers to visualize the wavefunctions.

        .. warning::

            It seems there's a bug with Vesta on MacOs if the user tries to open multiple wavefunctions
            as the tab in vesta is not updated!
        """
        def wfk_visualize(spin, kpoint, band, appname):
            kpoint = int(kpoint.split()[0])
            self.visualize_ur2(spin, kpoint, band, appname=appname)

        import ipywidgets as ipw
        return ipw.interact_manual(
                wfk_visualize,
                spin=list(range(self.nsppol)),
                kpoint=["%d %s" % (i, repr(kpt)) for i, kpt in enumerate(self.kpoints)],
                band=list(range(self.nband)),
                appname=[v.name for v in Visualizer.get_available()],
            )

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("wfk = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(wfk)"),
            nbv.new_code_cell("wfk.ebands.plot();"),
            nbv.new_code_cell("wfk.ebands.kpoints.plot();"),
            nbv.new_code_cell("wfk.ebands.plot();"),
            nbv.new_code_cell("""\
if wfk.ebands.kpoints.is_ibz:
    wfk.ebands.get_edos().plot();"""),
            nbv.new_code_cell("wfk.ipw_visualize_widget()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class WFK_Reader(ElectronsReader):
    """
    This object reads data from the WFK file.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: Wfk_Reader
    """

    def __init__(self, filepath):
        """Initialize the object from a filename."""
        super(WFK_Reader, self).__init__(filepath)

        self.kpoints = self.read_kpoints()
        self.nfft1 = self.read_dimvalue("number_of_grid_points_vector1")
        self.nfft2 = self.read_dimvalue("number_of_grid_points_vector2")
        self.nfft3 = self.read_dimvalue("number_of_grid_points_vector3")

        self.cplex_ug = self.read_dimvalue("real_or_complex_coefficients")
        assert self.cplex_ug == 2

        self.nspinor = self.read_dimvalue("number_of_spinor_components")
        self.nsppol = self.read_dimvalue("number_of_spins")
        self.nspden = self.read_dimvalue("number_of_components")

        self.ecut = self.read_value("kinetic_energy_cutoff")
        self.nband_sk = self.read_value("number_of_states")
        self.istwfk = self.read_value("istwfk")
        self.npwarr = self.read_value("number_of_coefficients")

        # Store G-vectors
        self._kg = self.read_value("reduced_coordinates_of_plane_waves")

    @lazy_property
    def basis_set(self):
        """String defining the basis set."""
        basis_set = self.read_value("basis_set")
        if six.PY2:
            return "".join(basis_set).strip()
        else:
            return "".join(str(basis_set, encoding='UTF-8')).strip()

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

        Accepts: :class:`Kpoint` instance or integer.
        """
        if duck.is_intlike(kpoint):
            return int(kpoint)
        else:
            return self.kpoints.index(kpoint)

    def read_gvecs_istwfk(self, kpoint):
        """
        Read the set of G-vectors and the value of istwfk for the given k-point.
        Accepts :class:`Kpoint` object or integer.
        """
        ik = self.kindex(kpoint)
        npw_k, istwfk = self.npwarr[ik], self.istwfk[ik]
        return self._kg[ik, :npw_k, :], istwfk

    def read_ug(self, spin, kpoint, band):
        """Read the Fourier components of the wavefunction."""
        ik = self.kindex(kpoint)
        npw_k, istwfk = self.npwarr[ik], self.istwfk[ik]
        if self.cplex_ug != 2:
            raise NotImplementedError("")

        # Read data from file (we don't store the full block full block in memory!).
        var = self.rootgrp.variables["coefficients_of_wavefunctions"]
        value = var[spin, ik, band, :, :npw_k, :]
        return value[..., 0] + 1j*value[..., 1]  # Build complex array
