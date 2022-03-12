# coding: utf-8
from __future__ import annotations

import functools
import numpy as np
import itertools
import pickle
import os
import json
import warnings
import pandas as pd
import abipy.core.abinit_units as abu

from collections import OrderedDict
from typing import Any, List
from monty.string import is_string, list_strings, marquee
from monty.collections import dict2namedtuple
from monty.functools import lazy_property
from monty.termcolor import cprint
from pymatgen.core.units import eV_to_Ha, Energy
from pymatgen.core.periodic_table import Element
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos as PmgCompletePhononDos, PhononDos as PmgPhononDos
from abipy.core.func1d import Function1D
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_PhononBands, NotebookWriter
from abipy.core.kpoints import Kpoint, Kpath, KpointList, kmesh_from_mpdivs
from abipy.core.structure import Structure
from abipy.abio.robots import Robot
from abipy.iotools import ETSF_Reader
from abipy.tools import duck
from abipy.tools.numtools import gaussian, sort_and_groupby
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims, get_axarray_fig_plt, set_visible,\
    set_ax_xylabels, get_figs_plotly, get_fig_plotly, add_plotly_fig_kwargs, plotlyfigs_to_browser,\
    push_to_chart_studio, PlotlyRowColDesc, plotly_klabels, plotly_set_xylabels, plotly_set_lims
from .phtk import match_eigenvectors, get_dyn_mat_eigenvec, open_file_phononwebsite, NonAnalyticalPh


__all__ = [
    "PhononBands",
    "PhononBandsPlotter",
    "PhbstFile",
    "PhononDos",
    "PhononDosPlotter",
    "PhdosReader",
    "PhdosFile",
]


@functools.total_ordering
class PhononMode:
    """
    A phonon mode has a q-point, a frequency, a cartesian displacement and a |Structure|.
    """

    __slots__ = [
        "qpoint",
        "freq",
        "displ_cart", # Cartesian displacement.
        "structure"
    ]

    def __init__(self, qpoint, freq, displ_cart, structure: Structure) -> None:
        """
        Args:
            qpoint: qpoint in reduced coordinates.
            freq: Phonon frequency in eV.
            displ_cart: Displacement (Cartesian coordinates in Angstrom)
            structure: |Structure| object.
        """
        self.qpoint = Kpoint.as_kpoint(qpoint, structure.reciprocal_lattice)
        self.freq = freq
        self.displ_cart = displ_cart
        self.structure = structure

    # Rich comparison support (ordered is based on the frequency).
    # Missing operators are automatically filled by total_ordering.
    def __eq__(self, other):
        return self.freq == other.freq

    def __lt__(self, other):
        return self.freq < other.freq

    def __str__(self):
        return self.to_string(with_displ=False)

    def to_string(self, with_displ=True, verbose=0) -> str:
        """
        String representation

        Args:
            verbose: Verbosity level.
            with_displ: True to print phonon displacement.
        """
        lines = ["%s: q-point %s, frequency %.5f (eV)" % (self.__class__.__name__, self.qpoint, self.freq)]
        app = lines.append

        if with_displ:
            app("Phonon displacement in cartesian coordinates [Angstrom]")
            app(str(self.displ_cart))

        return "\n".join(lines)

    #@property
    #def displ_red(self)
    #    return np.dot(self.xred, self.rprimd)

    #def export(self, path):
    #def visualize(self, visualizer):
    #def build_supercell(self):


class PhononBands:
    """
    Container object storing the phonon band structure.

    .. note::

        Frequencies are in eV. Cartesian displacements are in Angstrom.
    """

    @classmethod
    def from_file(cls, filepath: str) -> PhononBands:
        """Create the object from a netcdf_ file."""
        with PHBST_Reader(filepath) as r:
            structure = r.read_structure()

            # Build the list of q-points
            qpoints = Kpath(structure.reciprocal_lattice, frac_coords=r.read_qredcoords(),
                            weights=r.read_qweights(), names=None)

            for qpoint in qpoints:
                qpoint.set_name(structure.findname_in_hsym_stars(qpoint))

            # Read amu
            amu_list = r.read_amu()
            if amu_list is not None:
                atomic_numbers = r.read_value("atomic_numbers")
                amu = {at: a for at, a in zip(atomic_numbers, amu_list)}
            else:
                cprint("Warning: file %s does not contain atomic_numbers.\nParticular methods need them!" %
                       filepath, "red")
                amu = None

            non_anal_ph = None

            # Reading NonAnalyticalPh here is not 100% safe as it may happen that the netcdf file
            # does not contain all the directions required by AbiPy.
            # So we read NonAnalyticalPh only if we know that all directions are available.
            # The flag has_abipy_non_anal_ph is set at the Fortran level. See e.g ifc_mkphbs
            if ("non_analytical_directions" in r.rootgrp.variables and
                "has_abipy_non_anal_ph" in r.rootgrp.variables):
                #print("Found non_anal_ph term compatible with AbiPy plotter.")
                non_anal_ph = NonAnalyticalPh.from_file(filepath)

            epsinf, zcart = r.read_epsinf_zcart()

            return cls(structure=structure,
                       qpoints=qpoints,
                       phfreqs=r.read_phfreqs(),
                       phdispl_cart=r.read_phdispl_cart(),
                       amu=amu,
                       non_anal_ph=non_anal_ph,
                       epsinf=epsinf, zcart=zcart,
                       )

    @classmethod
    def as_phbands(cls, obj: Any) -> PhononBands:
        """
        Return an instance of |PhononBands| from a generic object ``obj``.
        Supports:

            - instances of cls
            - files (string) that can be open with ``abiopen`` and that provide a ``phbands`` attribute.
            - objects providing a ``phbands`` attribute.
        """
        if isinstance(obj, cls):
            return obj

        elif is_string(obj):
            # path?
            if obj.endswith(".pickle"):
                with open(obj, "rb") as fh:
                    return cls.as_phbands(pickle.load(fh))

            from abipy.abilab import abiopen
            with abiopen(obj) as abifile:
                return abifile.phbands

        elif hasattr(obj, "phbands"):
            # object with phbands
            return obj.phbands

        raise TypeError("Don't know how to extract a PhononBands from type %s" % type(obj))

    @staticmethod
    def phfactor_ev2units(units: str) -> float:
        """
        Return conversion factor eV --> units (case-insensitive)
        """
        return abu.phfactor_ev2units(units)

    def read_non_anal_from_file(self, filepath: str) -> None:
        """
        Reads the non analytical directions, frequencies and displacements from the anaddb.nc file
        specified and adds them to the object.
        """
        self.non_anal_ph = NonAnalyticalPh.from_file(filepath)

    def set_phonopy_obj_from_ananc(self, ananc, supercell_matrix, symmetrize_tensors=False,
                                   symprec=1e-5, set_masses=True):
        """
        Generates the Phonopy object from an anaddb.nc file that contains the interatomic force constants.
        Based on the converter implemented in abipy.dfpt.converters.

        Args:
            ananc: a string with the path to the anaddb.nc file or an instance of AnaddbNcFile.
            supercell_matrix: the supercell matrix used for phonopy. Any choice is acceptable, however
                the best agreement between the abinit and phonopy results is obtained if this is set to
                a diagonal matrix with on the diagonal the ngqpt used to generate the anaddb.nc.
            symmetrize_tensors: if True the tensors will be symmetrized in the Phonopy object and
                in the output files. This will apply to IFC, BEC and dielectric tensor.
            symprec: distance tolerance in Cartesian coordinates to find crystal symmetry in phonopy.
            set_masses: if True the atomic masses used by abinit will be added to the PhonopyAtoms
                and will be present in the returned Phonopy object. This should improve compatibility
                among abinit and phonopy results if frequencies needs to be calculated.
        """
        from abipy.dfpt.anaddbnc import AnaddbNcFile
        from abipy.dfpt.converters import abinit_to_phonopy

        if isinstance(ananc, str):
            with AnaddbNcFile(ananc) as anc:
                ph = abinit_to_phonopy(anc, supercell_matrix=supercell_matrix, symmetrize_tensors=symmetrize_tensors,
                                       symprec=symprec, set_masses=set_masses)
        else:
            ph = abinit_to_phonopy(ananc, supercell_matrix=supercell_matrix, symmetrize_tensors=symmetrize_tensors,
                                   symprec=symprec, set_masses=set_masses)
        self.phonopy_obj = ph

    def __init__(self, structure, qpoints, phfreqs, phdispl_cart, non_anal_ph=None, amu=None,
                 epsinf=None, zcart=None, linewidths=None, phonopy_obj=None):
        """
        Args:
            structure: |Structure| object.
            qpoints: |KpointList| instance.
            phfreqs: Phonon frequencies in eV.
            phdispl_cart: [nqpt, 3*natom, 3*natom] array with displacement in Cartesian coordinates in Angstrom.
                The last dimension stores the cartesian components.
            non_anal_ph: :class:`NonAnalyticalPh` with information of the non analytical contribution
                None if contribution is not present.
            amu: dictionary that associates the atomic species present in the structure to the values of the atomic
                mass units used for the calculation.
            epsinf: [3,3] matrix with electronic dielectric tensor in Cartesian coordinates.
                None if not avaiable.
            zcart: [natom, 3, 3] matrix with Born effective charges in Cartesian coordinates.
                None if not available.
            linewidths: Array-like object with the linewidths (eV) stored as [q, num_modes]
            phonopy_obj: an instance of a Phonopy object obtained from the same IFC used to generate the
                band structure.
        """
        self.structure = structure

        # KpointList with the q-points
        self.qpoints = qpoints
        self.num_qpoints = len(self.qpoints)

        # numpy array with phonon frequencies. Shape=(nqpt, 3*natom)
        self.phfreqs = phfreqs

        # phonon displacements in Cartesian coordinates.
        # `ndarray` of shape (nqpt, 3*natom, 3*natom).
        # The last dimension stores the cartesian components.
        self.phdispl_cart = phdispl_cart

        # Handy variables used to loop.
        self.num_atoms = structure.num_sites
        self.num_branches = 3 * self.num_atoms
        self.branches = range(self.num_branches)

        self.non_anal_ph = non_anal_ph
        self.amu = amu
        self.amu_symbol = None
        if amu is not None:
            self.amu_symbol = {}
            for z, m in amu.items():
                el = Element.from_Z(int(z))
                self.amu_symbol[el.symbol] = m

        self._linewidths = None
        if linewidths is not None:
            self._linewidths = np.reshape(linewidths, self.phfreqs.shape)

        self.epsinf = epsinf
        self.zcart = zcart
        self.phonopy_obj = phonopy_obj

        # Dictionary with metadata e.g. nkpt, tsmear ...
        self.params = OrderedDict()

    # TODO: Replace num_qpoints with nqpt, deprecate num_qpoints
    @property
    def nqpt(self) -> int:
        """An alias for num_qpoints."""
        return self.num_qpoints

    def __repr__(self):
        """String representation (short version)"""
        return "<%s, nk=%d, %s, id=%s>" % (
                self.__class__.__name__, self.num_qpoints, self.structure.formula, id(self))

    def __str__(self):
        return self.to_string()

    def to_string(self, title=None, with_structure=True, with_qpoints=False, verbose=0) -> str:
        """
        Human-readable string with useful information such as structure, q-points, ...

        Args:
            with_structure: False if structural info should not be displayed.
            with_qpoints: False if q-point info shoud not be displayed.
            verbose: Verbosity level.
        """
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))

        if with_structure:
            app(self.structure.to_string(verbose=verbose, title="Structure"))
            app("")

        #app(marquee("Phonon Bands", mark="="))
        app("Number of q-points: %d" % self.num_qpoints)
        app("Atomic mass units: %s" % str(self.amu))
        has_dipdip = self.non_anal_ph is not None
        app("Has non-analytical contribution for q --> 0: %s" % has_dipdip)
        if verbose and has_dipdip:
            app(str(self.non_anal_ph))

        if with_qpoints:
            app(self.qpoints.to_string(verbose=verbose, title="Q-points"))
            app("")

        return "\n".join(lines)

    def __add__(self, other: PhononBands) -> PhononBandsPlotter:
        """self + other returns a |PhononBandsPlotter| object."""
        if not isinstance(other, (PhononBands, PhononBandsPlotter)):
            raise TypeError("Cannot add %s to %s" % (type(self), type(other)))

        if isinstance(other, PhononBandsPlotter):
            self_key = repr(self)
            other.add_phbands(self_key, self)
            return other
        else:
            plotter = PhononBandsPlotter()
            self_key = repr(self)
            plotter.add_phbands(self_key, self)
            self_key = repr(self)
            other_key = repr(other)
            plotter.add_phbands(other_key, other)
            return plotter

    __radd__ = __add__

    @lazy_property
    def _auto_qlabels(self):
        # Find the q-point names in the pymatgen database.
        # We'll use _auto_qlabels to label the point in the matplotlib plot
        # if qlabels are not specified by the user.
        _auto_qlabels = OrderedDict()

        # If the first or the last q-point are not recognized in findname_in_hsym_stars
        # matplotlib won't show the full band structure along the k-path
        # because the labels are not defined. Here we make sure that
        # the labels for the extrema of the path are always defined.
        _auto_qlabels[0] = " "

        for idx, qpoint in enumerate(self.qpoints):
            name = qpoint.name if qpoint.name is not None else self.structure.findname_in_hsym_stars(qpoint)
            if name is not None:
                _auto_qlabels[idx] = name
                if qpoint.name is None: qpoint.set_name(name)

        last = len(self.qpoints) - 1
        if last not in _auto_qlabels: _auto_qlabels[last] = " "

        return _auto_qlabels

    @property
    def displ_shape(self):
        """The shape of phdispl_cart."""
        return self.phdispl_cart.shape

    @property
    def minfreq(self) -> float:
        """Minimum phonon frequency."""
        return self.get_minfreq_mode()

    @property
    def maxfreq(self) -> float:
        """Maximum phonon frequency in eV."""
        return self.get_maxfreq_mode()

    def get_minfreq_mode(self, mode=None):
        """Compute the minimum of the frequencies."""
        if mode is None:
            return np.min(self.phfreqs)
        else:
            return np.min(self.phfreqs[:, mode])

    def get_maxfreq_mode(self, mode=None):
        """Compute the minimum of the frequencies."""
        if mode is None:
            return np.max(self.phfreqs)
        else:
            return np.max(self.phfreqs[:, mode])

    @property
    def shape(self):
        """Shape of the array with the eigenvalues."""
        return self.num_qpoints, self.num_branches

    @property
    def linewidths(self):
        """linewidths in eV. |numpy-array| with shape [nqpt, num_branches]."""
        return self._linewidths

    @linewidths.setter
    def linewidths(self, linewidths):
        """Set the linewidths. Accept real array of shape [nqpt, num_branches] or None."""
        if linewidths is not None:
            linewidths = np.reshape(linewidths, self.shape)
        self._linewidths = linewidths

    @property
    def has_linewidths(self) -> bool:
        """True if bands with linewidths."""
        return getattr(self, "_linewidths", None) is not None

    @lazy_property
    def dyn_mat_eigenvect(self):
        """
        [nqpt, 3*natom, 3*natom] array with the orthonormal eigenvectors of the dynamical matrix.
        in Cartesian coordinates.
        """
        return get_dyn_mat_eigenvec(self.phdispl_cart, self.structure, amu=self.amu)

    @property
    def non_anal_directions(self):
        """Cartesian directions along which the non analytical frequencies and displacements are available"""
        if self.non_anal_ph:
            return self.non_anal_ph.directions
        else:
            return None

    @property
    def non_anal_phfreqs(self):
        """Phonon frequencies with non analytical contribution in eV along non_anal_directions"""
        if self.non_anal_ph:
            return self.non_anal_ph.phfreqs
        else:
            return None

    @property
    def non_anal_phdispl_cart(self):
        """Displacement in Cartesian coordinates with non analytical contribution along non_anal_directions"""
        if self.non_anal_ph:
            return self.non_anal_ph.phdispl_cart
        else:
            return None

    @property
    def non_anal_dyn_mat_eigenvect(self):
        """Eigenvalues of the dynamical matrix with non analytical contribution along non_anal_directions."""
        if self.non_anal_ph:
            return self.non_anal_ph.dyn_mat_eigenvect
        else:
            return None

    def to_xmgrace(self, filepath: str, units: str = "meV") -> str:
        """
        Write xmgrace_ file with phonon band structure energies and labels for high-symmetry q-points.

        Args:
            filepath: String with filename or stream.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
        """
        is_stream = hasattr(filepath, "write")

        if is_stream:
            f = filepath
        else:
            f = open(filepath, "wt")

        def w(s):
            f.write(s)
            f.write("\n")

        factor = abu.phfactor_ev2units(units)
        wqnu_units = self.phfreqs * factor

        import datetime
        w("# Grace project file with phonon band energies.")
        w("# Generated by AbiPy on: %s" % str(datetime.datetime.today()))
        w("# Crystalline structure:")
        for s in str(self.structure).splitlines():
            w("# %s" % s)
        w("# Energies are in %s." % units)
        w("# List of q-points and their index (C notation i.e. count from 0)")
        for iq, qpt in enumerate(self.qpoints):
            w("# %d %s" % (iq, str(qpt.frac_coords)))
        w("@page size 792, 612")
        w("@page scroll 5%")
        w("@page inout 5%")
        w("@link page off")
        w("@with g0")
        w("@world xmin 0.00")
        w('@world xmax %d' % (self.num_qpoints - 1))
        w('@world ymin %s' % wqnu_units.min())
        w('@world ymax %s' % wqnu_units.max())
        w('@default linewidth 1.5')
        w('@xaxis  tick on')
        w('@xaxis  tick major 1')
        w('@xaxis  tick major color 1')
        w('@xaxis  tick major linestyle 3')
        w('@xaxis  tick major grid on')
        w('@xaxis  tick spec type both')
        w('@xaxis  tick major 0, 0')

        qticks, qlabels = self._make_ticks_and_labels(qlabels=None)
        w('@xaxis  tick spec %d' % len(qticks))
        for iq, (qtick, qlabel) in enumerate(zip(qticks, qlabels)):
            w('@xaxis  tick major %d, %d' % (iq, qtick))
            w('@xaxis  ticklabel %d, "%s"' % (iq, qlabel))

        w('@xaxis  ticklabel char size 1.500000')
        w('@yaxis  tick major 10')
        w('@yaxis  label "Phonon %s"' % abu.phunit_tag(units))
        w('@yaxis  label char size 1.500000')
        w('@yaxis  ticklabel char size 1.500000')
        for nu in self.branches:
            w('@    s%d line color %d' % (nu, 1))

        # TODO: support LO-TO splitting (?)
        for nu in self.branches:
            w('@target G0.S%d' % nu)
            w('@type xy')
            for iq in range(self.num_qpoints):
                w('%d %.8E' % (iq, wqnu_units[iq, nu]))
            w('&')

        if not is_stream:
            f.close()

    # TODO
    #def to_bxsf(self, filepath):
    #    """
    #    Export the full band structure to `filepath` in BXSF format
    #    suitable for the visualization of isosurfaces with Xcrysden (xcrysden --bxsf FILE).
    #    Require q-points in IBZ and gamma-centered q-mesh.
    #    """
    #    self.get_phbands3d().to_bxsf(filepath)

    #def get_phbands3d(self):
    #    has_timrev, fermie = True, 0.0
    #    return PhononBands3D(self.structure, self.qpoints, has_timrev, self.phfreqs, fermie)

    def qindex(self, qpoint) -> int:
        """Returns the index of the qpoint. Accepts integer or reduced coordinates."""
        if duck.is_intlike(qpoint):
            return int(qpoint)
        else:
            return self.qpoints.index(qpoint)

    def qindex_qpoint(self, qpoint, is_non_analytical_direction=False):
        """
        Returns (qindex, qpoint) from an integer or a qpoint.

        Args:
            qpoint: integer, vector of reduced coordinates or |Kpoint| object.
            is_non_analytical_direction: True if qpoint should be interpreted as a fractional direction for q --> 0
                In this case qindex refers to the index of the direction in the :class:`NonAnalyticalPh` object.
        """
        if not is_non_analytical_direction:
            # Standard search in qpoints.
            qindex = self.qindex(qpoint)
            return qindex, self.qpoints[qindex]
        else:
            # Find index of direction given by qpoint.
            if self.non_anal_ph is None:
                raise ValueError("Phononbands does not contain non-analytical terms for q-->0")

            # Extract direction (assumed in fractional coordinates)
            if hasattr(qpoint, "frac_coords"):
                direction = qpoint.frac_coords
            elif duck.is_intlike(qpoint):
                direction = self.non_anal_ph.directions[qpoint]
            else:
                direction = qpoint

            qindex = self.non_anal_ph.index_direction(direction, cartesian=False)

            # Convert to fractional coords.
            cart_direc = self.non_anal_ph.directions[qindex]
            red_direc = self.structure.reciprocal_lattice.get_fractional_coords(cart_direc)
            qpoint = Kpoint(red_direc, self.structure.reciprocal_lattice, weight=None, name=None)

            return qindex, qpoint

    def get_unstable_modes(self, below_mev=-5.0):
        """
        Return a list of :class:`PhononMode` objects with the unstable modes.
        A mode is unstable if its frequency is < below_mev. Output list is sorted
        and modes with lowest frequency come first.
        """
        umodes = []

        for iq, qpoint in enumerate(self.qpoints):
            for nu in self.branches:
                freq = self.phfreqs[iq, nu]
                if freq < below_mev / 1000:
                    displ_cart = self.phdispl_cart[iq, nu, :]
                    umodes.append(PhononMode(qpoint, freq, displ_cart, self.structure))

        return sorted(umodes)

    # TODO
    #def find_irreps(self, qpoint, tolerance):
    #    """
    #    Find the irreducible representation at this q-point
    #    Raise: QIrrepsError if algorithm fails
    #    """
    #    qindex, qpoint = self.qindex_qpoint(qpoint)

    def get_dict4pandas(self, with_spglib=True) -> dict:
        """
        Return a :class:`OrderedDict` with the most important parameters:

            - Chemical formula and number of atoms.
            - Lattice lengths, angles and volume.
            - The spacegroup number computed by Abinit (set to None if not available).
            - The spacegroup number and symbol computed by spglib (set to None not `with_spglib`).

        Useful to construct pandas DataFrames

        Args:
            with_spglib: If True, spglib_ is invoked to get the spacegroup symbol and number
        """
        odict = OrderedDict([
            ("nqpt", self.num_qpoints), ("nmodes", self.num_branches),
            ("min_freq", self.minfreq), ("max_freq", self.maxfreq),
            ("mean_freq", self.phfreqs.mean()), ("std_freq", self.phfreqs.std())

        ])
        odict.update(self.structure.get_dict4pandas(with_spglib=with_spglib))

        return odict

    def get_phdos(self, method: str = "gaussian", step: float = 1.e-4, width: float = 4.e-4) -> PhononDos:
        """
        Compute the phonon DOS on a linear mesh.

        Args:
            method: String defining the method
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Returns:
            |PhononDos| object.

        .. warning::

            Requires a homogeneous sampling of the Brillouin zone.
        """
        if abs(self.qpoints.sum_weights() - 1) > 1.e-6:
            raise ValueError("Qpoint weights should sum up to one")

        # Compute the linear mesh for the DOS
        w_min = self.minfreq
        w_min -= 0.1 * abs(w_min)
        w_max = self.maxfreq
        w_max += 0.1 * abs(w_max)
        nw = 1 + (w_max - w_min) / step

        mesh, step = np.linspace(w_min, w_max, num=nw, endpoint=True, retstep=True)

        values = np.zeros(nw)
        if method == "gaussian":
            for q, qpoint in enumerate(self.qpoints):
                weight = qpoint.weight
                for nu in self.branches:
                    w = self.phfreqs[q, nu]
                    values += weight * gaussian(mesh, width, center=w)

        else:
            raise ValueError("Method %s is not supported" % str(method))

        return PhononDos(mesh, values)

    def create_xyz_vib(self, iqpt, filename, pre_factor=200, do_real=True, scale_matrix=None, max_supercell=None):
        """
        Create vibration XYZ file for visualization of phonons.

        Args:
            iqpt: index of qpoint.
            filename: name of the XYZ file that will be created.
            pre_factor: Multiplication factor of the displacements.
            do_real: True if we want only real part of the displacement, False means imaginary part.
            scale_matrix: Scaling matrix of the supercell.
            max_supercell: Maximum size of the supercell with respect to primitive cell.
        """
        if scale_matrix is None:
            if max_supercell is None:
                raise ValueError("If scale_matrix is None, max_supercell must be provided!")

            scale_matrix = self.structure.get_smallest_supercell(self.qpoints[iqpt].frac_coords,
                                                                 max_supercell=max_supercell)

        natoms = int(np.round(len(self.structure) * np.linalg.det(scale_matrix)))

        with open(filename, "wt") as xyz_file:
            for imode in np.arange(self.num_branches):
                xyz_file.write(str(natoms) + "\n")
                xyz_file.write("Mode " + str(imode) + " : " + str(self.phfreqs[iqpt, imode]) + "\n")
                self.structure.write_vib_file(
                    xyz_file, self.qpoints[iqpt].frac_coords,
                    pre_factor * np.reshape(self.phdispl_cart[iqpt, imode,:],(-1,3)),
                    do_real=True, frac_coords=False, max_supercell=max_supercell, scale_matrix=scale_matrix)

    def create_ascii_vib(self, iqpts, filename, pre_factor=1):
        """
        Create vibration ascii file for visualization of phonons.
        This format can be read with v_sim_ or ascii-phonons.

        Args:
            iqpts: an index or a list of indices of the qpoints in self. Note that at present only V_sim supports
                an ascii file with multiple qpoints.
            filename: name of the ascii file that will be created.
            pre_factor: Multiplication factor of the displacements.
        """
        if not isinstance(iqpts, (list, tuple)):
            iqpts = [iqpts]

        structure = self.structure
        a, b, c = structure.lattice.abc
        alpha, beta, gamma = (np.pi*a/180 for a in structure.lattice.angles)
        m = structure.lattice.matrix
        sign = np.sign(np.dot(np.cross(m[0], m[1]), m[2]))

        dxx = a
        dyx = b * np.cos(gamma)
        dyy = b * np.sin(gamma)
        dzx = c * np.cos(beta)
        dzy = c * (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma)
        # keep the same orientation
        dzz = sign*np.sqrt(c**2-dzx**2-dzy**2)

        lines = ["# ascii file generated with abipy"]
        lines.append("  {: 3.10f}  {: 3.10f}  {: 3.10f}".format(dxx, dyx, dyy))
        lines.append("  {: 3.10f}  {: 3.10f}  {: 3.10f}".format(dzx, dzy, dzz))

        # use reduced coordinates
        lines.append("#keyword: reduced")

        # coordinates
        for s in structure:
            lines.append("  {: 3.10f}  {: 3.10f}  {: 3.10f} {:>2}".format(s.a, s.b, s.c, s.specie.name))

        ascii_basis = [[dxx, 0, 0],
                       [dyx, dyy, 0],
                       [dzx, dzy, dzz]]

        for iqpt in iqpts:
            q = self.qpoints[iqpt].frac_coords

            displ_list = np.zeros((self.num_branches, self.num_atoms, 3), dtype=complex)
            for i in range(self.num_atoms):
                displ_list[:,i,:] = self.phdispl_cart[iqpt,:,3*i:3*(i+1)] * \
                    np.exp(-2*np.pi*1j*np.dot(structure[i].frac_coords, self.qpoints[iqpt].frac_coords))

            displ_list = np.dot(np.dot(displ_list, structure.lattice.inv_matrix), ascii_basis) * pre_factor

            for imode in np.arange(self.num_branches):
                lines.append("#metaData: qpt=[{:.6f};{:.6f};{:.6f};{:.6f} \\".format(
                    q[0], q[1], q[2], self.phfreqs[iqpt, imode]))

                for displ in displ_list[imode]:
                    line = "#; " + "; ".join("{:.6f}".format(i) for i in displ.real) + "; " \
                           + "; ".join("{:.6f}".format(i) for i in displ.imag) + " \\"
                    lines.append(line)

                lines.append(("# ]"))

        with open(filename, 'wt') as f:
            f.write("\n".join(lines))

    def view_phononwebsite(self, browser=None, verbose=0, dryrun=False, **kwargs):
        """
        Produce JSON_ file that can be parsed from the phononwebsite_ and open it in ``browser``.

        Args:
            browser: Open webpage in ``browser``. Use default $BROWSER if None.
            verbose: Verbosity level
            dryrun: Activate dryrun mode for unit testing purposes.
            kwargs: Passed to create_phononwebsite_json method

        Return: Exit status
        """
        # Create json in abipy_nbworkdir with relative path so that we can read it inside the browser.
        from abipy.core.globals import abinb_mkstemp
        prefix = self.structure.formula.replace(" ", "")
        _, rpath = abinb_mkstemp(force_abinb_workdir=not dryrun, use_relpath=True,
                                 prefix=prefix, suffix=".json", text=True)

        if verbose: print("Writing json file:", rpath)
        self.create_phononwebsite_json(rpath, indent=None, **kwargs)

        if dryrun: return 0
        return open_file_phononwebsite(rpath, browser=browser)

    def create_phononwebsite_json(self, filename, name=None, repetitions=None, highsym_qpts=None,
                                  match_bands=True, highsym_qpts_mode="std", indent=2):
        """
        Writes a JSON_ file that can be parsed from the phononwebsite_.

        Args:
            filename: name of the json file that will be created
            name: name associated with the data.
            repetitions: number of repetitions of the cell. List of three integers. Defaults to [3,3,3].
            highsym_qpts: list of tuples. The first element of each tuple should be a list with the coordinates
                of a high symmetry point, the second element of the tuple should be its label.
            match_bands: if True tries to follow the band along the path based on the scalar product of the eigenvectors.
            highsym_qpts_mode: if ``highsym_qpts`` is None, high symmetry q-points can be automatically determined.
                Accepts the following values:
                'split' will split the path based on points where the path changes direction in the Brillouin zone.
                Similar to what is done in phononwebsite. Only Gamma will be labeled.
                'std' uses the standard generation procedure for points and labels used in PhononBands.
                None does not set any point.
            indent: Indentation level, passed to json.dump
        """

        def split_non_collinear(qpts):
            r"""
            function that splits the list of qpoints at repetitions (only the first point will be considered as
            high symm) and where the direction changes. Also sets :math:`\Gamma` for [0, 0, 0].
            Similar to what is done in phononwebsite_.
            """
            h = []
            if np.array_equal(qpts[0], [0, 0, 0]):
                h.append((0, "\\Gamma"))
            for i in range(1, len(qpts)-1):
                if np.array_equal(qpts[i], [0,0,0]):
                    h.append((i, "\\Gamma"))
                elif np.array_equal(qpts[i], qpts[i+1]):
                    h.append((i, ""))
                else:
                    v1 = [a_i - b_i for a_i, b_i in zip(qpts[i+1], qpts[i])]
                    v2 = [a_i - b_i for a_i, b_i in zip(qpts[i-1], qpts[i])]
                    if not np.isclose(np.linalg.det([v1,v2,[1,1,1]]), 0):
                        h.append((i, ""))
            if np.array_equal(qpts[-1], [0, 0, 0]):
                h.append((len(qpts)-1, "\\Gamma"))

            return h

        def reasonable_repetitions(natoms):
            if (natoms < 4): return (3, 3, 3)
            if (4 < natoms < 50): return (2, 2, 2)
            if (50 < natoms): return (1, 1, 1)

        # http://henriquemiranda.github.io/phononwebsite/index.html
        data = {}
        data["name"] = name or self.structure.composition.reduced_formula
        data["natoms"] = self.num_atoms
        data["lattice"] = self.structure.lattice.matrix.tolist()
        data["atom_types"] = [e.name for e in self.structure.species]
        data["atom_numbers"] = self.structure.atomic_numbers
        data["formula"] = self.structure.formula.replace(" ", "")
        data["repetitions"] = repetitions or reasonable_repetitions(self.num_atoms)
        data["atom_pos_car"] = self.structure.cart_coords.tolist()
        data["atom_pos_red"] = self.structure.frac_coords.tolist()
        data["chemical_symbols"] = self.structure.symbol_set
        data["atomic_numbers"] = list(set(self.structure.atomic_numbers))

        qpoints = []
        for q_sublist in self.split_qpoints:
            qpoints.extend(q_sublist.tolist())

        if highsym_qpts is None:
            if highsym_qpts_mode is None:
                data["highsym_qpts"] = []
            elif highsym_qpts_mode == 'split':
                data["highsym_qpts"] = split_non_collinear(qpoints)
            elif highsym_qpts_mode == 'std':
                data["highsym_qpts"] = list(zip(*self._make_ticks_and_labels(None)))
        else:
            data["highsym_qpts"] = highsym_qpts

        distances = [0]
        for i in range(1, len(qpoints)):
            q_coord_1 = self.structure.reciprocal_lattice.get_cartesian_coords(qpoints[i])
            q_coord_2 = self.structure.reciprocal_lattice.get_cartesian_coords(qpoints[i-1])
            distances.append(distances[-1] + np.linalg.norm(q_coord_1-q_coord_2))

        eigenvalues = []
        for i, phfreqs_sublist in enumerate(self.split_phfreqs):
            phfreqs_sublist = phfreqs_sublist * eV_to_Ha * abu.Ha_cmm1
            if match_bands:
                ind = self.split_matched_indices[i]
                phfreqs_sublist = phfreqs_sublist[np.arange(len(phfreqs_sublist))[:, None], ind]
            eigenvalues.extend(phfreqs_sublist.tolist())

        vectors = []

        for i, (qpts, phdispl_sublist) in enumerate(zip(self.split_qpoints, self.split_phdispl_cart)):
            vect = np.array(phdispl_sublist)

            if match_bands:
                vect = vect[np.arange(vect.shape[0])[:, None, None],
                            self.split_matched_indices[i][...,None],
                            np.arange(vect.shape[2])[None, None,:]]
            v = vect.reshape((len(vect), self.num_branches,self.num_atoms, 3))
            norm = [np.linalg.norm(vi) for vi in v[0,0]]
            v /= max(norm)
            v = np.stack([v.real, v.imag], axis=-1)

            vectors.extend(v.tolist())

        data["qpoints"] = qpoints
        data["distances"] = distances
        data["eigenvalues"] = eigenvalues
        data["vectors"] = vectors
        #print("name", data["name"], "\nhighsym_qpts:", data["highsym_qpts"])

        with open(filename, 'wt') as json_file:
            json.dump(data, json_file, indent=indent)

    def make_isodistort_ph_dir(self, qpoint, select_modes=None, eta=1, workdir=None):
        """
        Compute ph-freqs for given q-point (default: Gamma),
        produce CIF files for unperturbed and distorded structure
        that can be used with ISODISTORT (https://stokes.byu.edu/iso/isodistort.php)
        to analyze the symmetry of phonon modes.
        See README.me file produced in output directory.

        Args:
            qpoint:
            wordir:
            select_modes:
            eta: Amplitude of the displacement to be applied to the system. Will correspond to the
                largest displacement of one atom in Angstrom.
            scale_matrix: the scaling matrix of the supercell. If None a scaling matrix suitable for
                the qpoint will be determined.
            max_supercell: mandatory if scale_matrix is None, ignored otherwise. Defines the largest
                supercell in the search for a scaling matrix suitable for the q point.
        """
        iq, qpoint = self.qindex_qpoint(qpoint)

        scale_matrix = np.eye(3, 3, dtype=int)
        important_fracs = (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
        for i in range(3):
            for comparison_frac in important_fracs:
                if abs(1 - qpoint.frac_coords[i] * comparison_frac) < 1e-4:
                    scale_matrix[i, i] = comparison_frac
                    break
        print(f"Using scale_matrix:\n {scale_matrix}")

        select_modes = self.branches if select_modes is None else select_modes
        if workdir is None:
            workdir = "%s_qpt%s" % (self.structure.formula, repr(qpoint))
            workdir = workdir.replace(" ", "_").replace("$", "").replace("\\", "").replace("[", "").replace("]", "")

        if os.path.isdir(workdir):
            cprint(f"Removing pre-existing directory: {workdir}", "yellow")
            import shutil
            shutil.rmtree(workdir)

        os.mkdir(workdir)

        print(f"\nCreating CIF files for ISODISTORT code in {workdir}. See README.md")
        self.structure.write_cif_with_spglib_symms(filename=os.path.join(workdir, "parent_structure.cif"))

        for imode in select_modes:
            # A namedtuple with a structure with the displaced atoms, a numpy array containing the
            # displacements applied to each atom and the scale matrix used to generate the supercell.
            r = self.get_frozen_phonons(qpoint, imode,
                                        eta=eta, scale_matrix=scale_matrix, max_supercell=None)

            print("after scale_matrix:", r.scale_matrix)
            r.structure.write_cif_with_spglib_symms(filename=os.path.join(workdir,
                                                    "distorted_structure_mode_%d.cif" % (imode + 1)))

        readme_string = """

Use Harold Stokes' code, [ISODISTORT](https://stokes.byu.edu/iso/isodistort.php),
loading in your structure that you did the DFPT calculation as the **parent**,
then, select mode decompositional analysis and upload the cif file from step (3).

Follow the on screen instructions.
You will then be presented with the mode irrep and other important symmetry information.

Thanks to Jack Baker for pointing out this approach.
See also <https://forum.abinit.org/viewtopic.php?f=10&t=545>
"""
        with open(os.path.join(workdir, "README.md"), "wt") as fh:
            fh.write(readme_string)

        return workdir

    def decorate_ax(self, ax, units: str = 'eV', **kwargs) -> None:
        """
        Add q-labels, title and unit name to axis ax.
        Use units="" to add k-labels without unit name.

        Args:
            title:
            fontsize
            qlabels:
            qlabel_size:
        """
        title = kwargs.pop("title", None)
        fontsize = kwargs.pop("fontsize", 12)
        if title is not None: ax.set_title(title, fontsize=fontsize)
        ax.grid(True)

        # Handle conversion factor.
        if units:
            ax.set_ylabel(abu.wlabel_from_units(units))

        ax.set_xlabel("Wave Vector")

        # Set ticks and labels.
        ticks, labels = self._make_ticks_and_labels(kwargs.pop("qlabels", None))
        if ticks:
            # Don't show label if previous k-point is the same.
            for il in range(1, len(labels)):
                if labels[il] == labels[il-1]: labels[il] = ""
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.pop("qlabel_size", "large"))
            #print("ticks", len(ticks), ticks)
            ax.set_xlim(ticks[0], ticks[-1])

    def decorate_plotly(self, fig, units: str = 'eV', **kwargs) -> None:
        """
        Add q-labels and unit name to figure ``fig``.
        Use units="" to add k-labels without unit name.

        Args:
            qlabels:
            qlabel_size:
            iax: An int, use iax=n to decorate the nth axis when the fig has subplots.
        """
        iax = kwargs.pop("iax", 1)
        xaxis = 'xaxis%u' % iax

        # Handle conversion factor.
        if units:
            fig.layout['yaxis%u' % iax].title.text = abu.wlabel_from_units(units, unicode=True)

        fig.layout[xaxis].title.text = "Wave Vector"

        # Set ticks and labels.
        ticks, labels = self._make_ticks_and_labels(kwargs.pop("qlabels", None))
        if ticks:
            labels = plotly_klabels(labels)

            fig.layout[xaxis].tickvals = ticks
            fig.layout[xaxis].ticktext = labels
            fig.layout[xaxis].tickfont.size = kwargs.pop("qlabel_size", 16)
            # print("ticks", len(ticks), ticks)
            fig.layout[xaxis].range = (ticks[0], ticks[-1])

    @add_fig_kwargs
    def plot(self, ax=None, units="eV", qlabels=None, branch_range=None, match_bands=False, temp=None,
             fontsize=12, **kwargs):
        r"""
        Plot the phonon band structure with matplotlib.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            branch_range: Tuple specifying the minimum and maximum branch index to plot (default: all branches are plotted).
            match_bands: if True the bands will be matched based on the scalar product between the eigenvectors.
            temp: Temperature in Kelvin. If not None, a scatter plot with the Bose-Einstein occupation factor
                at temperature `temp` is added.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        # Select the band range.
        branch_range = range(self.num_branches) if branch_range is None else \
                       range(branch_range[0], branch_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Decorate the axis (e.g. add ticks and labels).
        self.decorate_ax(ax, units=units, qlabels=qlabels)

        if "color" not in kwargs: kwargs["color"] = "black"
        if "linewidth" not in kwargs: kwargs["linewidth"] = 2.0

        # Plot the phonon branches.
        self.plot_ax(ax, branch_range, units=units, match_bands=match_bands, **kwargs)

        if temp is not None:
            # Scatter plot with Bose-Einstein occupation factors for T = temp
            factor = abu.phfactor_ev2units(units)
            if temp < 1: temp = 1
            ax.set_title("T = %.1f K" % temp, fontsize=fontsize)
            xs = np.arange(self.num_qpoints)
            for nu in self.branches:
                ws = self.phfreqs[:, nu]
                wkt = self.phfreqs[:, nu] / (abu.kb_eVK * temp)
                # 1 / (np.exp(1e-6) - 1)) ~ 999999.5
                wkt = np.where(wkt > 1e-6, wkt, 1e-6)
                occ = 1.0 / (np.exp(wkt) - 1.0)
                s = np.where(occ < 2, occ, 2) * 50
                ax.scatter(xs, ws * factor, s=s, marker="o", c="b", alpha=0.6)
                #ax.scatter(xs, ws, s=s, marker="o", c=occ, cmap="jet")

        return fig

    @add_plotly_fig_kwargs
    def plotly(self, units="eV", qlabels=None, branch_range=None, match_bands=False, temp=None,
               fig=None, rcd=None, fontsize=12, **kwargs):
        r"""
        Plot the phonon band structure with plotly.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            branch_range: Tuple specifying the minimum and maximum branch index to plot (default: all branches are plotted).
            match_bands: if True the bands will be matched based on the scalar product between the eigenvectors.
            temp: Temperature in Kelvin. If not None, a scatter plot with the Bose-Einstein occupation factor
                at temperature `temp` is added.
            fig: plotly figure or None if a new figure should be created.
            rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col) of the subplot in the grid.
            fontsize: Title fontsize.

        Returns: |plotly.graph_objects.Figure|
        """
        # Select the band range.
        branch_range = range(self.num_branches) if branch_range is None else \
                       range(branch_range[0], branch_range[1], 1)

        fig, _ = get_fig_plotly(fig=fig)

        # Decorate the axis (e.g. add ticks and labels).
        rcd = PlotlyRowColDesc.from_object(rcd)
        self.decorate_plotly(fig, units=units, qlabels=qlabels, iax=rcd.iax)

        if "color" not in kwargs: kwargs["color"] = "black"
        if "linewidth" not in kwargs: kwargs["linewidth"] = 2.0

        # Plot the phonon branches.
        self.plotly_traces(fig, branch_range, rcd=rcd, units=units, match_bands=match_bands, **kwargs)

        if temp is not None:
            # Scatter plot with Bose-Einstein occupation factors for T = temp
            factor = abu.phfactor_ev2units(units)
            if temp < 1: temp = 1
            fig.layout.annotations = [dict(text="T = %.1f K" % temp, font_size=fontsize, x=0.5, xref='paper',
                                           xanchor='center', y=1, yref='paper', yanchor='bottom', showarrow=False)]
            xs = np.arange(self.num_qpoints)
            for nu in self.branches:
                ws = self.phfreqs[:, nu]
                wkt = self.phfreqs[:, nu] / (abu.kb_eVK * temp)
                # 1 / (np.exp(1e-6) - 1)) ~ 999999.5
                wkt = np.where(wkt > 1e-6, wkt, 1e-6)
                occ = 1.0 / (np.exp(wkt) - 1.0)
                s = np.where(occ < 0.3, occ, 0.3) * 50
                #print("rcd", rcd)
                fig.add_scatter(x=xs, y=ws * factor, mode='markers', row=rcd.ply_row, col=rcd.ply_col, showlegend=False,
                                marker=dict(color='blue', size=s, opacity=0.6, line_width=0), name='')
                #               marker=dict(color=occ, colorscale='jet', size=s, opacity=0.6, line_width=0),
        return fig

    def plot_ax(self, ax, branch, units='eV', match_bands=False, **kwargs):
        """
        Plots the frequencies for the given branches indices as a function of the q-index on axis ``ax``.
        If ``branch`` is None, all phonon branches are plotted.

        Return: The list of matplotlib lines added.
        """
        if branch is None:
            branch_range = range(self.num_branches)
        elif isinstance(branch, (list, tuple, np.ndarray)):
            branch_range = branch
        else:
            branch_range = [branch]

        first_xx = 0
        lines = []

        factor = abu.phfactor_ev2units(units)

        for i, pf in enumerate(self.split_phfreqs):
            if match_bands:
                ind = self.split_matched_indices[i]
                pf = pf[np.arange(len(pf))[:, None], ind]
            pf = pf * factor
            xx = list(range(first_xx, first_xx + len(pf)))
            for branch in branch_range:
                lines.extend(ax.plot(xx, pf[:, branch], **kwargs))
            first_xx = xx[-1]

        return lines

    def plotly_traces(self, fig, branch, rcd=None, units='eV', name='', match_bands=False,
                      showlegend=False, **kwargs):
        """
        Plots the frequencies for the given branches indices as a function of the q-index on figure ``fig`` .
        If ``fig`` has subplots, ``rcd`` is used to add traces on these subplots.
        If ``branch`` is None, all phonon branches are plotted.
        kwargs: Passed to fig.add_scatter method.
        """
        linecolor = kwargs.pop("color", "black")
        linewidth = kwargs.pop("linewidth", 2.0)

        rcd = PlotlyRowColDesc.from_object(rcd)
        ply_row, ply_col = rcd.ply_row, rcd.ply_col

        if branch is None:
            branch_range = range(self.num_branches)
        elif isinstance(branch, (list, tuple, np.ndarray, range)):
            branch_range = branch
        else:
            branch_range = [branch]

        first_xx = 0
        factor = abu.phfactor_ev2units(units)

        for i, pf in enumerate(self.split_phfreqs):
            if match_bands:
                ind = self.split_matched_indices[i]
                pf = pf[np.arange(len(pf))[:, None], ind]
            pf = pf * factor
            xx = list(range(first_xx, first_xx + len(pf)))
            for branch in branch_range:
                fig.add_scatter(x=xx, y=pf[:, branch], mode='lines', name=name, legendgroup=name, showlegend=False,
                                   line=dict(color=linecolor, width=linewidth), **kwargs, row=ply_row, col=ply_col)
            first_xx = xx[-1]

        if showlegend:
            fig.data[-1].showlegend = True

    @add_fig_kwargs
    def plot_colored_matched(self, ax=None, units="eV", qlabels=None, branch_range=None,
                             colormap="rainbow", max_colors=None, **kwargs):
        r"""
        Plot the phonon band structure with different colors for each line.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            branch_range: Tuple specifying the minimum and maximum branch_i index to plot
                (default: all branches are plotted).
            colormap: matplotlib colormap to determine the colors available. The colors will be chosen not in a
                sequential order to avoid difficulties in distinguishing the lines.
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            max_colors: maximum number of colors to be used. If max_colors < num_braches the colors will be reapeated.
                It may be useful to better distinguish close bands when the number of branches is large.

        Returns: |matplotlib-Figure|
        """
        # Select the band range.
        if branch_range is None:
            branch_range = range(self.num_branches)
        else:
            branch_range = range(branch_range[0], branch_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, units=units, qlabels=qlabels)

        first_xx = 0
        lines = []
        factor = abu.phfactor_ev2units(units)

        if max_colors is None:
            max_colors = len(branch_range)

        colormap = plt.get_cmap(colormap)

        for i, pf in enumerate(self.split_phfreqs):
            ind = self.split_matched_indices[i]
            pf = pf[np.arange(len(pf))[:, None], ind]
            pf = pf * factor
            xx = range(first_xx, first_xx + len(pf))
            colors = itertools.cycle(colormap(np.linspace(0, 1, max_colors)))
            for branch_i in branch_range:
                kwargs = dict(kwargs)
                kwargs['color'] = next(colors)
                lines.extend(ax.plot(xx, pf[:, branch_i], **kwargs))
            first_xx = xx[-1]

        return fig

    @add_fig_kwargs
    def plot_lt_character(self, units="eV", qlabels=None, ax=None, xlims=None, ylims=None, scale_size=50,
                          use_becs=True, colormap="jet", fontsize=12, **kwargs):
        r"""
        Plot the phonon band structure with colored lines. The color of the lines indicates
        the degree to which the mode is longitudinal.
        Red corresponds to longitudinal modes and black to purely transverse modes.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: y-axis limits.
            scale_size: Scaling factor for marker size. Increase this value to enlarge the markers.
            use_becs: True to compute polar strength: q . Z[atom] . disp[q, nu, atom]
                False to use: q . disp[q, nu, atom]. Useful to highlight LA/TA modes.
            colormap: Matplotlib colormap.
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        if use_becs and self.zcart is None:
            cprint("Bandstructure does not have Born effective charges", "yellow")
            return None

        factor = abu.phfactor_ev2units(units)
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)

        if "color" not in kwargs: kwargs["color"] = "black"
        if "linewidth" not in kwargs: kwargs["linewidth"] = 2.0

        first_xx = 0
        scatt_x, scatt_y, scatt_s = [], [], []
        for p_qpts, p_freqs, p_dcart in zip(self.split_qpoints, self.split_phfreqs, self.split_phdispl_cart):
            xx = list(range(first_xx, first_xx + len(p_freqs)))

            for iq, (qpt, ws, dis) in enumerate(zip(p_qpts, p_freqs, p_dcart)):
                qcart = self.structure.reciprocal_lattice.get_cartesian_coords(qpt)
                qnorm = np.linalg.norm(qcart)
                inv_qepsq = 0.0
                if qnorm > 1e-3:
                    qvers = qcart / qnorm
                    inv_qepsq = 1.0 / np.dot(qvers, np.dot(self.epsinf, qvers))

                # We are not interested in the amplitudes so normalize all displacements to one.
                dis = dis.reshape(self.num_branches, self.num_atoms, 3)
                for nu in range(self.num_branches):
                    if use_becs:
                        # q . Z[atom] . disp[q, nu, atom]
                        v = sum(np.dot(qcart, np.dot(self.zcart[iatom], dis[nu, iatom])) for iatom in range(self.num_atoms))
                    else:
                        v = sum(np.dot(qcart, dis[nu, iatom]) for iatom in range(self.num_atoms))
                    scatt_x.append(xx[iq])
                    scatt_y.append(ws[nu])
                    scatt_s.append(v * inv_qepsq)

            p_freqs = p_freqs * factor
            ax.plot(xx, p_freqs, **kwargs)
            first_xx = xx[-1]

        scatt_y = np.array(scatt_y) * factor
        scatt_s = np.abs(np.array(scatt_s))
        scatt_s /= scatt_s.max()
        scatt_s *= scale_size
        #print("scatt_s", scatt_s, "min", scatt_s.min(), "max", scatt_s.max())

        ax.scatter(scatt_x, scatt_y, s=scatt_s,
            #c=None, marker=None, cmap=None, norm=None, vmin=None, vmax=None, alpha=None,
            #linewidths=None, verts=None, edgecolors=None, *, data=None
        )
        self.decorate_ax(ax, units=units, qlabels=qlabels)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")

        return fig

    @property
    def split_qpoints(self):
        try:
            return self._split_qpoints
        except AttributeError:
            self._set_split_intervals()
            return self._split_qpoints

    @property
    def split_phfreqs(self):
        try:
            return self._split_phfreqs
        except AttributeError:
            self._set_split_intervals()
            return self._split_phfreqs

    @property
    def split_phdispl_cart(self):
        # prepare the splitted phdispl_cart as a separate internal variable only when explicitely requested and
        # not at the same time as split_qpoints and split_phfreqs as it requires a larger array and not used
        # most of the times.
        try:
            return self._split_phdispl_cart
        except AttributeError:
            self.split_phfreqs
            split_phdispl_cart = [np.array(self.phdispl_cart[self._split_indices[i]:self._split_indices[i + 1] + 1])
                                  for i in range(len(self._split_indices) - 1)]
            if self.non_anal_ph is not None:
                for i, q in enumerate(self.split_qpoints):
                    if np.array_equal(q[0], (0, 0, 0)):
                        if self.non_anal_ph.has_direction(q[1]):
                            split_phdispl_cart[i][0, :] = self._get_non_anal_phdispl(q[1])
                    if np.array_equal(q[-1], (0, 0, 0)):
                        if self.non_anal_ph.has_direction(q[-2]):
                            split_phdispl_cart[i][-1, :] = self._get_non_anal_phdispl(q[-2])

            self._split_phdispl_cart = split_phdispl_cart
            return self._split_phdispl_cart

    def _set_split_intervals(self):
        # Calculations available for LO-TO splitting
        # Split the lines at each Gamma to handle possible discontinuities
        if self.non_anal_phfreqs is not None and self.non_anal_directions is not None:
            end_points_indices = [0]

            end_points_indices.extend(
                [i for i in range(1, self.num_qpoints - 1) if np.array_equal(self.qpoints.frac_coords[i], [0, 0, 0])])
            end_points_indices.append(self.num_qpoints - 1)

            # split the list of qpoints and frequencies at each end point. The end points are in both the segments.
            # Lists since the array contained have different shapes
            split_qpoints = [np.array(self.qpoints.frac_coords[end_points_indices[i]:end_points_indices[i + 1] + 1])
                             for i in range(len(end_points_indices) - 1)]
            split_phfreqs = [np.array(self.phfreqs[end_points_indices[i]:end_points_indices[i + 1] + 1])
                             for i in range(len(end_points_indices) - 1)]

            for i, q in enumerate(split_qpoints):
                if np.array_equal(q[0], (0, 0, 0)):
                    split_phfreqs[i][0, :] = self._get_non_anal_freqs(q[1])
                if np.array_equal(q[-1], (0, 0, 0)):
                    split_phfreqs[i][-1, :] = self._get_non_anal_freqs(q[-2])
        else:
            split_qpoints = [self.qpoints.frac_coords]
            split_phfreqs = [self.phfreqs]
            end_points_indices = [0, self.num_qpoints-1]

        self._split_qpoints = split_qpoints
        self._split_phfreqs = split_phfreqs
        self._split_indices = end_points_indices
        return split_phfreqs, split_qpoints

    @property
    def split_matched_indices(self):
        """
        A list of numpy arrays containing the indices in which each band should be sorted in order to match the
        scalar product of the eigenvectors. The shape is the same as that of split_phfreqs.
        Lazy property.
        """
        try:
            return self._split_matched_indices
        except AttributeError:

            split_matched_indices = []
            last_eigenvectors = None

            # simpler method based just on the matching with the previous point
            #TODO remove after verifying the other method currently in use
            # for i, displ in enumerate(self.split_phdispl_cart):
            #     eigenvectors = get_dyn_mat_eigenvec(displ, self.structure, amu=self.amu)
            #     ind_block = np.zeros((len(displ), self.num_branches), dtype=int)
            #     # if it's not the first block, match with the last of the previous block. Should give a match in case
            #     # of LO-TO splitting
            #     if i == 0:
            #         ind_block[0] = range(self.num_branches)
            #     else:
            #         match = match_eigenvectors(last_eigenvectors, eigenvectors[0])
            #         ind_block[0] = [match[m] for m in split_matched_indices[-1][-1]]
            #     for j in range(1, len(displ)):
            #         match = match_eigenvectors(eigenvectors[j-1], eigenvectors[j])
            #         ind_block[j] = [match[m] for m in ind_block[j-1]]
            #
            #     split_matched_indices.append(ind_block)
            #     last_eigenvectors = eigenvectors[-1]

            # The match is applied between subsequent qpoints, except that right after a high symmetry point.
            # In that case the first point after the high symmetry point will be matched with the one immediately
            # before. This should avoid exchange of lines due to degeneracies.
            # The code will assume that there is a high symmetry point if the points are not collinear (change in the
            # direction in the path).
            def collinear(a, b, c):
                v1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]]
                v2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]]
                d = [v1, v2, [1, 1, 1]]
                return np.isclose(np.linalg.det(d), 0, atol=1e-5)

            for i, displ in enumerate(self.split_phdispl_cart):
                eigenvectors = get_dyn_mat_eigenvec(displ, self.structure, amu=self.amu)
                ind_block = np.zeros((len(displ), self.num_branches), dtype=int)
                # if it's not the first block, match the first two points with the last of the previous block.
                # Should give a match in case of LO-TO splitting
                if i == 0:
                    ind_block[0] = range(self.num_branches)
                    match = match_eigenvectors(eigenvectors[0], eigenvectors[1])
                    ind_block[1] = [match[m] for m in ind_block[0]]
                else:
                    match = match_eigenvectors(last_eigenvectors, eigenvectors[0])
                    ind_block[0] = [match[m] for m in split_matched_indices[-1][-2]]
                    match = match_eigenvectors(last_eigenvectors, eigenvectors[1])
                    ind_block[1] = [match[m] for m in split_matched_indices[-1][-2]]
                for j in range(2, len(displ)):
                    k = j-1
                    if not collinear(self.split_qpoints[i][j-2], self.split_qpoints[i][j-1], self.split_qpoints[i][j]):
                        k = j-2
                    match = match_eigenvectors(eigenvectors[k], eigenvectors[j])
                    ind_block[j] = [match[m] for m in ind_block[k]]

                split_matched_indices.append(ind_block)
                last_eigenvectors = eigenvectors[-2]

            self._split_matched_indices = split_matched_indices

            return self._split_matched_indices

    def _get_non_anal_freqs(self, frac_direction):
        # directions for the qph2l in anaddb are given in cartesian coordinates
        cart_direction = self.structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(frac_direction)
        cart_direction = cart_direction / np.linalg.norm(cart_direction)

        for i, d in enumerate(self.non_anal_directions):
            d = d / np.linalg.norm(d)
            if np.allclose(cart_direction, d):
                return self.non_anal_phfreqs[i]

        raise ValueError("Non analytical contribution has not been calculated for reduced direction {0} ".format(frac_direction))

    def _get_non_anal_phdispl(self, frac_direction):
        # directions for the qph2l in anaddb are given in cartesian coordinates
        cart_direction = self.structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(frac_direction)
        cart_direction = cart_direction / np.linalg.norm(cart_direction)

        for i, d in enumerate(self.non_anal_directions):
            d = d / np.linalg.norm(d)
            if np.allclose(cart_direction, d):
                return self.non_anal_phdispl_cart[i]

        raise ValueError("Non analytical contribution has not been calcolated for reduced direction {0} ".format(frac_direction))

    def _make_ticks_and_labels(self, qlabels):
        """Return ticks and labels from the mapping {qred: qstring} given in qlabels."""
        # TODO should be modified in order to handle the "split" list of qpoints
        if qlabels is not None:
            d = OrderedDict()

            for qcoord, qname in qlabels.items():
                # Build Kpoint instancee
                qtick = Kpoint(qcoord, self.structure.reciprocal_lattice)
                for q, qpoint in enumerate(self.qpoints):
                    if qtick == qpoint:
                        d[q] = qname
        else:
            d = self._auto_qlabels

        # Return ticks, labels
        return list(d.keys()), list(d.values())

    # TODO: fatbands along x, y, z
    @add_fig_kwargs
    def plot_fatbands(self, use_eigvec=True, units="eV", colormap="jet", phdos_file=None,
                      alpha=0.6, max_stripe_width_mev=5.0, width_ratios=(2, 1),
                      qlabels=None, ylims=None, fontsize=12, **kwargs):
        r"""
        Plot phonon fatbands and, optionally, atom-projected phonon DOSes with matplotlib.
        The width of the band is given by ||v_{type}||
        where v is the (complex) phonon displacement (eigenvector) in cartesian coordinates and
        v_{type} selects only the terms associated to the atomic type.

        Args:
            use_eigvec: True if the width of the phonon branch should be computed from the eigenvectors.
                False to use phonon displacements. Note that the PHDOS is always decomposed in
                terms of (orthonormal) eigenvectors.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            colormap: Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            phdos_file: Used to activate fatbands + PJDOS plot.
                Accept string with path of PHDOS.nc file or |PhdosFile| object.
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            max_stripe_width_mev: The maximum width of the stripe in meV. Will be rescaled according to ``units``.
            width_ratios: Ratio between the width of the fatbands plots and the DOS plots.
                Used if `phdos_file` is not None
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        lw = kwargs.pop("lw", 2)
        factor = abu.phfactor_ev2units(units)
        ntypat = self.structure.ntypesp

        # Prepare PJDOS.
        close_phdos_file = False
        if phdos_file is not None:
            if is_string(phdos_file):
                phdos_file = PhdosFile(phdos_file)
                close_phdos_file = True
            else:
                if not isinstance(phdos_file, PhdosFile):
                    raise TypeError("Expecting string or PhdosFile, got %s" % type(phdos_file))

        # Grid with [ntypat] plots if fatbands only or [ntypat, 2] if fatbands + PJDOS
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        fig = plt.figure()
        nrows, ncols = (ntypat, 1) if phdos_file is None else (ntypat, 2)
        gspec = GridSpec(nrows=nrows, ncols=ncols, width_ratios=width_ratios if ncols == 2 else None,
                         wspace=0.05, hspace=0.1)

        cmap = plt.get_cmap(colormap)
        qq = list(range(self.num_qpoints))

        # phonon_displacements are in cartesian coordinates and stored in an array with shape
        # (nqpt, 3*natom, 3*natom) where the last dimension stores the cartesian components.
        # PJDoses are in cartesian coordinates and are computed by anaddb using the the
        # phonon eigenvectors that are orthonormal.

        # Precompute normalization factor:
        # Here I use d2[q, nu] = \sum_{i=0}^{3*Nat-1) |d^{q\nu}_i|**2
        # it makes sense only for displacements
        d2_qnu = np.ones((self.num_qpoints, self.num_branches))
        if not use_eigvec:
            for iq in range(self.num_qpoints):
                for nu in self.branches:
                    cvect = self.phdispl_cart[iq, nu, :]
                    d2_qnu[iq, nu] = np.vdot(cvect, cvect).real

        # Plot fatbands: one plot per atom type.
        ax00 = None
        for ax_row, symbol in enumerate(self.structure.symbol_set):
            last_ax = (ax_row == len(self.structure.symbol_set) - 1)
            ax = plt.subplot(gspec[ax_row, 0], sharex=ax00, sharey=ax00)
            if ax_row == 0: ax00 = ax
            self.decorate_ax(ax, units=units, qlabels=qlabels)
            color = cmap(float(ax_row) / max(1, ntypat - 1))

            # dir_indices lists the coordinate indices for the atoms of the same type.
            atom_indices = self.structure.indices_from_symbol(symbol)
            dir_indices = []
            for aindx in atom_indices:
                start = 3 * aindx
                dir_indices.extend([start, start + 1, start + 2])
            dir_indices = np.array(dir_indices)

            for nu in self.branches:
                yy_qq = self.phfreqs[:, nu] * factor

                # Exctract the sub-vector associated to this atom type (eigvec or diplacement).
                if use_eigvec:
                    v_type = self.dyn_mat_eigenvect[:, nu, dir_indices]
                else:
                    v_type = self.phdispl_cart[:, nu, dir_indices]

                v2_type = np.empty(self.num_qpoints)
                for iq in range(self.num_qpoints):
                    v2_type[iq] = np.vdot(v_type[iq], v_type[iq]).real

                # Normalize and scale by max_stripe_width_mev taking into account units.
                # The stripe is centered on the phonon branch hence the factor 2
                stype_qq = (factor * max_stripe_width_mev * 1.e-3 / 2) * np.sqrt(v2_type / d2_qnu[:, nu])

                # Plot the phonon branch with the stripe.
                if nu == 0:
                    ax.plot(qq, yy_qq, lw=lw, color=color, label=symbol)
                else:
                    ax.plot(qq, yy_qq, lw=lw, color=color)

                ax.fill_between(qq, yy_qq + stype_qq, yy_qq - stype_qq, facecolor=color, alpha=alpha, linewidth=0)

            set_axlims(ax, ylims, "y")
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        # Type projected DOSes (always computed from eigenvectors in anaddb).
        if phdos_file is not None:
            ax01 = None
            for ax_row, symbol in enumerate(self.structure.symbol_set):
                color = cmap(float(ax_row) / max(1, ntypat - 1))
                ax = plt.subplot(gspec[ax_row, 1], sharex=ax01, sharey=ax00)
                if ax_row == 0: ax01 = ax

                # Get PJDOS: Dictionary symbol --> partial PhononDos
                pjdos = phdos_file.pjdos_symbol[symbol]
                x, y = pjdos.mesh * factor, pjdos.values / factor

                ax.plot(y, x, lw=lw, color=color)
                ax.grid(True)
                ax.yaxis.set_ticks_position("right")
                ax.yaxis.set_label_position("right")
                set_axlims(ax, ylims, "y")

            if close_phdos_file:
                phdos_file.close()

        return fig

    # TODO: fatbands along x, y, z
    @add_plotly_fig_kwargs
    def plotly_fatbands(self, use_eigvec=True, units="eV", colormap="G10", phdos_file=None,
                        alpha=0.6, max_stripe_width_mev=5.0, width_ratios=(2, 1),
                        qlabels=None, ylims=None, fontsize=16, **kwargs):
        r"""
        Plot phonon fatbands and, optionally, atom-projected phonon DOSes with plotly.
        The width of the band is given by ||v_{type}||
        where v is the (complex) phonon displacement (eigenvector) in cartesian coordinates and
        v_{type} selects only the terms associated to the atomic type.

        Args:
            use_eigvec: True if the width of the phonon branch should be computed from the eigenvectors.
                False to use phonon displacements. Note that the PHDOS is always decomposed in
                terms of (orthonormal) eigenvectors.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            colormap: Have a look at the colormaps here and decide which one you like:
                https://plotly.com/python/discrete-color/
            phdos_file: Used to activate fatbands + PJDOS plot.
                Accept string with path of PHDOS.nc file or |PhdosFile| object.
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            max_stripe_width_mev: The maximum width of the stripe in meV. Will be rescaled according to ``units``.
            width_ratios: Ratio between the width of the fatbands plots and the DOS plots.
                Used if `phdos_file` is not None
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            fontsize: Legend and title fontsize.

        Returns: |plotly.graph_objects.Figure|
        """
        lw = kwargs.pop("lw", 2)
        factor = abu.phfactor_ev2units(units)
        ntypat = self.structure.ntypesp

        # Prepare PJDOS.
        close_phdos_file = False
        if phdos_file is not None:
            if is_string(phdos_file):
                phdos_file = PhdosFile(phdos_file)
                close_phdos_file = True
            else:
                if not isinstance(phdos_file, PhdosFile):
                    raise TypeError("Expecting string or PhdosFile, got %s" % type(phdos_file))

        # Fig with [ntypat] plots if fatbands only or [ntypat, 2] if fatbands + PJDOS
        nrows, ncols = (ntypat, 1) if phdos_file is None else (ntypat, 2)
        fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, sharex=True, sharey=True, vertical_spacing=0.05,
                                 horizontal_spacing=0.02, column_widths=width_ratios if ncols == 2 else None)

        import plotly.express as px
        color_l = getattr(px.colors.qualitative, colormap)
        if len(color_l) < len(self.structure.symbol_set):
            raise ValueError(f"Colormap {colormap} is not enough, please provide more than %d colors"
                             % len(self.structure.symbol_set))
        qq = list(range(self.num_qpoints))

        # phonon_displacements are in cartesian coordinates and stored in an array with shape
        # (nqpt, 3*natom, 3*natom) where the last dimension stores the cartesian components.
        # PJDoses are in cartesian coordinates and are computed by anaddb using the the
        # phonon eigenvectors that are orthonormal.

        # Precompute normalization factor:
        # Here I use d2[q, nu] = \sum_{i=0}^{3*Nat-1) |d^{q\nu}_i|**2
        # it makes sense only for displacements
        d2_qnu = np.ones((self.num_qpoints, self.num_branches))
        if not use_eigvec:
            for iq in range(self.num_qpoints):
                for nu in self.branches:
                    cvect = self.phdispl_cart[iq, nu, :]
                    d2_qnu[iq, nu] = np.vdot(cvect, cvect).real

        # Plot fatbands: one plot per atom type.
        for row, symbol in enumerate(self.structure.symbol_set):
            color = color_l[row]
            rcd = PlotlyRowColDesc(row, 0, nrows, ncols)
            iax=rcd.iax
            self.decorate_plotly(fig, units=units, qlabels=qlabels, iax=iax)
            if row != len(self.structure.symbol_set):
                xaxis = 'xaxis%u' % iax
                fig.layout[xaxis].title.text = ""

            # dir_indices lists the coordinate indices for the atoms of the same type.
            atom_indices = self.structure.indices_from_symbol(symbol)
            dir_indices = []
            for aindx in atom_indices:
                start = 3 * aindx
                dir_indices.extend([start, start + 1, start + 2])
            dir_indices = np.array(dir_indices)

            for nu in self.branches:
                yy_qq = self.phfreqs[:, nu] * factor

                # Ectract the sub-vector associated to this atom type (eigvec or displacement).
                if use_eigvec:
                    v_type = self.dyn_mat_eigenvect[:, nu, dir_indices]
                else:
                    v_type = self.phdispl_cart[:, nu, dir_indices]

                v2_type = np.empty(self.num_qpoints)
                for iq in range(self.num_qpoints):
                    v2_type[iq] = np.vdot(v_type[iq], v_type[iq]).real

                # Normalize and scale by max_stripe_width_mev taking into account units.
                # The stripe is centered on the phonon branch hence the factor 2
                stype_qq = (factor * max_stripe_width_mev * 1.e-3 / 2) * np.sqrt(v2_type / d2_qnu[:, nu])

                # Plot the phonon branch with the stripe.
                ply_row, ply_col = rcd.ply_row, rcd.ply_col
                if nu == 0:
                    fig.add_scatter(x=qq, y=yy_qq, mode='lines', line=dict(width=lw, color=color), name=symbol,
                                    legendgroup=row ,row=ply_row, col=ply_col)
                else:
                    fig.add_scatter(x=qq, y=yy_qq, mode='lines', line=dict(width=lw, color=color), name=symbol,
                                    showlegend=False, legendgroup=row ,row=ply_row, col=ply_col)

                fig.add_scatter(x=qq, y=yy_qq-stype_qq, mode='lines', line=dict(width=0, color=color), name='',
                                showlegend=False, legendgroup=row ,row=ply_row, col=ply_col)
                fig.add_scatter(x=qq, y=yy_qq+stype_qq, mode='lines', line=dict(width=0, color=color), name='',
                                showlegend=False, legendgroup=row, fill='tonexty', row=ply_row, col=ply_col)

            plotly_set_lims(fig, ylims, "y", iax=iax)

        # Type projected DOSes (always computed from eigenvectors in anaddb).
        if phdos_file is not None:
            for row, symbol in enumerate(self.structure.symbol_set):
                color = color_l[row]
                rcd = PlotlyRowColDesc(row, 1, nrows, ncols)

                # Get PJDOS: Dictionary symbol --> partial PhononDos
                pjdos = phdos_file.pjdos_symbol[symbol]
                x, y = pjdos.mesh * factor, pjdos.values / factor

                ply_row, ply_col = rcd.ply_row, rcd.ply_col
                fig.add_scatter(x=y, y=x, mode='lines', line=dict(width=lw, color=color), name='', showlegend=False,
                                legendgroup=row ,row=ply_row, col=ply_col)
                plotly_set_lims(fig, ylims, "y", iax=rcd.iax)

            if close_phdos_file:
                phdos_file.close()

        fig.layout.legend.font.size = fontsize
        fig.layout.title.font.size = fontsize + 2

        return fig

    @add_fig_kwargs
    def plot_with_phdos(self, phdos, units="eV", qlabels=None, ax_list=None, width_ratios=(2, 1), **kwargs):
        r"""
        Plot the phonon band structure with the phonon DOS.

        Args:
            phdos: An instance of |PhononDos| or a netcdf file providing a PhononDos object.
            units: Units for plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            ax_list: The axes for the bandstructure plot and the DOS plot. If ax_list is None, a new figure
                is created and the two axes are automatically generated.
            width_ratios: Ratio between the width of the bands plots and the DOS plots.
                Used if ``ax_list`` is None

        Returns: |matplotlib-Figure|
        """
        phdos = PhononDos.as_phdos(phdos, phdos_kwargs=None)

        import matplotlib.pyplot as plt
        if ax_list is None:
            # Build axes and align bands and DOS.
            from matplotlib.gridspec import GridSpec
            fig = plt.figure()
            gspec = GridSpec(1, 2, width_ratios=width_ratios, wspace=0.05)
            ax1 = plt.subplot(gspec[0])
            ax2 = plt.subplot(gspec[1], sharey=ax1)
        else:
            # Take them from ax_list.
            ax1, ax2 = ax_list
            fig = plt.gcf()

        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the phonon band structure.
        self.plot_ax(ax1, branch=None, units=units, **kwargs)
        self.decorate_ax(ax1, units=units, qlabels=qlabels)

        factor = abu.phfactor_ev2units(units)
        emin = np.min(self.minfreq)
        emin -= 0.05 * abs(emin)
        emin *= factor
        emax = np.max(self.maxfreq)
        emax += 0.05 * abs(emax)
        emax *= factor
        ax1.set_ylim(emin, emax)

        # Plot Phonon DOS
        phdos.plot_dos_idos(ax2, what="d", units=units, exchange_xy=True, **kwargs)

        ax2.grid(True)
        ax2.yaxis.set_ticks_position("right")
        #ax2.yaxis.set_label_position("right")

        return fig

    @add_plotly_fig_kwargs
    def plotly_with_phdos(self, phdos, units="eV", qlabels=None, fig=None, rcd_phbands=None, rcd_phdos=None,
                          width_ratios=(2, 1), fontsize=12, **kwargs):
        r"""
        Plot the phonon band structure with the phonon DOS with plotly.

        Args:
            phdos: An instance of |PhononDos| or a netcdf file providing a PhononDos object.
            units: Units for plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            fig: The fig for the bandstructure plot and the DOS plot. If fig is None, a new figure
                is created.
            width_ratios: Ratio between the width of the bands plots and the DOS plots.
                Used if ``fig`` is None
            fontsize: Title fontsize.

        Returns: |plotly.graph_objects.Figure|
        """
        phdos = PhononDos.as_phdos(phdos, phdos_kwargs=None)

        if fig is None:
            # build fig and align bands and DOS.
            fig, _ = get_figs_plotly(nrows=1, ncols=2, subplot_titles=[], sharex=False, sharey=True,
                                     horizontal_spacing=0.02, column_widths=width_ratios)

        if not kwargs:
            kwargs = {"line_color": "black", "line_width": 2.0}

        # Plot the phonon band structure.
        rcd_phbands = PlotlyRowColDesc.from_object(rcd_phbands)
        self.plotly_traces(fig, branch=None, rcd=rcd_phbands, units=units, **kwargs)
        self.decorate_plotly(fig, units=units, qlabels=qlabels, iax=rcd_phbands.iax)

        factor = abu.phfactor_ev2units(units)
        emin = np.min(self.minfreq)
        emin -= 0.05 * abs(emin)
        emin *= factor
        emax = np.max(self.maxfreq)
        emax += 0.05 * abs(emax)
        emax *= factor
        fig.layout.yaxis.range = (emin, emax)

        # Plot Phonon DOS
        if rcd_phdos is not None:
            rcd_phdos = PlotlyRowColDesc.from_object(rcd_phdos)
        else:
            rcd_phdos = PlotlyRowColDesc(0, 1, 1, 2)

        phdos.plotly_dos_idos(fig, rcd=rcd_phdos, what="d", units=units, exchange_xy=True, showlegend=False, **kwargs)
        fig.update_layout(font_size=fontsize)

        return fig

    @add_fig_kwargs
    def plot_phdispl(self, qpoint, cart_dir=None, use_reduced_coords=False, ax=None, units="eV",
                     is_non_analytical_direction=False, use_eigvec=False,
                     colormap="viridis", hatches="default", atoms_index=None, labels_groups=None,
                     normalize=True, use_sqrt=False, fontsize=12, branches=None, format_w="%.3f", **kwargs):
        """
        Plot vertical bars with the contribution of the different atoms or atomic types to all the phonon modes
        at a given ``qpoint``. The contribution is given by ||v_{type}||
        where v is the (complex) phonon displacement (eigenvector) in cartesian coordinates and
        v_{type} selects only the terms associated to the atomic type.
        Options allow to specify which atoms should be taken into account and how should be reparted.

        Args:
            qpoint: integer, vector of reduced coordinates or |Kpoint| object.
            cart_dir: "x", "y", or "z" to select a particular Cartesian directions. or combinations separated by "+".
                Example: "x+y". None if no projection is wanted.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon frequencies. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            is_non_analytical_direction: If True, the ``qpoint`` is interpreted as a direction in q-space
                and the phonon (displacements/eigenvectors) for q --> 0 along this direction are used.
                Requires band structure with :class:`NonAnalyticalPh` object.
            use_eigvec: True if eigenvectors should be used instead of displacements (eigenvectors
                are orthonormal, unlike diplacements)
            colormap: Matplotlib colormap used for atom type.
            hatches: List of strings with matplotlib hatching patterns. None or empty list to disable hatching.
            fontsize: Legend and title fontsize.
            normalize: if True divides by the square norm of the total eigenvector
            use_sqrt: if True the square root of the sum of the components will be taken
            use_reduced_coords: if True coordinates will be converted to reduced coordinates. So the values will be
                fraction of a,b,c rather than x,y,z.
            atoms_index: list of lists. Each list contains the indices of atoms in the structure that will be
                summed on a separate group. if None all the atoms will be considered and grouped by type.
            labels_groups: If atoms_index is not None will provide the labels for each of the group in atoms_index.
                Should have the same length of atoms_index or be None. If None automatic labelling will be used.
            branches: list of indices for the modes that should be represented. If None all the modes will be shown.
            format_w: string used to format the values of the frequency. Default "%.3f".

        Returns: |matplotlib-Figure|
        """
        factor = abu.phfactor_ev2units(units)

        dxyz = {"x": 0, "y": 1, "z": 2, None: None}

        if cart_dir is None:
            icart = None
        else:
            icart = [dxyz[c] for c in cart_dir.split("+")]

        iq, qpoint = self.qindex_qpoint(qpoint, is_non_analytical_direction=is_non_analytical_direction)

        if use_sqrt:
            f_sqrt = np.sqrt
        else:
            f_sqrt = lambda x: x

        if branches is None:
            branches = self.branches
        elif not isinstance(branches, (list, tuple)):
            branches = [branches]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)
        ntypat = self.structure.ntypesp

        if is_non_analytical_direction:
            ax.set_title("q-direction = %s" % repr(qpoint), fontsize=fontsize)
        else:
            ax.set_title("qpoint = %s" % repr(qpoint), fontsize=fontsize)
        ax.set_xlabel('Frequency %s' % abu.phunit_tag(units))

        what = r"\epsilon" if use_eigvec else "d"
        if icart is None:
            ax.set_ylabel(r"${|\vec{%s}_{type}|} (stacked)$" % what, fontsize=fontsize)
        else:
            ax.set_ylabel(r"${|\vec{%s}_{%s,type}|} (stacked)$" % (what, cart_dir), fontsize=fontsize)

        symbol2indices = self.structure.get_symbol2indices()

        width, pad = 4, 1
        pad = width + pad
        xticks, xticklabels = [], []
        if hatches == "default":
            hatches = ["/", "\\", "'", "|", "-", "+", "x", "o", "O", ".", "*"]
        else:
            hatches = list_strings(hatches) if hatches is not None else []

        x = 0
        for inu, nu in enumerate(branches):
            # Select frequencies and cartesian displacements/eigenvectors
            if is_non_analytical_direction:
                w_qnu = self.non_anal_phfreqs[iq, nu] * factor
                if use_eigvec:
                    vcart_qnu = np.reshape(self.non_anal_ph.dyn_mat_eigenvect[iq, nu], (len(self.structure), 3))
                else:
                    vcart_qnu = np.reshape(self.non_anal_phdispl_cart[iq, nu], (len(self.structure), 3))
            else:
                w_qnu = self.phfreqs[iq, nu] * factor
                if use_eigvec:
                    vcart_qnu = np.reshape(self.dyn_mat_eigenvect[iq, nu], (len(self.structure), 3))
                else:
                    vcart_qnu = np.reshape(self.phdispl_cart[iq, nu], (len(self.structure), 3))

            if use_reduced_coords:
                vcart_qnu = np.dot(vcart_qnu, self.structure.lattice.inv_matrix)

            if normalize:
                vnorm2 = f_sqrt(sum(np.linalg.norm(d) ** 2 for d in vcart_qnu))
            else:
                vnorm2 = 1.0

            # Make a bar plot with rectangles bounded by (x - width/2, x + width/2, bottom, bottom + height)
            # The align keyword controls if x is interpreted as the center or the left edge of the rectangle.
            bottom, height = 0.0, 0.0
            if atoms_index is None:
                for itype, (symbol, inds) in enumerate(symbol2indices.items()):
                    if icart is None:
                        height = f_sqrt(sum(np.linalg.norm(d) ** 2 for d in vcart_qnu[inds]) / vnorm2)
                    else:
                        height = f_sqrt(
                            sum(np.linalg.norm(d) ** 2 for ic in icart for d in vcart_qnu[inds, ic]) / vnorm2)

                    ax.bar(x, height, width, bottom, align="center",
                           color=cmap(float(itype) / max(1, ntypat - 1)),
                           label=symbol if inu == 0 else None, edgecolor='black',
                           hatch=hatches[itype % len(hatches)] if hatches else None,
                           )
                    bottom += height
            else:
                for igroup, inds in enumerate(atoms_index):
                    inds = np.array(inds)

                    if labels_groups:
                        symbol = labels_groups[igroup]
                    else:
                        symbol = "+".join("{}{}".format(self.structure[ia].specie.name, ia) for ia in inds)

                    if icart is None:
                        height = f_sqrt(sum(np.linalg.norm(d) ** 2 for d in vcart_qnu[inds]) / vnorm2)
                    else:
                        height = f_sqrt(
                            sum(np.linalg.norm(d) ** 2 for ic in icart for d in vcart_qnu[inds, ic]) / vnorm2)

                    ax.bar(x, height, width, bottom, align="center",
                           color=cmap(float(igroup) / max(1, len(atoms_index) - 1)),
                           label=symbol if inu == 0 else None, edgecolor='black',
                           hatch=hatches[igroup % len(hatches)] if hatches else None,
                           )
                    bottom += height

            xticks.append(x)
            xticklabels.append(format_w % w_qnu)
            x += (width + pad) / 2

        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_phdispl_cartdirs(self, qpoint, cart_dirs=("x", "y", "z"), units="eV",
                              is_non_analytical_direction=False, use_eigvec=False,
                              colormap="viridis", hatches="default", atoms_index=None, labels_groups=None,
                              normalize=True, use_sqrt=False, fontsize=8, branches=None, format_w="%.3f", **kwargs):
        """
        Plot three panels. Each panel shows vertical bars with the contribution of the different atomic types
        to all the phonon displacements at the given ``qpoint`` along on the Cartesian directions in ``cart_dirs``.

        Args:
            qpoint: integer, vector of reduced coordinates or |Kpoint| object.
            cart_dirs: List of strings defining the Cartesian directions. "x", "y", or "z" to select a particular
                Cartesian directions. or combinations separated by "+". Example: "x+y".
            units: Units for phonon frequencies. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            is_non_analytical_direction: If True, the ``qpoint`` is interpreted as a direction in q-space
                and the phonon (displacements/eigenvectors) for q --> 0 along this direction are used.
                Requires band structure with :class:`NonAnalyticalPh` object.
            use_eigvec: True if eigenvectors should be used instead of displacements (eigenvectors
                are orthonormal, unlike diplacements)
            colormap: Matplotlib colormap used for atom type.
            hatches: List of strings with matplotlib hatching patterns. None or empty list to disable hatching.
            fontsize: Legend and title fontsize.
            normalize: if True divides by the square norm of the total eigenvector
            use_sqrt: if True the square root of the sum of the components will be taken
                fraction of a,b,c rather than x,y,z.
            atoms_index: list of lists. Each list contains the indices of atoms in the structure that will be
                summed on a separate group. if None all the atoms will be considered and grouped by type.
            labels_groups: If atoms_index is not None will provide the labels for each of the group in atoms_index.
                Should have the same length of atoms_index or be None. If None automatic labelling will be used.
            branches: list of indices for the modes that should be represented. If None all the modes will be shown.
            format_w: string used to format the values of the frequency. Default "%.3f".

        See plot_phdispl for the meaning of the other arguments.
        """
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=len(cart_dirs), ncols=1,
                                                sharex=True, sharey=True, squeeze=False)

        for i, (cart_dir, ax) in enumerate(zip(cart_dirs, ax_list.ravel())):
            self.plot_phdispl(qpoint, cart_dir=cart_dir, ax=ax, units=units, colormap=colormap,
                              is_non_analytical_direction=is_non_analytical_direction, use_eigvec=use_eigvec,
                              fontsize=fontsize, hatches=hatches, atoms_index=atoms_index, labels_groups=labels_groups,
                              normalize=normalize, use_sqrt=use_sqrt, branches=branches, show=False, format_w=format_w)
            # Disable artists.
            if i != 0:
                set_visible(ax, False, "legend", "title")
            #if len(cart_dirs) == 3 and i != 1:
            #    set_visible(ax, False, "ylabel")
            if i != len(cart_dirs) - 1:
                set_visible(ax, False, "xlabel")

        return fig

    def get_dataframe(self, mode_range=None) -> pd.DataFrame:
        """
        Return a |pandas-DataFrame| with the following columns:

            ['qidx', 'mode', 'freq', 'qpoint']

        where:

        ==============  ==========================
        Column          Meaning
        ==============  ==========================
        qidx            q-point index.
        mode            phonon branch index.
        freq            Phonon frequency in eV.
        qpoint          |Kpoint| object
        ==============  ==========================

        Args:
            mode_range: Only modes such as `mode_range[0] <= mode_index < mode_range[1]`.
        """
        rows = []
        for iq, qpoint in enumerate(self.qpoints):
            for nu in self.branches:
                if mode_range is not None and (nu < mode_range[0] or nu >= mode_range[1]):
                    continue

                rows.append(OrderedDict([
                           ("qidx", iq),
                           ("mode", nu),
                           ("freq", self.phfreqs[iq, nu]),
                           ("qpoint", self.qpoints[iq]),
                        ]))

        df = pd.DataFrame(rows, columns=list(rows[0].keys()))
        return df

    @add_fig_kwargs
    def boxplot(self, ax=None, units="eV", mode_range=None, swarm=False, **kwargs):
        """
        Use seaborn_ to draw a box plot showing the distribution of the phonon
        frequencies with respect to the mode index.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            mode_range: Only modes such as `mode_range[0] <= mode_index < mode_range[1]` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            kwargs: Keyword arguments passed to seaborn boxplot.

        Return: |matplotlib-Figure|
        """
        df = self.get_dataframe(mode_range=mode_range)
        factor = abu.phfactor_ev2units(units)
        yname = "freq %s" % abu.phunit_tag(units)
        df[yname] = factor * df["freq"]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)
        import seaborn as sns
        hue = None
        ax = sns.boxplot(x="mode", y=yname, data=df, hue=hue, ax=ax, **kwargs)
        if swarm:
            sns.swarmplot(x="mode", y=yname, data=df, hue=hue, color=".25", ax=ax)

        return fig

    @add_plotly_fig_kwargs
    def boxplotly(self, units="eV", mode_range=None, swarm=False, fig=None, rcd=None, **kwargs):
        """
        Use plotly to draw a box plot to show the distribution of the phonon
        frequencies with respect to the mode index.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            mode_range: Only modes such as `mode_range[0] <= mode_index < mode_range[1]` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            fig: plotly figure or None if a new figure should be created.
            rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col)
                of the subplot in the grid.
            kwargs: Keyword arguments passed to plotly px.box.

        Returns: |plotly.graph_objects.Figure|
        """
        df = self.get_dataframe(mode_range=mode_range)
        factor = abu.phfactor_ev2units(units)
        yname = "freq %s" % abu.phunit_tag(units)
        df[yname] = factor * df["freq"]

        import plotly.express as px
        hue = None
        points = 'outliers' if not swarm else "all"
        px_fig = px.box(df, x="mode", y=yname, color=hue, points=points, **kwargs)

        if rcd is None: return px_fig

        # Add px_fig traces to input fig with subplot.
        rcd = PlotlyRowColDesc.from_object(rcd)
        ply_row, ply_col, iax = rcd.ply_row, rcd.ply_col, rcd.iax
        for trace in px_fig.data:
            fig.add_trace(trace, row=ply_row, col=ply_col)

        return fig

    def to_pymatgen(self, qlabels=None) -> PhononBandStructureSymmLine:
        r"""
        Creates a pymatgen :class:`PhononBandStructureSymmLine` object.

        Args:
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels e.g. ``qlabels = {(0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L"}``.
                If None labels will be determined automatically.
        """
        # pymatgen labels dict is inverted
        if qlabels is None:
            qlabels = self._auto_qlabels
            # the indices in qlabels are without the split
            labels_dict = {v: self.qpoints[k].frac_coords for k, v in qlabels.items()}
        else:
            labels_dict = {v: k for k, v in qlabels.items()}

        labelled_q_list = list(labels_dict.values())

        ph_freqs, qpts, displ = [], [], []
        for split_q, split_phf, split_phdispl in zip(self.split_qpoints, self.split_phfreqs, self.split_phdispl_cart):
            # if the qpoint has a label it needs to be repeated. If it is one of the extrema either it should
            # not be repeated (if they are the real first or last point) or they will be already repeated due
            # to the split. Also they should not be repeated in case there are two consecutive labelled points.
            # So first determine which ones have a label.
            labelled = [any(np.allclose(q, labelled_q) for labelled_q in labelled_q_list) for q in split_q]

            for i, (q, phf, d, l) in enumerate(zip(split_q, split_phf, split_phdispl, labelled)):
                ph_freqs.append(phf)
                qpts.append(q)
                d = d.reshape(self.num_branches, self.num_atoms, 3)
                displ.append(d)

                if 0 < i < len(split_q) - 1 and l and not labelled[i-1] and not labelled[i+1]:
                    ph_freqs.append(phf)
                    qpts.append(q)
                    displ.append(d)

        ph_freqs = np.transpose(ph_freqs) * abu.eV_to_THz
        qpts = np.array(qpts)
        displ = np.transpose(displ, (1, 0, 2, 3))

        return PhononBandStructureSymmLine(qpoints=qpts, frequencies=ph_freqs,
                                           lattice=self.structure.reciprocal_lattice,
                                           has_nac=self.non_anal_ph is not None, eigendisplacements=displ,
                                           labels_dict=labels_dict, structure=self.structure)

    @classmethod
    def from_pmg_bs(cls, pmg_bs: PhononBandStructureSymmLine, structure=None) -> PhononBands:
        """
        Creates an instance of the object from a :class:`PhononBandStructureSymmLine` object.

        Args:
            pmg_bs: the instance of PhononBandStructureSymmLine.
            structure: a |Structure| object. Should be present if the structure attribute is not set in pmg_bs.
        """
        structure = structure or pmg_bs.structure
        if not structure:
            raise ValueError("The structure is needed to create the abipy object.")

        structure = Structure.from_sites(structure)
        structure.spgset_abi_spacegroup(has_timerev=False)

        qpoints = []
        phfreqs = []
        phdispl_cart = []
        names = []

        prev_q = None
        for b in pmg_bs.branches:
            qname1, qname2 = b["name"].split("-")
            start_index = b["start_index"]
            if prev_q is not None and qname1 == prev_q:
                start_index += 1

            # it can happen depending on how the object was generated
            if start_index >= b["end_index"]:
                continue

            prev_q = qname2

            if start_index == b["start_index"]:
                names.append(qname1)

            names.extend([None] * (b["end_index"] - b["start_index"] - 1))
            names.append(qname2)

            for i in range(start_index, b["end_index"] + 1):
                qpoints.append(pmg_bs.qpoints[i].frac_coords)
            phfreqs.extend(pmg_bs.bands.T[start_index:b["end_index"] + 1])
            if pmg_bs.has_eigendisplacements:
                e = pmg_bs.eigendisplacements[:, start_index:b["end_index"] + 1]
                e = np.transpose(e, [0, 1, 2, 3])
                e = np.reshape(e, e.shape[:-2] + (-1,))
                phdispl_cart.extend(e)

        #print(len(names), len(phfreqs))
        qpoints_list = KpointList(reciprocal_lattice=structure.reciprocal_lattice,
                                  frac_coords=qpoints, names=names)

        phfreqs = np.array(phfreqs) / abu.eV_to_THz
        n_modes = 3 * len(structure)
        if not phdispl_cart:
            phdispl_cart = np.zeros((len(phfreqs), n_modes, n_modes))
        else:
            phdispl_cart = np.array(phdispl_cart)

        na = None
        if pmg_bs.has_nac:
            directions = []
            nac_phreqs = []
            nac_phdispl = []

            for t in pmg_bs.nac_frequencies:
                # directions in NonAnalyticalPh are given in cartesian coordinates
                cart_direction = structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(t[0])
                cart_direction = cart_direction / np.linalg.norm(cart_direction)

                directions.append(cart_direction)
                nac_phreqs.append(t[1])

            nac_phreqs = np.array(nac_phreqs) / abu.eV_to_THz

            for t in pmg_bs.nac_eigendisplacements:
                nac_phdispl.append(t[1].reshape(n_modes, n_modes))

            na = NonAnalyticalPh(structure=structure, directions=np.array(directions),
                                 phfreqs=nac_phreqs, phdispl_cart=np.array(nac_phdispl))

        phb = cls(structure=structure, qpoints=qpoints_list, phfreqs=phfreqs, phdispl_cart=phdispl_cart,
                  non_anal_ph=na)

        return phb

    def acoustic_indices(self, qpoint, threshold=0.95, raise_on_no_indices=True):
        """
        Extract the indices of the three acoustic modes for a qpoint.
        Acoustic modes could be reasonably identified for Gamma and points close to Gamma.

        Args:
            qpoint: the qpoint. Accepts integer or reduced coordinates
            threshold: fractional value allowed for the matching of the displacements to identify acoustic modes.
            raise_on_no_indices: if True a RuntimeError will be raised if the acoustic mode will not be
                correctly identified. If False [0, 1, 2] will be returned.
        """
        qindex = self.qindex(qpoint)
        phdispl = self.phdispl_cart[qindex]

        indices = []
        for mode, displ_mode in enumerate(phdispl):
            displ_mode = np.reshape(displ_mode, (-1, 3))
            a = displ_mode[0] / np.linalg.norm(displ_mode[0])
            for d in displ_mode[1:]:
                b = d / np.linalg.norm(d)
                if np.dot(a, b) < threshold:
                    break
            else:
                indices.append(mode)

        if len(indices) != 3:
            if raise_on_no_indices:
                raise RuntimeError('wrong number of indices: {}'.format(indices))
            else:
                indices = [0, 1, 2]

        return indices

    def asr_breaking(self, units='eV', threshold=0.95, raise_on_no_indices=True):
        """
        Calculates the breaking of the acoustic sum rule.
        Requires the presence of Gamma.

        Args:
            units: Units for the output. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            threshold: fractional value allowed for the matching of the displacements to identify acoustic modes.
            raise_on_no_indices: if True a RuntimeError will be raised if the acoustic mode will not be
                correctly identified

        Returns:
            A namedtuple with:
                - the three breaking of the acoustic modes
                - the maximum breaking with sign
                - the absolute value of the maximum breaking
        """
        gamma_ind = self.qpoints.index((0, 0, 0))
        ind = self.acoustic_indices(gamma_ind, threshold=threshold, raise_on_no_indices=raise_on_no_indices)
        asr_break = self.phfreqs[0, ind] * abu.phfactor_ev2units(units)

        imax = np.argmax(asr_break)

        return dict2namedtuple(breakings=asr_break, max_break=asr_break[imax], absmax_break=abs(asr_break[imax]))

    def get_frozen_phonons(self, qpoint, nmode, eta=1, scale_matrix=None, max_supercell=None):
        """
        Creates a supercell with displaced atoms for the specified q-point and mode.

        Args:
            qpoint: q vector in reduced coordinate in reciprocal space or index of the qpoint.
            nmode: number of the mode.
            eta: pre-factor multiplying the displacement. Gives the value in Angstrom of the
                largest displacement.
            scale_matrix: the scaling matrix of the supercell. If None a scaling matrix suitable for
                the qpoint will be determined.
            max_supercell: mandatory if scale_matrix is None, ignored otherwise. Defines the largest
                supercell in the search for a scaling matrix suitable for the q point.

        Returns:
            A namedtuple with a Structure with the displaced atoms, a numpy array containing the
            displacements applied to each atom and the scale matrix used to generate the supercell.
        """
        qind = self.qindex(qpoint)
        displ = self.phdispl_cart[qind, nmode].reshape((-1, 3))

        return self.structure.frozen_phonon(qpoint=self.qpoints[qind].frac_coords, displ=displ, eta=eta,
                                            frac_coords=False, scale_matrix=scale_matrix, max_supercell=max_supercell)

    def get_longitudinal_fraction(self, qpoint, idir=None):
        """
        Calculates "longitudinal" fraction of the eigendisplacements.

        Args:
            qpoint: q vector in reduced coordinate in reciprocal space or index of the qpoint.
            idir: an integer with the index of the non analytical direction if qpoint is gamma.
                If None all will be given.

        Returns:
            A numpy array with the longitudinal fractions for each mode of the specified q point.
            If qpoint is gamma and idir is None it will be a numpy array with all the non analytical
            directions.
        """
        qind = self.qindex(qpoint)
        qpoint = self.qpoints[qind]

        def get_fraction(direction, displ):
            displ = np.real(displ)
            # Normalization. Such that \sum_i dot(q, displ[i]) <= 1
            # and = 1 if q is parallel to displ[i] for each i.
            displ_norm = np.sum(np.linalg.norm(displ, axis=-1), axis=-1)
            displ = displ / displ_norm[:, None, None]
            versor = direction / np.linalg.norm(direction)
            return np.absolute(np.dot(displ, versor)).sum(axis=-1)

        if qpoint.is_gamma():
            if self.non_anal_phdispl_cart is None:
                raise RuntimeError("Cannot calculate the lo/to fraction at Gamma if the non analytical"
                                   "contributions have not been calculated.")
            phdispl = self.non_anal_phdispl_cart.reshape((len(self.non_anal_directions), self.num_branches, self.num_atoms, 3))
            if idir is None:
                fractions = []
                for non_anal_dir, phd in zip(self.non_anal_directions, phdispl):
                    fractions.append(get_fraction(non_anal_dir, phd))
                return np.array(fractions)
            else:
                return get_fraction(self.non_anal_directions[idir], phdispl[idir])
        else:
            phdispl = self.phdispl_cart[qind].reshape((self.num_branches, self.num_atoms, 3))
            return get_fraction(qpoint.cart_coords, phdispl)

    @add_fig_kwargs
    def plot_longitudinal_fraction(self, qpoint, idir=None, ax_list=None, units="eV", branches=None,
                                   format_w="%.3f", fontsize=10, **kwargs):
        """
        Plots an histogram "longitudinal" fraction of the eigendisplacements.

        Args:
            qpoint: q vector in reduced coordinate in reciprocal space or index of the qpoint.
            idir: an integer with the index of the non analytical direction if qpoint is gamma.
                If None all will be plot.
            ax_list: The axes for the plot. If ax_list is None, a new figure is created and
                the axes are automatically generated.
            units: Units for the output. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            branches: list of indices for the modes that should be represented. If None all the modes will be shown.
            format_w: string used to format the values of the frequency. Default "%.3f".
            fontsize: Labels and title fontsize.

        Returns:
            |matplotlib-Figure|

        """
        qind = self.qindex(qpoint)
        qpoint = self.qpoints[qind]
        fractions = self.get_longitudinal_fraction(qind, idir)

        factor = abu.phfactor_ev2units(units)

        if branches is None:
            branches = self.branches
        elif not isinstance(branches, (list, tuple)):
            branches = [branches]

        is_non_anal = qpoint.is_gamma()

        # if non analytical directions at gamma the
        if len(fractions.shape) == 1:
            fractions = [fractions]

        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(fractions), ncols=1,
                                                sharex=False, sharey=False, squeeze=False)

        width, pad = 4, 1
        pad = width + pad

        for i, ax in enumerate(ax_list.ravel()):
            xticks, xticklabels = [], []
            x = 0
            if idir is not None:
                i_ref = idir
            else:
                i_ref = i
            for inu, nu in enumerate(branches):
                height = fractions[i][nu]
                ax.bar(x, height, width, 0, align="center",
                       color="r", edgecolor='black')

                xticks.append(x)
                if is_non_anal:
                    w_qnu = self.non_anal_phfreqs[i_ref, nu] * factor
                else:
                    w_qnu = self.phfreqs[qind, nu] * factor
                xticklabels.append(format_w % w_qnu)

                x += (width + pad) / 2

            if is_non_anal:
                # no title for multiple axes, not enough space.
                if idir is not None:
                    ax.set_title(f"q-direction = {self.non_anal_directions[i_ref]}", fontsize=fontsize)
            else:
                ax.set_title(f"qpoint = {repr(qpoint)}", fontsize=fontsize)

            ax.set_ylabel(r"Longitudinal fraction", fontsize=fontsize)
            ax.set_ylim(0, 1)

            ax.set_xticks(xticks)
            ax.set_xticklabels((xticklabels))

            if i == len(fractions) - 1:
                ax.set_xlabel(f'Frequency {abu.phunit_tag(units)}')

        return fig

    @add_fig_kwargs
    def plot_longitudinal_fatbands(self, ax=None, units="eV", qlabels=None, branch_range=None, match_bands=False,
                                   sum_degenerate=False, factor=1, **kwargs):
        r"""
        Plot the phonon band structure with width representing the longitudinal fraction of the fatbands.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            branch_range: Tuple specifying the minimum and maximum branch index to plot (default: all branches are plotted).
            match_bands: if True the bands will be matched based on the scalar product between the eigenvectors.
            sum_degenerate: if True modes with similar frequencies will be considered as degenerated and their
                contributions will be summed (squared sum). Notice that this may end up summing contributions
                from modes that are just accidentally degenerated.
            factor: a float that will used to scale the width of the fatbands.

        Returns:
            |matplotlib-Figure|
        """
        # Select the band range.
        if branch_range is None:
            branch_range = range(self.num_branches)
        else:
            branch_range = range(branch_range[0], branch_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, units=units, qlabels=qlabels)

        if "color" not in kwargs: kwargs["color"] = "black"
        if "linewidth" not in kwargs: kwargs["linewidth"] = 1.0

        first_xx = 0

        units_factor = abu.phfactor_ev2units(units)

        for i, (q_l, pf_l) in enumerate(zip(self.split_qpoints, self.split_phfreqs)):
            if match_bands:
                ind = self.split_matched_indices[i]
                pf_l = pf_l[np.arange(len(pf_l))[:, None], ind]
            pf_l = pf_l * units_factor
            xx = list(range(first_xx, first_xx + len(pf_l)))
            for branch in branch_range:
                ax.plot(xx, pf_l[:, branch], **kwargs)
            first_xx = xx[-1]

            width = []
            for iq, (q, pf) in enumerate(zip(q_l, pf_l)):

                #print(q)
                if np.allclose(np.mod(q, 1), [0, 0, 0]):
                    if self.non_anal_ph is not None:
                        if iq == 0:
                            direction = q_l[iq+1]
                        else:
                            direction = q_l[iq-1]
                        idir = self.non_anal_ph.index_direction(direction)
                        frac = self.get_longitudinal_fraction(q, idir)
                    else:
                        frac = np.zeros(self.num_branches)
                else:
                    frac = self.get_longitudinal_fraction(q)

                # sum the contributions from degenerate modes
                if sum_degenerate:
                    pf_round = pf.round(decimals=int(6 * units_factor))
                    partitioned_pf = [np.where(pf_round == element)[0].tolist() for element in np.unique(pf_round)]
                    for group in partitioned_pf:
                        if len(group) > 1:
                            frac[group[0]] = np.linalg.norm(frac[group])
                            frac[group[1:]] = 0

                if match_bands:
                    ind = self.split_matched_indices[i]
                    frac = frac[ind[iq]]

                width.append(frac * units_factor * factor / 600)

            width = np.array(width)
            for branch in branch_range:
                ax.fill_between(xx, pf_l[:, branch] + width[:, branch], pf_l[:, branch] - width[:, branch],
                                facecolor="r", alpha=0.4, linewidth=0)

        return fig

    @add_fig_kwargs
    def plot_qpt_distance(self, qpt_list=None, ngqpt=None, shiftq=(0, 0, 0), plot_distances=False,
                          units="eV", qlabels=None, branch_range=None, colormap="viridis_r",
                          match_bands=False, log_scale=False, **kwargs):
        r"""
        Plot the phonon band structure coloring the point according to the minimum distance of
        the qpoints of the path from a list of qpoints. This can be for example defined as the
        q-points effectively calculated in DFPT.
        Optionally plot the explicit values.

        Args:
            qpt_list: list of fractional coordinates or KpointList of the qpoints from which the minimum
                distance will be calculated.
            ngqpt: the division of a regular grid of qpoints. Used to automatically fill in the qpt_list
                based on abipy.core.kpoints.kmesh_from_mpdivs.
            shiftq: the shifts of a regular grid of qpoints. Used to automatically fill in the qpt_list
                based on abipy.core.kpoints.kmesh_from_mpdivs.
            plot_distances: if True a second plot will be added with the explicit values of the distances.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            branch_range: Tuple specifying the minimum and maximum branch_i index to plot
                (default: all branches are plotted).
            colormap: matplotlib colormap to determine the colors available.
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            match_bands: if True the bands will be matched based on the scalar product between the eigenvectors.
            log_scale: if True the values will be plotted in a log scale.

        Returns: |matplotlib-Figure|
        """
        from matplotlib.collections import LineCollection

        if qpt_list is None:
            if ngqpt is None:
                raise ValueError("at least one among qpt_list and ngqpt should be provided")
            qpt_list = kmesh_from_mpdivs(ngqpt, shiftq, pbc=False, order="bz")

        if isinstance(qpt_list, KpointList):
            qpt_list = qpt_list.frac_coords

        # Select the band range.
        if branch_range is None:
            branch_range = range(self.num_branches)
        else:
            branch_range = range(branch_range[0], branch_range[1], 1)

        nrows = 2 if plot_distances else 1
        ncols = 1
        ax_list, fig, plt = get_axarray_fig_plt(ax_array=None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=True)

        # make a list in case of only one plot
        if not plot_distances:
            ax_list = [ax_list]

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax_list[-1], units=units, qlabels=qlabels)

        first_xx = 0
        factor = abu.phfactor_ev2units(units)

        linewidth = 2
        if "lw" in kwargs:
            linewidth = kwargs.pop("lw")
        elif "linewidth" in kwargs:
            linewidth = kwargs.pop("linewidth")

        rec_latt = self.structure.reciprocal_lattice

        # calculate all the value to set the color normalization
        split_min_dist = []
        for i, q_l in enumerate(self.split_qpoints):
            all_dist = rec_latt.get_all_distances(q_l, qpt_list)
            split_min_dist.append(np.min(all_dist, axis=-1))

        if log_scale:
            import matplotlib
            # find the minimum value larger than zero and set the 0 to that value
            min_value = np.min([v for l in split_min_dist for v in l if v > 0])
            for min_list in split_min_dist:
                min_list[min_list == 0] = min_value
            norm = matplotlib.colors.LogNorm(min_value, np.max(split_min_dist), clip=True)
        else:
            norm = plt.Normalize(np.min(split_min_dist), np.max(split_min_dist))

        segments = []
        total_min_dist = []

        for i, (pf, min_dist) in enumerate(zip(self.split_phfreqs, split_min_dist)):
            if match_bands:
                ind = self.split_matched_indices[i]
                pf = pf[np.arange(len(pf))[:, None], ind]
            pf = pf * factor
            xx = range(first_xx, first_xx + len(pf))

            for branch_i in branch_range:
                points = np.array([xx, pf[:, branch_i]]).T.reshape(-1, 1, 2)
                segments.append(np.concatenate([points[:-1], points[1:]], axis=1))
                total_min_dist.extend(min_dist[:-1])

            first_xx = xx[-1]

        segments = np.concatenate(segments)
        total_min_dist = np.array(total_min_dist)

        lc = LineCollection(segments, cmap=colormap, norm=norm)
        lc.set_array(total_min_dist)
        lc.set_linewidth(linewidth)

        line = ax_list[-1].add_collection(lc)

        # line collection does not autoscale the plot
        ax_list[-1].set_ylim(np.min(self.split_phfreqs), np.max(self.split_phfreqs))

        fig.colorbar(line, ax=ax_list)

        if plot_distances:
            first_xx = 0
            for i, (q_l, min_dist) in enumerate(zip(self.split_qpoints, split_min_dist)):
                xx = list(range(first_xx, first_xx + len(q_l)))
                ax_list[0].plot(xx, min_dist, linewidth=linewidth, **kwargs)

                first_xx = xx[-1]
            ax_list[0].grid(True)
        return fig


class PHBST_Reader(ETSF_Reader):
    """
    This object reads data from PHBST.nc file produced by anaddb.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PHBST_Reader
    """

    def read_qredcoords(self):
        """Array with the reduced coordinates of the q-points."""
        return self.read_value("qpoints")

    def read_qweights(self):
        """The weights of the q-points"""
        return self.read_value("qweights")

    def read_phfreqs(self):
        """|numpy-array| with the phonon frequencies in eV."""
        return self.read_value("phfreqs")

    def read_phdispl_cart(self):
        """
        Complex array with the Cartesian displacements in **Angstrom**
        shape is [num_qpoints,  mu_mode,  cart_direction].
        """
        return self.read_value("phdispl_cart", cmode="c")

    def read_amu(self):
        """The atomic mass units"""
        return self.read_value("atomic_mass_units", default=None)

    def read_epsinf_zcart(self):
        """
        Read and return electronic dielectric tensor and Born effective charges in Cartesian coordinates
        Return (None, None) if data is not available.
        """
        # nctkarr_t('emacro_cart', "dp", 'number_of_cartesian_directions, number_of_cartesian_directions')
        # nctkarr_t('becs_cart', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, number_of_atoms")]
        epsinf = self.read_value("emacro_cart", default=None)
        if epsinf is not None: epsinf = epsinf.T.copy()
        zcart = self.read_value("becs_cart", default=None)
        if zcart is not None: zcart = zcart.transpose(0, 2, 1).copy()
        return epsinf, zcart


class PhbstFile(AbinitNcFile, Has_Structure, Has_PhononBands, NotebookWriter):
    """
    Object used to access data stored in the PHBST.nc file produced by ABINIT.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PhbstFile
    """

    def __init__(self, filepath: str):
        """
        Args:
            path: path to the file
        """
        super().__init__(filepath)
        self.reader = PHBST_Reader(filepath)

        # Initialize Phonon bands and add metadata from ncfile
        self._phbands = PhononBands.from_file(filepath)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """
        String representation

        Args:
            verbose: verbosity level.
        """
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")

        app(self.phbands.to_string(title=None, with_structure=True, with_qpoints=False, verbose=verbose))

        return "\n".join(lines)

    @property
    def structure(self) -> Structure:
        """|Structure| object"""
        return self.phbands.structure

    @property
    def qpoints(self):
        """List of q-point objects."""
        return self.phbands.qpoints

    @property
    def phbands(self) -> PhononBands:
        """|PhononBands| object"""
        return self._phbands

    def close(self) -> None:
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self) -> dict:
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_phbands_params()
        return od

    def qindex(self, qpoint):
        """
        Returns the index of the qpoint in the PhbstFile.
        Accepts integer, vector with reduced coordinates or |Kpoint|.
        """
        return self.phbands.qindex(qpoint)

    def qindex_qpoint(self, qpoint, is_non_analytical_direction=False):
        """
        Returns (qindex, qpoint).
        Accepts integer, vector with reduced coordinates or |Kpoint|.
        """
        return self.phbands.qindex_qpoint(qpoint, is_non_analytical_direction=is_non_analytical_direction)

    def get_phframe(self, qpoint, with_structure=True):
        """
        Return a |pandas-DataFrame| with the phonon frequencies at the given q-point and
        information on the crystal structure (used for convergence studies).

        Args:
            qpoint: integer, vector of reduced coordinates or |Kpoint| object.
            with_structure: True to add structural parameters.
        """
        qindex, qpoint = self.qindex_qpoint(qpoint)
        phfreqs = self.phbands.phfreqs

        d = dict(
            omega=phfreqs[qindex, :],
            branch=list(range(3 * len(self.structure))),
        )

        # Add geometrical information
        if with_structure:
            d.update(self.structure.get_dict4pandas(with_spglib=True))

        # Build the pandas Frame and add the q-point as attribute.
        df = pd.DataFrame(d, columns=list(d.keys()))
        df.qpoint = qpoint

        return df

    def get_phmode(self, qpoint, branch):
        """
        Returns the :class:`PhononMode` with the given qpoint and branch nu.

        Args:
            qpoint: Either a vector with the reduced components of the q-point
                or an integer giving the sequential index (C-convention).
            branch: branch index (C-convention)

        Returns:
            :class:`PhononMode` instance.
        """
        qindex, qpoint = self.qindex_qpoint(qpoint)

        return PhononMode(qpoint=qpoint,
                          freq=self.phbands.phfreqs[qindex, branch],
                          displ_cart=self.phbands.phdispl_cart[qindex, branch, :],
                          structure=self.structure)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        return self.yield_phbands_figs(**kwargs)

    def yield_plotly_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of plotly figures with minimal input from the user.
        """
        return self.yield_phbands_plotly_figs(**kwargs)

    def write_notebook(self, nbpath=None):
        """
        Write an jupyter_ notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.phbands.plot();"),
            nbv.new_code_cell("ncfile.phbands.qpoints.plot();"),
            #nbv.new_code_cell("ncfile.phbands.get_phdos().plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


_THERMO_YLABELS = {  # [name][units] --> latex string
    "internal_energy": {"eV": "$U(T)$ (eV/cell)", "Jmol": "$U(T)$ (J/mole)"},
    "free_energy": {"eV": "$F(T) + ZPE$ (eV/cell)", "Jmol": "$F(T) + ZPE$ (J/mole)"},
    "entropy": {"eV": "$S(T)$ (eV/cell)", "Jmol": "$S(T)$ (J/mole)"},
    "cv": {"eV": "$C_V(T)$ (eV/cell)", "Jmol": "$C_V(T)$ (J/mole)"},
}

_PLOTLY_THERMO_YLABELS = {  # [name][units] --> string (no latex allowed here!)
            "internal_energy": {"eV": "U(T) (eV/cell)", "Jmol": "U(T) (J/mole)"},
            "free_energy": {"eV": "F(T) + ZPE (eV/cell)", "Jmol": "F(T) + ZPE (J/mole)"},
            "entropy": {"eV": "S(T) (eV/cell)", "Jmol": "S(T) (J/mole)"},
            "cv": {"eV": "C_V(T) (eV/cell)", "Jmol": "C_V(T) (J/mole)"},
        }


class PhononDos(Function1D):
    """
    This object stores the phonon density of states.
    An instance of ``PhononDos`` has a ``mesh`` (numpy array with the points of the mesh)
    and another numpy array, ``values``, with the DOS on the mesh.

    .. note::

        mesh is given in eV, values are in states/eV.
    """

    @classmethod
    def as_phdos(cls, obj: Any, phdos_kwargs=None) -> PhononDos:
        """
        Return an instance of |PhononDos| from a generic obj. Supports::

            - instances of cls
            - files (string) that can be open with abiopen and that provide one of the following attributes: [`phdos`, `phbands`]
            - instances of |PhononBands|.
            - objects providing a ``phbands`` attribute.

        Args:
            phdos_kwargs: optional dictionary with the options passed to ``get_phdos`` to compute the phonon DOS.
            Used when obj is not already an instance of `cls` or when we have to compute the DOS from obj.
        """
        if phdos_kwargs is None: phdos_kwargs = {}

        if isinstance(obj, cls):
            return obj

        elif is_string(obj):
            # path? (pickle or file supported by abiopen)
            if obj.endswith(".pickle"):
                with open(obj, "rb") as fh:
                    return cls.as_phdos(pickle.load(fh), phdos_kwargs)

            from abipy.abilab import abiopen
            with abiopen(obj) as abifile:
                if hasattr(abifile, "phdos"):
                    return abifile.phdos
                elif hasattr(abifile, "phbands"):
                    return abifile.phbands.get_phdos(**phdos_kwargs)
                else:
                    raise TypeError("Don't know how to create `PhononDos` from type: %s" % type(abifile))

        elif isinstance(obj, PhononBands):
            return obj.get_phdos(**phdos_kwargs)

        elif hasattr(obj, "phbands"):
            return obj.phbands.get_phdos(**phdos_kwargs)

        elif hasattr(obj, "phdos"):
            return obj.phdos

        raise TypeError("Don't know how to create PhononDos object from type: `%s`" % type(obj))

    @lazy_property
    def iw0(self) -> int:
        """
        Index of the first point in the mesh whose value is >= 0
        """
        iw0 = self.find_mesh_index(0.0)
        if iw0 == -1:
            raise ValueError("Cannot find zero in energy mesh")
        return iw0

    @lazy_property
    def idos(self):
        """Integrated DOS."""
        return self.integral()

    @lazy_property
    def zero_point_energy(self):
        """Zero point energy in eV per unit cell."""
        iw0 = self.iw0
        return Energy(0.5 * np.trapz(self.mesh[iw0:] * self.values[iw0:], x=self.mesh[iw0:]), "eV")

    def plot_dos_idos(self, ax, what="d", exchange_xy=False, units="eV", **kwargs):
        """
        Helper function to plot DOS/IDOS on the axis ``ax``.

        Args:
            ax: |matplotlib-Axes|
            what: string selecting the quantity to plot:
                "d" for DOS, "i" for IDOS. chars can be concatenated
                hence what="id" plots both IDOS and DOS. (default "d").
            exchange_xy: True to exchange axis
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            kwargs: Options passed to matplotlib plot method.

        Return:
            list of lines added to the plot.
        """
        opts = [c.lower() for c in what]
        lines = []

        for c in opts:
            f = {"d": self, "i": self.idos}[c]
            xfactor = abu.phfactor_ev2units(units)
            # Don't rescale IDOS
            yfactor = 1 / xfactor if c == "d" else 1

            ls = f.plot_ax(ax, exchange_xy=exchange_xy, xfactor=xfactor, yfactor=yfactor, **kwargs)
            lines.extend(ls)

        return lines

    def plotly_dos_idos(self, fig, what="d", trace_name=None, exchange_xy=False, units="eV", rcd=None, **kwargs):
        """
        Helper function to plot DOS/IDOS on the figure ``fig``.

        Args:
            fig: |plotly.graph_objects.Figure|
            what: string selecting the quantity to plot:
                "d" for DOS, "i" for IDOS. chars can be concatenated
                hence what="id" plots both IDOS and DOS. (default "d").
            exchange_xy: True to exchange axis
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col) of the subplot in the grid.
            kwargs: Passed to fig.add_scatter method.
        """
        opts = [c.lower() for c in what]

        for c in opts:
            f = {"d": self, "i": self.idos}[c]
            if trace_name is None:
                trace_name = {"d": 'DOS', "i": 'IDOS'}[c]
            xfactor = abu.phfactor_ev2units(units)
            # Don't rescale IDOS
            yfactor = 1 / xfactor if c == "d" else 1

            f.plotly_traces(fig, rcd=rcd, exchange_xy=exchange_xy, xfactor=xfactor, yfactor=yfactor,
                            name=trace_name, **kwargs)

    # TODO: This should be called plot_dos_idos!
    @add_fig_kwargs
    def plot(self, units="eV", **kwargs):
        """
        Plot Phonon DOS and IDOS on two distinct plots.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            kwargs: Keyword arguments passed to :mod:`matplotlib`.

        Returns: |matplotlib-Figure|
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        fig = plt.figure()
        gspec = GridSpec(2, 1, height_ratios=[1, 2], wspace=0.05)
        ax1 = plt.subplot(gspec[0])
        ax2 = plt.subplot(gspec[1])

        for ax in (ax1, ax2):
            ax.grid(True)

        ax2.set_xlabel('Energy %s' % abu.phunit_tag(units))
        ax1.set_ylabel("IDOS (states)")
        ax2.set_ylabel("DOS %s" % abu.phdos_label_from_units(units))

        self.plot_dos_idos(ax1, what="i", units=units, **kwargs)
        self.plot_dos_idos(ax2, what="d", units=units, **kwargs)

        return fig

    # TODO: This should be called plotly_dos_idos!
    @add_plotly_fig_kwargs
    def plotly(self, units="eV", **kwargs):
        """
        Plot Phonon DOS and IDOS on two distinct plots.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            kwargs: Keyword arguments passed to mod:`plotly`.

        Returns: |plotly.graph_objects.Figure|
        """
        fig, _ = get_figs_plotly(nrows=2, ncols=1, subplot_titles=[], sharex=True, sharey=False,
                                 vertical_spacing=0.05, row_heights=[1, 2])

        fig.layout['xaxis2'].title = {'text': 'Energy %s' % abu.phunit_tag(units, unicode=True)}
        fig.layout['yaxis1'].title = {'text': "IDOS (states)"}
        fig.layout['yaxis2'].title = {'text': "DOS %s" % abu.phdos_label_from_units(units, unicode=True)}

        rcd = PlotlyRowColDesc(0, 0, 2, 1)
        self.plotly_dos_idos(fig, rcd=rcd, what="i", units=units, **kwargs)
        rcd = PlotlyRowColDesc(1, 0, 2, 1)
        self.plotly_dos_idos(fig, rcd=rcd, what="d", units=units, **kwargs)

        return fig

    def get_internal_energy(self, tstart=5, tstop=300, num=50) -> Function1D:
        """
        Returns the internal energy, in eV, in the harmonic approximation for different temperatures
        Zero point energy is included.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num (int): optional Number of samples to generate. Default is 50.

        Return: |Function1D| object with U(T) + ZPE.
        """
        tmesh = np.linspace(tstart, tstop, num=num)
        w, gw = self.mesh[self.iw0:], self.values[self.iw0:]
        if w[0] < 1e-12:
            w, gw = self.mesh[self.iw0+1:], self.values[self.iw0+1:]
        coth = lambda x: 1.0 / np.tanh(x)

        vals = np.empty(len(tmesh))
        for it, temp in enumerate(tmesh):
            if temp == 0:
                vals[it] = self.zero_point_energy
            else:
                wd2kt = w / (2 * abu.kb_eVK * temp)
                vals[it] = 0.5 * np.trapz(w * coth(wd2kt) * gw, x=w)
            #print(vals[it])

        return Function1D(tmesh, vals)

    def get_entropy(self, tstart=5, tstop=300, num=50) -> Function1D:
        """
        Returns the entropy, in eV/K, in the harmonic approximation for different temperatures

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num (int): optional Number of samples to generate. Default is 50.

        Return: |Function1D| object with S(T).
        """
        tmesh = np.linspace(tstart, tstop, num=num)
        w, gw = self.mesh[self.iw0:], self.values[self.iw0:]
        if w[0] < 1e-12:
            w, gw = self.mesh[self.iw0+1:], self.values[self.iw0+1:]
        coth = lambda x: 1.0 / np.tanh(x)

        vals = np.empty(len(tmesh))
        for it, temp in enumerate(tmesh):
            if temp == 0:
                vals[it] = 0
            else:
                wd2kt = w / (2 * abu.kb_eVK * temp)
                vals[it] = np.trapz((wd2kt * coth(wd2kt) - np.log(2 * np.sinh(wd2kt))) * gw, x=w)

        return Function1D(tmesh, abu.kb_eVK * vals)

    def get_free_energy(self, tstart=5, tstop=300, num=50) -> Function1D:
        """
        Returns the free energy, in eV, in the harmonic approximation for different temperatures
        Zero point energy is included.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num (int): optional Number of samples to generate. Default is 50.

        Return: |Function1D| object with F(T) = U(T) + ZPE - T x S(T)
        """
        uz = self.get_internal_energy(tstart=tstart, tstop=tstop, num=num)
        s = self.get_entropy(tstart=tstart, tstop=tstop, num=num)

        return Function1D(uz.mesh, uz.values - s.mesh * s.values)

    def get_cv(self, tstart=5, tstop=300, num=50) -> Function1D:
        """
        Returns the constant-volume specific heat, in eV/K, in the harmonic approximation
        for different temperatures

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num (int): optional Number of samples to generate. Default is 50.

        Return: |Function1D| object with C_v(T).
        """
        tmesh = np.linspace(tstart, tstop, num=num)
        w, gw = self.mesh[self.iw0:], self.values[self.iw0:]
        if w[0] < 1e-12:
            w, gw = self.mesh[self.iw0+1:], self.values[self.iw0+1:]
        csch2 = lambda x: 1.0 / (np.sinh(x) ** 2)

        vals = np.empty(len(tmesh))
        for it, temp in enumerate(tmesh):
            if temp == 0:
                vals[it] = 0
            else:
                wd2kt = w / (2 * abu.kb_eVK * temp)
                vals[it] = np.trapz(wd2kt ** 2 * csch2(wd2kt) * gw, x=w)

        return Function1D(tmesh, abu.kb_eVK * vals)

    @add_fig_kwargs
    def plot_harmonic_thermo(self, tstart=5, tstop=300, num=50, units="eV", formula_units=None,
                             quantities="all", fontsize=8, **kwargs):
        """
        Plot thermodynamic properties from the phonon DOSes within the harmonic approximation.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 50.
            quantities: List of strings specifying the thermodynamic quantities to plot.
                Possible values: ["internal_energy", "free_energy", "entropy", "c_v"].
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            formula_units: the number of formula units per unit cell. If unspecified, the
                thermodynamic quantities will be given on a per-unit-cell basis.
            fontsize: Label and title fontsize.

        Returns: |matplotlib-Figure|
        """
        quantities = list_strings(quantities) if quantities != "all" else \
            ["internal_energy", "free_energy", "entropy", "cv"]

        # Build grid of plots.
        ncols, nrows = 1, 1
        num_plots = len(quantities)
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_mat[-1, -1].axis("off")

        for iax, (qname, ax) in enumerate(zip(quantities, ax_mat.flat)):
            irow, icol = divmod(iax, ncols)
            # Compute thermodynamic quantity associated to qname.
            f1d = getattr(self, "get_" + qname)(tstart=tstart, tstop=tstop, num=num)
            ys = f1d.values
            if formula_units is not None: ys /= formula_units
            if units == "Jmol": ys = ys * abu.e_Cb * abu.Avogadro
            ax.plot(f1d.mesh, ys)

            ax.set_title(qname, fontsize=fontsize)
            ax.grid(True)
            ax.set_xlabel("T (K)", fontsize=fontsize)
            ax.set_ylabel(_THERMO_YLABELS[qname][units], fontsize=fontsize)
            #ax.legend(loc="best", fontsize=fontsize, shadow=True)

            if irow != nrows-1:
                set_visible(ax, False, "xlabel")

        return fig

    @add_plotly_fig_kwargs
    def plotly_harmonic_thermo(self, tstart=5, tstop=300, num=50, units="eV", formula_units=None,
                               quantities="all", fontsize=12, **kwargs):
        """
        Plot thermodynamic properties from the phonon DOS within the harmonic approximation.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 50.
            quantities: List of strings specifying the thermodynamic quantities to plot.
                Possible values: ["internal_energy", "free_energy", "entropy", "c_v"].
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            formula_units: the number of formula units per unit cell. If unspecified, the
                thermodynamic quantities will be given on a per-unit-cell basis.
            fontsize: Label and title fontsize.

        Returns |plotly.graph_objects.Figure|
        """
        quantities = list_strings(quantities) if quantities != "all" else \
            ["internal_energy", "free_energy", "entropy", "cv"]

        # Build grid of plots.
        ncols, nrows = 1, 1
        num_plots = len(quantities)
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=quantities, sharex=True, sharey=False)

        for iq, qname in enumerate(quantities):
            irow, icol = divmod(iq, ncols)
            # Compute thermodynamic quantity associated to qname.
            f1d = getattr(self, "get_" + qname)(tstart=tstart, tstop=tstop, num=num)
            ys = f1d.values
            if formula_units is not None: ys /= formula_units
            if units == "Jmol": ys = ys * abu.e_Cb * abu.Avogadro
            fig.add_scatter(x=f1d.mesh, y=ys, mode="lines", name=qname, row=irow + 1, col=icol + 1)
            fig.layout.annotations[iq].font.size = fontsize
            iax = iq + 1
            fig.layout['yaxis%u' % iax].title = {'text': _PLOTLY_THERMO_YLABELS[qname][units], 'font_size': fontsize}

            if irow == nrows - 1:
                fig.layout['xaxis%u' % iax].title = {'text': 'T (K)', 'font_size': fontsize}

        return fig

    def to_pymatgen(self) -> PmgPhononDos:
        """
        Creates a pymatgen :class:`PmgPhononDos` object
        """
        factor = abu.phfactor_ev2units("thz")

        return PmgPhononDos(self.mesh * factor, self.values / factor)

    @property
    def debye_temp(self) -> float:
        """
        Debye temperature in K.
        """
        integrals = (self * self.mesh ** 2).spline_integral() / self.spline_integral()
        t_d = np.sqrt(5 / 3 * integrals) / abu.kb_eVK

        return t_d

    def get_acoustic_debye_temp(self, nsites) -> float:
        """
        Acoustic Debye temperature in K, i.e. the Debye temperature divided by nsites**(1/3).

        Args:
            nsites: the number of sites in the cell.
        """
        return self.debye_temp / nsites ** (1 / 3)


class PhdosReader(ETSF_Reader):
    """
    This object reads data from the PHDOS.nc file produced by anaddb.

    .. note::

            Frequencies are in eV, DOSes are in states/eV per unit cell.
    """
    @lazy_property
    def structure(self):
        """|Structure| object."""
        return self.read_structure()

    @lazy_property
    def wmesh(self):
        """The frequency mesh for the PH-DOS in eV."""
        return self.read_value("wmesh")

    def read_pjdos_type(self):
        """[ntypat, nomega] array with Phonon DOS projected over atom types."""
        return self.read_value("pjdos_type")

    def read_pjdos_atdir(self):
        """
        Return [natom, three, nomega] array with Phonon DOS projected over atoms and cartesian directions.
        """
        return self.read_value("pjdos")

    def read_phdos(self) -> PhononDos:
        """Return |PhononDos| object with the total phonon DOS"""
        return PhononDos(self.wmesh, self.read_value("phdos"))

    def read_pjdos_symbol_xyz_dict(self):
        """
        Return :class:`OrderedDict` mapping element symbol --> [3, nomega] array
        with the the phonon DOSes summed over atom-types and decomposed along
        the three cartesian directions.
        """
        # The name is a bit confusing: rc stands for "real-space cartesian"
        # phdos_rc_type[ntypat, 3, nomega]
        values = self.read_value("pjdos_rc_type")

        od = OrderedDict()
        for symbol in self.chemical_symbols:
            type_idx = self.typeidx_from_symbol(symbol)
            od[symbol] = values[type_idx]

        return od

    def read_pjdos_symbol_dict(self):
        """
        Ordered dictionary mapping element symbol --> |PhononDos|
        where PhononDos is the contribution to the total DOS summed over atoms
        with chemical symbol ``symbol``.
        """
        # [ntypat, nomega] array with PH-DOS projected over atom types.
        values = self.read_pjdos_type()

        od = OrderedDict()
        for symbol in self.chemical_symbols:
            type_idx = self.typeidx_from_symbol(symbol)
            od[symbol] = PhononDos(self.wmesh, values[type_idx])

        return od

    def read_msq_dos(self):
        """
        Read generalized DOS with MSQ displacement tensor in cartesian coords.

        Return: |MsqDos| object.
        """
        if "msqd_dos_atom" not in self.rootgrp.variables:
            raise RuntimeError("PHBST file does not contain `msqd_dos_atom` variable.\n" +
                               "Please use a more recent Abinit version >= 9")

        # nctkarr_t('msqd_dos_atom', "dp", 'number_of_frequencies, three, three, number_of_atoms') &
        # symmetric tensor still transpose (3,3) to be consistent.
        values = self.read_value("msqd_dos_atom").transpose([0, 2, 1, 3]).copy()

        # Read atomic masses and build dictionary element_symbol --> amu
        amu_symbol = self.read_amu_symbol()

        from abipy.dfpt.msqdos import MsqDos
        return MsqDos(self.structure, self.wmesh, values, amu_symbol)


class PhdosFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    Container object storing the different DOSes stored in the
    PHDOS.nc file produced by anaddb.
    Provides helper function to visualize/extract data.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PhdosFile
    """

    def __init__(self, filepath: str):
        # Open the file, read data and create objects.
        super().__init__(filepath)

        self.reader = r = PhdosReader(filepath)
        self.wmesh = r.wmesh

    def close(self) -> None:
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self) -> dict:
        """
        :class:`OrderedDict` with the convergence parameters
        Used to construct |pandas-DataFrames|.
        """
        return {}
        #od = OrderedDict([
        #    ("nsppol", self.nsppol),
        #])
        #return od

    def __str__(self):
        """Invoked by str"""
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """
        Human-readable string with useful information such as structure...

        Args:
            verbose: Verbosity level.
        """
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")

        return "\n".join(lines)

    @lazy_property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.reader.structure

    @lazy_property
    def phdos(self) -> PhononDos:
        """|PhononDos| object."""
        return self.reader.read_phdos()

    @lazy_property
    def pjdos_symbol(self):
        """
        Ordered dictionary mapping element symbol --> `PhononDos`
        where PhononDos is the contribution to the total DOS summed over atoms
        with chemical symbol `symbol`.
        """
        return self.reader.read_pjdos_symbol_dict()

    @lazy_property
    def msqd_dos(self):
        """
        |MsqDos| object with Mean square displacement tensor in cartesian coords.
        Allows one to calculate Debye Waller factors by integration with 1/omega and the Bose-Einstein factor.
        """
        return self.reader.read_msq_dos()

    @add_fig_kwargs
    def plot_pjdos_type(self, units="eV", stacked=True, colormap="jet", alpha=0.7, exchange_xy=False,
                        ax=None, xlims=None, ylims=None, fontsize=12, **kwargs):
        """
        Plot type-projected phonon DOS with matplotlib.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            colormap: Have a look at the colormaps
                `here <http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html>`_
                and decide which one you'd like:
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque).
            exchange_xy: True to exchange x-y axis.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: y-axis limits.
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        lw = kwargs.pop("lw", 2)
        factor = abu.phfactor_ev2units(units)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)

        ax.grid(True)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        xlabel, ylabel = 'Frequency %s' % abu.phunit_tag(units), 'PJDOS %s' % abu.phdos_label_from_units(units)
        set_ax_xylabels(ax, xlabel, ylabel, exchange_xy)

        # Type projected DOSes.
        num_plots = len(self.pjdos_symbol)
        cumulative = np.zeros(len(self.wmesh))

        for i, (symbol, pjdos) in enumerate(self.pjdos_symbol.items()):
            x, y = pjdos.mesh * factor, pjdos.values / factor
            if exchange_xy: x, y = y, x
            if num_plots != 1:
                color = cmap(float(i) / (num_plots - 1))
            else:
                color = cmap(0.0)

            if not stacked:
                ax.plot(x, y, lw=lw, label=symbol, color=color)
            else:
                if not exchange_xy:
                    ax.plot(x, cumulative + y, lw=lw, label=symbol, color=color)
                    ax.fill_between(x, cumulative, cumulative + y, facecolor=color, alpha=alpha)
                    cumulative += y
                else:
                    ax.plot(cumulative + x, y, lw=lw, label=symbol, color=color)
                    ax.fill_betweenx(y, cumulative, cumulative + x, facecolor=color, alpha=alpha)
                    cumulative += x

        # Total PHDOS
        x, y = self.phdos.mesh * factor, self.phdos.values / factor
        if exchange_xy: x, y = y, x
        ax.plot(x, y, lw=lw, label="Total PHDOS", color='black')
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_plotly_fig_kwargs
    def plotly_pjdos_type(self, units="eV", stacked=True, exchange_xy=False,
                        fig=None, xlims=None, ylims=None, fontsize=12, **kwargs):
        """
        Plot type-projected phonon DOS with plotly.

        Args:
            fig: plotly figure or None if a new figure should be created.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            exchange_xy: True to exchange x-y axis.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``.
            fontsize: legend and title fontsize.

        Returns: |plotly.graph_objects.Figure|
        """
        lw = kwargs.pop("lw", 2)
        factor = abu.phfactor_ev2units(units)

        fig, _ = get_fig_plotly(fig=fig)

        plotly_set_lims(fig, xlims, "x")
        plotly_set_lims(fig, ylims, "y")

        xlabel = 'Frequency %s' % abu.phunit_tag(units, unicode=True)
        ylabel = 'PJDOS %s' % abu.phdos_label_from_units(units, unicode=True)
        plotly_set_xylabels(fig, xlabel, ylabel, exchange_xy)

        # Type projected DOSes.
        cumulative = np.zeros(len(self.wmesh))

        for i, (symbol, pjdos) in enumerate(self.pjdos_symbol.items()):
            x, y = pjdos.mesh * factor, pjdos.values / factor
            if exchange_xy: x, y = y, x

            if not stacked:
                fig.add_scatter(x=x, y=y, mode='lines', name=symbol, line=dict(width=lw))
            else:
                if not exchange_xy:
                    fig.add_scatter(x=x, y=cumulative + y, mode='lines', name=symbol,
                                    line=dict(width=lw), fill='tonextx')
                    cumulative += y
                else:
                    fig.add_scatter(x=cumulative + x, y=y, mode='lines', name=symbol,
                                    line=dict(width=lw), fill='tonexty')
                    cumulative += x

        # Total PHDOS
        x, y = self.phdos.mesh * factor, self.phdos.values / factor
        if exchange_xy: x, y = y, x
        fig.add_scatter(x=x, y=y, mode='lines', line=dict(width=lw, color='black'), name="Total PHDOS")
        fig.layout.legend.font.size = fontsize
        fig.layout.title.font.size = fontsize

        return fig

    @add_fig_kwargs
    def plot_pjdos_cartdirs_type(self, units="eV", stacked=True, colormap="jet", alpha=0.7,
                                 xlims=None, ylims=None, ax_list=None, fontsize=8, **kwargs):
        """
        Plot type-projected phonon DOS decomposed along the three cartesian directions.
        Three rows for each cartesian direction. Each row shows the contribution of each atomic type + Total Phonon DOS.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            colormap: Have a look at the colormaps
                `here <http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html>`_
                and decide which one you'd like:
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            ylims: y-axis limits.
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and label fontsize.

        Returns: |matplotlib-Figure|.
        """
        lw = kwargs.pop("lw", 2)
        ntypat = self.structure.ntypesp
        factor = abu.phfactor_ev2units(units)

        # Three rows for each direction.
        # Each row shows the contribution of each atomic type + Total PH DOS.
        nrows, ncols = 3, 1
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=True, squeeze=True)
        ax_list = np.reshape(ax_list, (nrows, ncols)).ravel()
        cmap = plt.get_cmap(colormap)

        # symbol --> [three, number_of_frequencies] in cart dirs
        pjdos_symbol_rc = self.reader.read_pjdos_symbol_xyz_dict()

        xx = self.phdos.mesh * factor
        for idir, ax in enumerate(ax_list):
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")

            ax.set_ylabel(r'PJDOS along %s' % {0: "x", 1: "y", 2: "z"}[idir])
            if idir == 2:
                ax.set_xlabel('Frequency %s' % abu.phunit_tag(units))

            # Plot Type projected DOSes along cartesian direction idir
            cumulative = np.zeros(len(self.wmesh))
            for itype, symbol in enumerate(self.reader.chemical_symbols):
                color = cmap(float(itype) / max(1, ntypat - 1))
                yy = pjdos_symbol_rc[symbol][idir] / factor

                if not stacked:
                    ax.plot(xx, yy, label=symbol, color=color)
                else:
                    ax.plot(xx, cumulative + yy, lw=lw, label=symbol, color=color)
                    ax.fill_between(xx, cumulative, cumulative + yy, facecolor=color, alpha=alpha)
                    cumulative += yy

            # Add Total PHDOS
            ax.plot(xx, self.phdos.values / factor, lw=lw, label="Total PHDOS", color='black')
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_pjdos_cartdirs_site(self, view="inequivalent", units="eV", stacked=True, colormap="jet", alpha=0.7,
                                 xlims=None, ylims=None, ax_list=None, fontsize=8, verbose=0, **kwargs):
        """
        Plot phonon PJDOS for each atom in the unit cell. By default, only "inequivalent" atoms are shown.

        Args:
            view: "inequivalent" to show only inequivalent atoms. "all" for all sites.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            colormap: matplotlib colormap.
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and title fontsize.
            verbose: Verbosity level.

        Returns: |matplotlib-Figure|
        """
        # Define num_plots and ax2atom depending on view.
        factor = abu.phfactor_ev2units(units)
        #natom, ntypat = len(self.structure), self.structure.ntypesp
        lw = kwargs.pop("lw", 2)

        # Select atoms.
        aview = self._get_atomview(view, verbose=verbose)

        # Three rows for each cartesian direction.
        # Each row shows the contribution of each site + Total PH DOS.
        nrows, ncols = 3, 1
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=True, squeeze=True)
        ax_list = np.reshape(ax_list, (nrows, ncols)).ravel()
        cmap = plt.get_cmap(colormap)

        # [natom, three, nomega] array with PH-DOS projected over atoms and cartesian directions
        pjdos_atdir = self.reader.read_pjdos_atdir()

        xx = self.phdos.mesh * factor
        for idir, ax in enumerate(ax_list):
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")

            ax.set_ylabel(r'PJDOS along %s' % {0: "x", 1: "y", 2: "z"}[idir])
            if idir == 2:
                ax.set_xlabel('Frequency %s' % abu.phunit_tag(units))

            # Plot Type projected DOSes along cartesian direction idir
            cumulative = np.zeros(len(self.wmesh))
            for iatom in aview.iatom_list:
                site = self.structure[iatom]
                symbol = str(site)
                color = cmap(float(iatom) / max((len(aview.iatom_list) - 1), 1))
                yy = pjdos_atdir[iatom, idir] / factor

                if not stacked:
                    ax.plot(xx, yy, label=symbol, color=color)
                else:
                    ax.plot(xx, cumulative + yy, lw=lw, label=symbol, color=color)
                    ax.fill_between(xx, cumulative, cumulative + yy, facecolor=color, alpha=alpha)
                    cumulative += yy

            # Add Total PHDOS
            ax.plot(xx, self.phdos.values / factor, lw=lw, label="Total PHDOS", color='black')
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        units = kwargs.get("units", "mev")
        yield self.phdos.plot(units=units, show=False)
        yield self.plot_pjdos_type(units=units, show=False)
        # Old formats do not have MSQDOS arrays.
        try:
            msqd_dos = self.msqd_dos
        except Exception:
            msqd_dos = None
        if msqd_dos is not None:
            yield msqd_dos.plot(units=units, show=False)
            yield msqd_dos.plot_tensor(show=False)

    def yield_plotly_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of plotly figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        units = kwargs.get("units", "mev")
        yield self.phdos.plotly(units=units, show=False)
        yield self.plotly_pjdos_type(units=units, show=False)
        # Old formats do not have MSQDOS arrays.
        #try:
        #    msqd_dos = self.msqd_dos
        #except Exception:
        #    msqd_dos = None
        #if msqd_dos is not None:
        #    yield msqd_dos.plot(units=units, show=False)
        #    yield msqd_dos.plot_tensor(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If ``nbpath`` is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.phdos.plot();"),
            nbv.new_code_cell("ncfile.plot_pjdos_type();"),
            nbv.new_code_cell("ncfile.plot_pjdos_cartdirs_type(units='meV', stacked=True);"),
            nbv.new_code_cell("ncfile.plot_pjdos_cartdirs_site(view='inequivalent', units='meV', stacked=True);"),
            # TODO
            #msqd_dos = self.msqd_dos
            #msqd_dos.plot(units=self.units, show=False)
            #msqd_dos.plot_tensor(show=False)
        ])

        return self._write_nb_nbpath(nb, nbpath)

    def to_pymatgen(self) -> PmgCompletePhononDos:
        """
        Creates a pymatgen :class:`PmgCompletePhononDos` object.
        """
        total_dos = self.phdos.to_pymatgen()

        # [natom, three, nomega] array with PH-DOS projected over atoms and cartesian directions"""
        pjdos_atdir = self.reader.read_pjdos_atdir()

        factor = abu.phfactor_ev2units("thz")
        summed_pjdos = np.sum(pjdos_atdir, axis=1) / factor

        pdoss = {site: pdos for site, pdos in zip(self.structure, summed_pjdos)}

        return PmgCompletePhononDos(self.structure, total_dos, pdoss)


# FIXME: Remove. Use PhononBandsPlotter API.
@add_fig_kwargs
def phbands_gridplot(phb_objects, titles=None, phdos_objects=None, phdos_kwargs=None,
                     units="eV", width_ratios=(2, 1), fontsize=8, **kwargs):
    """
    Plot multiple phonon bandstructures and optionally DOSes on a grid.

    Args:
        phb_objects: List of objects from which the phonon band structures are extracted.
            Each item in phb_objects is either a string with the path of the netcdf file,
            or one of the abipy object with an ``phbands`` attribute or a |PhononBands| object.
        phdos_objects: List of objects from which the phonon DOSes are extracted.
            Accept filepaths or |PhononDos| objects. If phdos_objects is not None,
            each subplot in the grid contains a band structure with DOS else a simple bandstructure plot.
        titles: List of strings with the titles to be added to the subplots.
        phdos_kwargs: optional dictionary with the options passed to ``get_phdos`` to compute the phonon DOS.
            Used only if ``phdos_objects`` is not None.
        units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
            Case-insensitive.
        width_ratios: Ratio between the width of the phonon band plots and the DOS plots.
            Used if `phdos_objects` is not None
        fontsize: legend and title fontsize.

    Returns: |matplotlib-Figure|
    """
    # Build list of PhononBands objects.
    phbands_list = [PhononBands.as_phbands(obj) for obj in phb_objects]

    # Build list of PhononDos objects.
    phdos_list = []
    if phdos_objects is not None:
        if phdos_kwargs is None: phdos_kwargs = {}
        phdos_list = [PhononDos.as_phdos(obj, phdos_kwargs) for obj in phdos_objects]
        if len(phdos_list) != len(phbands_list):
            raise ValueError("The number of objects for DOS must equal be to the number of bands")

    import matplotlib.pyplot as plt
    nrows, ncols = 1, 1
    numeb = len(phbands_list)
    if numeb > 1:
        ncols = 2
        nrows = numeb // ncols + numeb % ncols

    if not phdos_list:
        # Plot grid with phonon bands only.
        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()
        # don't show the last ax if numeb is odd.
        if numeb % ncols != 0: ax_list[-1].axis("off")

        for i, (phbands, ax) in enumerate(zip(phbands_list, ax_list)):
            phbands.plot(ax=ax, units=units, show=False)
            if titles is not None: ax.set_title(titles[i], fontsize=fontsize)
            if i % ncols != 0:
                ax.set_ylabel("")

    else:
        # Plot grid with phonon bands + DOS
        # see http://matplotlib.org/users/gridspec.html
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure()
        gspec = GridSpec(nrows, ncols)

        for i, (phbands, phdos) in enumerate(zip(phbands_list, phdos_list)):
            subgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=gspec[i], width_ratios=width_ratios, wspace=0.05)
            # Get axes and align bands and DOS.
            ax1 = plt.subplot(subgrid[0])
            ax2 = plt.subplot(subgrid[1], sharey=ax1)
            phbands.plot_with_phdos(phdos, ax_list=(ax1, ax2), units=units, show=False)

            if titles is not None: ax1.set_title(titles[i], fontsize=fontsize)
            if i % ncols != 0:
                for ax in (ax1, ax2):
                    ax.set_ylabel("")

    return fig


def dataframe_from_phbands(phbands_objects, index=None, with_spglib=True) -> pd.DataFrame:
    """
    Build pandas dataframe with the most important results available in a list of band structures.

    Args:
        phbands_objects: List of objects that can be converted to phonon bands objects..
            Support netcdf filenames or |PhononBands| objects
            See ``PhononBands.as_phbands`` for the complete list.
        index: Index of the dataframe.
        with_spglib: If True, spglib is invoked to get the spacegroup symbol and number.

    Return: |pandas-DataFrame|
    """
    phbands_list = [PhononBands.as_phbands(obj) for obj in phbands_objects]
    # Use OrderedDict to have columns ordered nicely.
    odict_list = [(phbands.get_dict4pandas(with_spglib=with_spglib)) for phbands in phbands_list]

    return pd.DataFrame(odict_list, index=index,
                        columns=list(odict_list[0].keys()) if odict_list else None)


class PhononBandsPlotter(NotebookWriter):
    """
    Class for plotting phonon band structure and DOSes.
    Supports plots on the same graph or separated plots.

    Usage example:

    .. code-block:: python

        plotter = PhononBandsPlotter()
        plotter.add_phbands("foo bands", "foo_PHBST.nc")
        plotter.add_phbands("bar bands", "bar_PHBST.nc")
        plotter.gridplot()
    """
    # Used in iter_lineopt to generate matplotlib linestyles.
    _LINE_COLORS = ["blue", "red", "green", "magenta", "yellow", "black"]
    _LINE_STYLES = ["-", ":", "--", "-.",]
    _LINE_STYLES_PLOTLY = ['solid', "dot", 'dash', 'dashdot',]
    _LINE_WIDTHS = [2, ]

    def __init__(self, key_phbands=None, key_phdos=None, phdos_kwargs=None):
        """
        Args:
            key_phbands: List of (label, phbands) tuples.
                phbands is any object that can be converted into |PhononBands| e.g. ncfile, path.
            key_phdos: List of (label, phdos) tuples.
                phdos is any object that can be converted into |PhononDos|.
        """
        if key_phbands is None: key_phbands = []
        key_phbands = [(k, PhononBands.as_phbands(v)) for k, v in key_phbands]
        self._bands_dict = OrderedDict(key_phbands)

        if key_phdos is None: key_phdos = []
        key_phdos = [(k, PhononDos.as_phdos(v, phdos_kwargs)) for k, v in key_phdos]
        self._phdoses_dict = OrderedDict(key_phdos)
        if key_phdos:
            if not key_phbands:
                raise ValueError("key_phbands must be specifed when key_dos is not None")
            if len(key_phbands) != len(key_phdos):
                raise ValueError("key_phbands and key_phdos must have the same number of elements.")

    def __len__(self):
        return len(self._bands_dict)

    def __repr__(self):
        """Invoked by repr"""
        return self.to_string(func=repr)

    def __str__(self):
        """Invoked by str"""
        return self.to_string(func=str)

    def append_plotter(self, other):
        """Append phbands and phdos from other plotter to self."""
        for label, phbands in other._bands_dict.items():
            phdos = self._phdoses_dict.get(label, None)
            self.add_phbands(label, phbands, phdos=phdos)

    def add_plotter(self, other: PhononBandsPlotter) -> PhononBandsPlotter:
        """Merge two plotters, return new plotter."""
        if not isinstance(other, self.__class__):
            raise TypeError("Don't know to to add %s to %s" % (other.__class__, self.__class__))

        key_phbands = list(self._bands_dict.items()) + list(other._bands_dict.items())
        key_phdos = list(self._phdoses_dict.items()) + list(other._phdoses_dict.items())

        return self.__class__(key_phbands=key_phbands, key_phdos=key_phdos)

    def to_string(self, func=str, verbose: int = 0) -> str:
        """String representation."""
        lines = []
        app = lines.append
        for i, (label, phbands) in enumerate(self.phbands_dict.items()):
            app("[%d] %s --> %s" % (i, label, func(phbands)))

        if self.phdoses_dict and verbose:
            for i, (label, phdos) in enumerate(self.phdoses_dict.items()):
                app("[%d] %s --> %s" % (i, label, func(phdos)))

        return "\n".join(lines)

    def has_same_formula(self) -> bool:
        """
        True of plotter contains structures with same chemical formula.
        """
        structures = [phbands.structure for phbands in self.phbands_dict.values()]
        if structures and any(s.formula != structures[0].formula for s in structures): return False
        return True

    def get_phbands_frame(self, with_spglib=True) -> pd.DataFrame:
        """
        Build a |pandas-DataFrame| with the most important results available in the band structures.
        """
        return dataframe_from_phbands(list(self.phbands_dict.values()),
                                      index=list(self.phbands_dict.keys()), with_spglib=with_spglib)

    @property
    def phbands_dict(self) -> dict:
        """Dictionary with the mapping label --> phbands."""
        return self._bands_dict

    # TODO: Just an alias. To be removed in 0.4
    bands_dict = phbands_dict

    @property
    def phdoses_dict(self) -> dict:
        """Dictionary with the mapping label --> phdos."""
        return self._phdoses_dict

    @property
    def phbands_list(self) -> List[PhononBands]:
        """"List of |PhononBands| objects."""
        return list(self._bands_dict.values())

    @property
    def phdoses_list(self) -> List[PhononDos]:
        """"List of |PhononDos|."""
        return list(self._phdoses_dict.values())

    def iter_lineopt(self):
        """Generates matplotlib linestyles."""
        for o in itertools.product(self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def iter_lineopt_plotly(self):
        """Generates plotly linestyles."""
        for o in itertools.product(self._LINE_WIDTHS,  self._LINE_STYLES_PLOTLY, self._LINE_COLORS):
            yield {"line_width": o[0], "line_dash": o[1], "line_color": o[2]}

    def add_phbands(self, label, bands, phdos=None, dos=None, phdos_kwargs=None) -> None:
        """
        Adds a band structure for plotting.

        Args:
            label: label for the bands. Must be unique.
            bands: |PhononBands| object.
            phdos: |PhononDos| object.
            phdos_kwargs: optional dictionary with the options passed to ``get_phdos`` to compute the phonon DOS.
              Used only if ``phdos`` is not None.
        """
        if dos is not None:
            warnings.warn("dos has been renamed phdos. The argument will removed in abipy 0.4")
            if phdos is not None:
                raise ValueError("phdos and dos are mutually exclusive")
            phdos = dos

        if label in self._bands_dict:
            raise ValueError("label %s is already in %s" % (label, list(self._bands_dict.keys())))

        self._bands_dict[label] = PhononBands.as_phbands(bands)

        if phdos is not None:
            self.phdoses_dict[label] = PhononDos.as_phdos(phdos, phdos_kwargs)

    @add_fig_kwargs
    def combiplot(self, qlabels=None, units='eV', ylims=None, width_ratios=(2, 1), fontsize=8,
                  linestyle_dict=None, **kwargs):
        r"""
        Plot the band structure and the DOS on the same figure with matplotlib.
        Use ``gridplot`` to plot band structures on different figures.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the k-points.
                The values are the labels e.g. ``klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            width_ratios: Ratio between the width of the phonon bands plots and the DOS plots.
                Used if plotter has DOSes.
            fontsize: fontsize for legend.
            linestyle_dict: Dictionary mapping labels to matplotlib linestyle options.

        Returns: |matplotlib-Figure|
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        # Build grid of plots.
        fig = plt.figure()
        if self.phdoses_dict:
            gspec = GridSpec(1, 2, width_ratios=width_ratios, wspace=0.05)
            ax1 = plt.subplot(gspec[0])
            # Align bands and DOS.
            ax2 = plt.subplot(gspec[1], sharey=ax1)
            ax_list = [ax1, ax2]

        else:
            ax1 = fig.add_subplot(111)
            ax_list = [ax1]

        for ax in ax_list:
            ax.grid(True)

        if ylims is not None:
            for ax in ax_list:
                set_axlims(ax, ylims, "y")

        # Plot phonon bands.
        lines, legends = [], []
        my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        nqpt_list = [phbands.nqpt for phbands in self._bands_dict.values()]
        if any(nq != nqpt_list[0] for nq in nqpt_list):
            cprint("WARNING combiblot: Bands have different number of k-points:\n%s" % str(nqpt_list), "yellow")

        for (label, phbands), lineopt in zip(self._bands_dict.items(), self.iter_lineopt()):
            i += 1
            if linestyle_dict is not None and label in linestyle_dict:
                my_kwargs.update(linestyle_dict[label])
            else:
                my_kwargs.update(lineopt)
            opts_label[label] = my_kwargs.copy()

            l = phbands.plot_ax(ax1, branch=None, units=units, **my_kwargs)
            lines.append(l[0])

            # Use relative paths if label is a file.
            if os.path.isfile(label):
                legends.append("%s" % os.path.relpath(label))
            else:
                legends.append("%s" % label)

            # Set ticks and labels, legends.
            if i == 0:
                phbands.decorate_ax(ax1, qlabels=qlabels, units=units)

        ax1.legend(lines, legends, loc='best', fontsize=fontsize, shadow=True)

        # Add DOSes
        if self.phdoses_dict:
            ax = ax_list[1]
            for label, dos in self.phdoses_dict.items():
                dos.plot_dos_idos(ax, exchange_xy=True, units=units, **opts_label[label])

        return fig

    @add_plotly_fig_kwargs
    def combiplotly(self, qlabels=None, units='eV', ylims=None, width_ratios=(2, 1), fontsize=12,
                  linestyle_dict=None, **kwargs):
        r"""
        Plot the band structure and the DOS on the same figure with plotly.
        Use ``gridplotply`` to plot band structures on different figures.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the k-points.
                The values are the labels e.g. ``klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``.
            width_ratios: Ratio between the width of the phonon bands plots and the DOS plots.
                Used if plotter has DOSes.
            fontsize: fontsize for titles and legend.
            linestyle_dict: Dictionary mapping labels to linestyle options passed to |plotly.graph_objects.scatter|.

        Returns: |plotly.graph_objects.Figure|
        """
        if self.phdoses_dict:
            nrows, ncols = (1, 2)
            fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=[], sharex=False, sharey=True,
                                     horizontal_spacing=0.02, column_widths=width_ratios)
        else:
            nrows, ncols = (1, 1)
            fig, _ = get_fig_plotly()

        plotly_set_lims(fig, ylims, 'y')

        # Plot phonon bands.
        my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        nqpt_list = [phbands.nqpt for phbands in self._bands_dict.values()]
        if any(nq != nqpt_list[0] for nq in nqpt_list):
            cprint("WARNING combiblot: Bands have different number of k-points:\n%s" % str(nqpt_list), "yellow")

        for (label, phbands), lineopt in zip(self._bands_dict.items(), self.iter_lineopt_plotly()):
            i += 1
            if linestyle_dict is not None and label in linestyle_dict:
                my_kwargs.update(linestyle_dict[label])
            else:
                my_kwargs.update(lineopt)
            opts_label[label] = my_kwargs.copy()

            # Use relative paths if label is a file.
            if os.path.isfile(label): label = os.path.relpath(label)

            rcd = PlotlyRowColDesc(0, 0, nrows, ncols)
            phbands.plotly_traces(fig, branch=None, rcd=rcd, units=units, name=label, showlegend=True, **my_kwargs)

            # Set ticks and labels, legends.
            if i == 0:
                phbands.decorate_plotly(fig, qlabels=qlabels, units=units, iax=rcd.iax)

        fig.layout.legend.font.size = fontsize
        fig.layout.title.font.size = fontsize

        # Add DOSes
        if self.phdoses_dict:
            rcd = PlotlyRowColDesc(0, 1, nrows, ncols)
            for label, dos in self.phdoses_dict.items():
                dos.plotly_dos_idos(fig, rcd=rcd, exchange_xy=True, units=units, trace_name=label, legendgroup=label,
                                    showlegend=False, **opts_label[label])

        return fig

    def plot(self, *args, **kwargs):
        """An alias for combiplot."""
        return self.combiplot(*args, **kwargs)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """This function *generates* a predefined list of matplotlib figures with minimal input from the user."""
        yield self.gridplot(show=False)
        yield self.boxplot(show=False)
        if self.has_same_formula():
            yield self.combiplot(show=False)
            yield self.combiboxplot(show=False)

    def yield_plotly_figs(self, **kwargs):  # pragma: no cover
        """This function *generates* a predefined list of matplotlib figures with minimal input from the user."""
        yield self.gridplotly(show=False)
        #yield self.boxplotly(show=False)
        if self.has_same_formula():
            yield self.combiplotly(show=False)
            #yield self.combiboxplotly(show=False)

    @add_fig_kwargs
    def gridplot(self, with_dos=True, units="eV", fontsize=8, **kwargs):
        """
        Plot multiple phonon bandstructures and optionally DOSes on a grid with matplotlib.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            with_dos: True to plot phonon DOS (if available).
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        titles = list(self._bands_dict.keys())
        phb_objects = list(self._bands_dict.values())
        phdos_objects = None
        if self.phdoses_dict and with_dos:
            phdos_objects = list(self.phdoses_dict.values())

        return phbands_gridplot(phb_objects, titles=titles, phdos_objects=phdos_objects,
                                units=units, fontsize=fontsize, show=False)

    @add_plotly_fig_kwargs
    def gridplotly(self, with_dos=True, units="eV", fontsize=12, **kwargs):
        """
        Plot multiple phonon bandstructures and optionally DOSes on a grid with plotly.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            with_dos: True to plot phonon DOS (if available).
            fontsize: legend and title fontsize.

        Returns: |plotly.graph_objects.Figure|
        """
        titles = list(self._bands_dict.keys())
        phb_objects = list(self._bands_dict.values())
        phdos_objects = None
        plot_with_phdos = False
        if self.phdoses_dict and with_dos:
            phdos_objects = list(self.phdoses_dict.values())
            plot_with_phdos = True

        # Build grid of plots.
        ncols, nrows = 1, 1
        num_plots = len(phb_objects)

        if plot_with_phdos:
            # Special treatment required for phbands with DOS.
            num_plots *= 2
            titles = []
            for k in self._bands_dict.keys():
                titles.extend([k, ""])

        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        sharex = not plot_with_phdos
        fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=titles, sharex=sharex, sharey=False)

        if plot_with_phdos:
            #print("Warning: plot_with_phdos is still under development!!!!!!!!!!")
            jj = 0
            for i, phbands in enumerate(phb_objects):
                phdos = phdos_objects[i]
                row, col = divmod(jj, ncols)
                jj += 2
                rcd_phbands = PlotlyRowColDesc(row, col, nrows, ncols)
                rcd_phdos = PlotlyRowColDesc(row, col + 1, nrows, ncols)
                phbands.plotly_with_phdos(phdos, fig=fig, rcd_phbands=rcd_phbands, rcd_phdos=rcd_phdos,
                                          units=units, fontsize=fontsize,
                                          width_ratios=(2, 1), show=False)
        else:
            for i, phbands in enumerate(phb_objects):
                row, col = divmod(i, ncols)
                rcd = PlotlyRowColDesc(row, col, nrows, ncols)
                phbands.plotly(fig=fig, rcd=rcd, units=units, fontsize=fontsize, show=False)

        return fig

    @add_fig_kwargs
    def gridplot_with_hue(self, hue, with_dos=False, units="eV", width_ratios=(2, 1),
                          ylims=None, fontsize=8, **kwargs):
        """
        Plot multiple phonon bandstructures and optionally DOSes on a grid.
        Group results by ``hue``.

        Example:

            plotter.gridplot_with_hue("tsmear")

        Args:
            hue: Variable that define subsets of the phonon bands, which will be drawn on separate plots.
                Accepts callable or string
                If string, it's assumed that the phbands has an attribute with the same name and getattr is invoked.
                Dot notation is also supported e.g. hue="structure.formula" --> abifile.structure.formula
                If callable, the output of hue(phbands) is used.
            with_dos: True to plot phonon DOS (if available).
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            width_ratios: Ratio between the width of the fatbands plots and the DOS plots.
                Used if plotter has PH DOSes is not None
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                or scalar e.g. `left`. If left (right) is None, default values are used
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        # Extract all quantities available in the plotter to prepare grouping.
        all_labels = list(self._bands_dict.keys())
        all_phb_objects = list(self._bands_dict.values())
        all_phdos_objects = None
        if self.phdoses_dict and with_dos:
            all_phdos_objects = list(self.phdoses_dict.values())

        # Need index to handle all_phdos_objects if DOSes are wanted.
        if callable(hue):
            items = [(hue(phb), phb, i, label) for i, (phb, label) in enumerate(zip(all_phb_objects, all_labels))]
        else:
            # Assume string. Either phbands.hue or phbands.params[hue].
            if duck.hasattrd(all_phb_objects[0], hue):
                items = [(duck.getattrd(phb, hue), phb, i, label)
                        for i, (phb, label) in enumerate(zip(all_phb_objects, all_labels))]
            else:
                items = [(phb.params[hue], phb, i, label)
                        for i, (phb, label) in enumerate(zip(all_phb_objects, all_labels))]

        # Group items by hue value.
        hvalues, groups = sort_and_groupby(items, key=lambda t: t[0], ret_lists=True)
        nrows, ncols = len(groups), 1

        if not all_phdos_objects:
            # Plot grid with phonon bands only.
            ax_phbands, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                       sharex=True, sharey=True, squeeze=False)
            ax_phbands = ax_phbands.ravel()

            # Loop over groups
            for ax, hvalue, grp in zip(ax_phbands, hvalues, groups):
                # Unzip items
                # See https://stackoverflow.com/questions/19339/transpose-unzip-function-inverse-of-zip
                _, phb_list, indices, labels = tuple(map(list, zip(*grp)))
                assert len(phb_list) == len(indices) and len(phb_list) == len(labels)
                ax.grid(True)
                sh = str(hue) if not callable(hue) else str(hue.__doc__)
                ax.set_title("%s = %s" % (sh, hvalue), fontsize=fontsize)

                nqpt_list = [phbands.nqpt for phbands in phb_list]
                if any(nq != nqpt_list[0] for nq in nqpt_list):
                    cprint("WARNING: Bands have different number of k-points:\n%s" % str(nqpt_list), "yellow")

                # Plot all bands in grups on the same axis.
                for i, (phbands, lineopts) in enumerate(zip(phb_list, self.iter_lineopt())):
                    # Plot all branches with lineopts and set the label of the last line produced.
                    phbands.plot_ax(ax, branch=None, units=units, **lineopts)
                    ax.lines[-1].set_label(labels[i])

                    if i == 0:
                        # Set ticks and labels
                        phbands.decorate_ax(ax, qlabels=None, units=units)

                # Set legends.
                ax.legend(loc='best', fontsize=fontsize, shadow=True)
                set_axlims(ax, ylims, "y")

        else:
            # Plot grid with phonon bands + DOS (grouped by hue)
            # see http://matplotlib.org/users/gridspec.html
            import matplotlib.pyplot as plt
            from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
            fig = plt.figure()
            gspec = GridSpec(nrows, ncols)

            # Loop over groups
            for i, (hvalue, grp) in enumerate(zip(hvalues, groups)):
                # Unzip items
                _, phb_list, indices, labels = tuple(map(list, zip(*grp)))
                assert len(phb_list) == len(indices) and len(phb_list) == len(labels)

                subgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=gspec[i], width_ratios=width_ratios, wspace=0.05)
                # Get axes and align bands and DOS.
                ax1 = plt.subplot(subgrid[0])
                ax2 = plt.subplot(subgrid[1], sharey=ax1)

                sh = str(hue) if not callable(hue) else str(hue.__doc__)
                ax1.set_title("%s = %s" % (sh, hvalue), fontsize=fontsize)

                # Plot all bands in grups on the same axis.
                nqpt_list = [phbands.nqpt for phbands in phb_list]
                if any(nq != nqpt_list[0] for nq in nqpt_list):
                    cprint("WARNING: Bands have different number of k-points:\n%s" % str(nqpt_list), "yellow")

                phdos_list = [all_phdos_objects[j] for j in indices]
                for j, (phbands, phdos, lineopts) in enumerate(zip(phb_list, phdos_list, self.iter_lineopt())):
                    # Plot all branches with DOS and lineopts and set the label of the last line produced
                    phbands.plot_with_phdos(phdos, ax_list=(ax1, ax2), units=units, show=False, **lineopts)
                    ax1.lines[-1].set_label(labels[j])

                # Set legends on ax1
                ax1.legend(loc='best', fontsize=fontsize, shadow=True)

                for ax in (ax1, ax2):
                    set_axlims(ax, ylims, "y")

        return fig

    @add_fig_kwargs
    def boxplot(self, mode_range=None, units="eV", swarm=False, **kwargs):
        """
        Use seaborn_ to draw a box plot to show distribution of eigenvalues with respect to the band index.
        Band structures are drawn on different subplots.

        Args:
            mode_range: Only bands such as ``mode_range[0] <= nu_index < mode_range[1]`` are included in the plot.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            swarm: True to show the datapoints on top of the boxes
            kwargs: Keywork arguments passed to seaborn_ boxplot.
        """
        # Build grid of plots.
        num_plots, ncols, nrows = len(self.phbands_dict), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # don't show the last ax if numeb is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        for (label, phbands), ax in zip(self.phbands_dict.items(), ax_list):
            phbands.boxplot(ax=ax, units=units, mode_range=mode_range, show=False)
            ax.set_title(label)

        return fig

    @add_fig_kwargs
    def combiboxplot(self, mode_range=None, units="eV", swarm=False, ax=None, **kwargs):
        """
        Use seaborn_ to draw a box plot comparing the distribution of the frequencies.
        Phonon Band structures are drawn on the same plot.

        Args:
            mode_range: Only bands such as ``mode_range[0] <= nu_index < mode_range[1]`` are included in the plot.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            swarm: True to show the datapoints on top of the boxes
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: Keyword arguments passed to seaborn_ boxplot.
        """
        df_list = []
        for label, phbands in self.phbands_dict.items():
            # Get the dataframe, select bands and add column with label
            frame = phbands.get_dataframe()
            if mode_range is not None:
                frame = frame[(frame["mode"] >= mode_range[0]) & (frame["mode"] < mode_range[1])]
            frame["label"] = label
            df_list.append(frame)

        # Merge df_list ignoring index (not meaningful here)
        data = pd.concat(df_list, ignore_index=True)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        # Create column with frequencies in `units`.
        factor = abu.phfactor_ev2units(units)
        yname = "freq %s" % abu.phunit_tag(units)
        data[yname] = factor * data["freq"]

        import seaborn as sns
        sns.boxplot(x="mode", y=yname, data=data, hue="label", ax=ax, **kwargs)
        if swarm:
            sns.swarmplot(x="mode", y=yname, data=data, hue="label", color=".25", ax=ax)

        return fig

    @add_fig_kwargs
    def plot_phdispl(self, qpoint, **kwargs):
        """
        Plot vertical bars with the contribution of the different atomic types to the phonon displacements
        at a given q-point. One panel for all |PhononBands| stored in the plotter.

        Args:
            qpoint: integer, vector of reduced coordinates or |Kpoint| object.
            kwargs: keyword arguments passed to phbands.plot_phdispl

        Returns: |matplotlib-Figure|
        """
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=len(self.phbands_dict), ncols=1,
                                                sharex=False, sharey=False, squeeze=False)

        for i, (ax, (label, phbands)) in enumerate(zip(ax_list.ravel(), self.phbands_dict.items())):
            phbands.plot_phdispl(qpoint, cart_dir=None, ax=ax, show=False, **kwargs)
            # Disable artists.
            if i != 0:
                #set_visible(ax, False, "title")
                ax.set_title(label, fontsize=kwargs.get("fontsize", 8))
            if i != len(self.phbands_dict) - 1:
                set_visible(ax, False, "xlabel")

        return fig

    def animate(self, interval=500, savefile=None, units="eV", width_ratios=(2, 1), show=True):
        """
        Use matplotlib to animate a list of band structure plots (with or without DOS).

        Args:
            interval: draws a new frame every interval milliseconds.
            savefile: Use e.g. 'myanimation.mp4' to save the animation in mp4 format.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            width_ratios: Ratio between the band structure plot and the dos plot.
                Used when there are DOS stored in the plotter.
            show: True if the animation should be shown immediately.

        Returns: Animation object.

        .. Seealso::

            http://matplotlib.org/api/animation_api.html
            http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

        .. Note::

            It would be nice to animate the title of the plot, unfortunately
            this feature is not available in the present version of matplotlib.
            See: http://stackoverflow.com/questions/17558096/animated-title-in-matplotlib
        """
        phbands_list, phdos_list = self.phbands_list, self.phdoses_list
        if phdos_list and len(phdos_list) != len(phbands_list):
            raise ValueError("The number of objects for DOS must be equal to the number of bands")
        #titles = list(self.phbands_dict.keys())

        import matplotlib.pyplot as plt
        fig = plt.figure()
        plotax_kwargs = {"color": "black", "linewidth": 2.0}

        artists = []
        if not phdos_list:
            # Animation with band structures
            ax = fig.add_subplot(1, 1, 1)
            phbands_list[0].decorate_ax(ax, units=units)
            for i, phbands in enumerate(phbands_list):
                lines = phbands.plot_ax(ax=ax, branch=None, units=units, **plotax_kwargs)
                #if titles is not None: lines += [ax.set_title(titles[i])]
                artists.append(lines)
        else:
            # Animation with band structures + DOS.
            from matplotlib.gridspec import GridSpec
            gspec = GridSpec(1, 2, width_ratios=width_ratios, wspace=0.05)
            ax1 = plt.subplot(gspec[0])
            ax2 = plt.subplot(gspec[1], sharey=ax1)
            phbands_list[0].decorate_ax(ax1)
            ax2.grid(True)
            ax2.yaxis.set_ticks_position("right")
            ax2.yaxis.set_label_position("right")

            for i, (phbands, phdos) in enumerate(zip(phbands_list, phdos_list)):
                phbands_lines = phbands.plot_ax(ax=ax1, branch=None, units=units, **plotax_kwargs)
                phdos_lines = phdos.plot_dos_idos(ax=ax2, units=units, exchange_xy=True, **plotax_kwargs)
                lines = phbands_lines + phdos_lines
                #if titles is not None: lines += [ax.set_title(titles[i])]
                artists.append(lines)

        import matplotlib.animation as animation
        anim = animation.ArtistAnimation(fig, artists, interval=interval,
                                         blit=False, # True is faster but then the movie starts with an empty frame!
                                         #repeat_delay=1000
                                         )

        if savefile is not None: anim.save(savefile)
        if show: plt.show()

        return anim

    def ipw_select_plot(self): # pragma: no cover
        """
        Return an ipython widget with controllers to select the plot.
        """
        def plot_callback(plot_type, units):
            r = getattr(self, plot_type)(units=units, show=True)
            if plot_type == "animate": return r

        import ipywidgets as ipw
        return ipw.interact_manual(
                plot_callback,
                plot_type=["combiplot", "gridplot", "boxplot", "combiboxplot", "animate"],
                units=["eV", "cm-1", "Ha"],
            )

    def _repr_html_(self):
        """Integration with jupyter_ notebooks."""
        return self.ipw_select_plot()

    def get_panel(self):
        """Return tabs with widgets to interact with the |PhononBandsPlotter| file."""
        from abipy.panels.phonons import PhononBandsPlotterPanel
        return PhononBandsPlotterPanel(self).get_panel()

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        # Use pickle files for data persistence.
        tmpfile = self.pickle_dump()

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("plotter = abilab.PhononBandsPlotter.pickle_load('%s')" % tmpfile),
            nbv.new_code_cell("print(plotter)"),
            nbv.new_code_cell("frame = plotter.get_phbands_frame()\ndisplay(frame)"),
            nbv.new_code_cell("plotter.ipw_select_plot()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class PhononDosPlotter(NotebookWriter):
    """
    Class for plotting multiple phonon DOSes.

    Usage example:

    .. code-block:: python

        plotter = PhononDosPlotter()
        plotter.add_phdos("foo dos", "foo.nc")
        plotter.add_phdos("bar dos", "bar.nc")
        plotter.gridplot()
    """
    def __init__(self, key_phdos=None, phdos_kwargs=None):
        self._phdoses_dict = OrderedDict()
        if key_phdos is None: key_phdos = []
        for label, phdos in key_phdos:
            self.add_phdos(label, phdos, phdos_kwargs=phdos_kwargs)

    @property
    def phdos_list(self) -> List[PhononDos]:
        """List of phonon DOSes"""
        return list(self._phdoses_dict.values())

    def add_phdos(self, label, phdos, phdos_kwargs=None) -> None:
        """
        Adds a DOS for plotting.

        Args:
            label: label for the phonon DOS. Must be unique.
            phdos: |PhononDos| object.
            phdos_kwargs: optional dictionary with the options passed to `get_phdos` to compute the phonon DOS.
                Used when phdos is not already an instance of `cls` or when we have to compute the DOS from obj.
        """
        if label in self._phdoses_dict:
            raise ValueError("label %s is already in %s" % (label, list(self._phdoses_dict.keys())))

        self._phdoses_dict[label] = PhononDos.as_phdos(phdos, phdos_kwargs)

    #def has_same_formula(self):
    #    """
    #    True of plotter contains structures with the same chemical formula.
    #    """
    #    structures = [phdos.structure for phdos in self._phdoses_dict.values()]
    #    if structures and any(s.formula != structures[0].formula for s in structures): return False
    #    return True

    @add_fig_kwargs
    def combiplot(self, ax=None, units="eV", xlims=None, ylims=None, fontsize=8, **kwargs):
        """
        Plot DOSes on the same figure. Use ``gridplot`` to plot DOSes on different figures.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            ylims: y-axis limits.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        ax.set_xlabel('Energy %s' % abu.phunit_tag(units))
        ax.set_ylabel('DOS %s' % abu.phdos_label_from_units(units))

        lines, legends = [], []
        for label, dos in self._phdoses_dict.items():
            l = dos.plot_dos_idos(ax, units=units, **kwargs)[0]
            lines.append(l)
            legends.append("DOS: %s" % label)

        # Set legends.
        ax.legend(lines, legends, loc='best', fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def combiplotly(self, fig=None, units="eV", xlims=None, ylims=None, fontsize=8, **kwargs):
        """
        Plot DOSes on the same plotly figure. Use ``gridplotly`` to plot DOSes on different figures.

        Args:
            fig: plotly figure or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``.
            fontsize: Legend and title fontsize.

        Returns: |plotly.graph_objects.Figure|
        """
        fig, _ = get_fig_plotly(fig=fig)
        plotly_set_lims(fig, xlims, "x")
        plotly_set_lims(fig, ylims, "y")

        fig.layout['xaxis1'].title = {'text': 'Energy %s' % abu.phunit_tag(units, unicode=True)}
        fig.layout['yaxis1'].title = {"text": 'DOS %s' % abu.phdos_label_from_units(units, unicode=True)}

        for label, dos in self._phdoses_dict.items():
            dos.plotly_dos_idos(fig, units=units, trace_name="DOS: %s" % label)

        # Set legends.
        #ax.legend(lines, legends, loc='best', fontsize=fontsize, shadow=True)

        return fig

    def plot(self, **kwargs):
        """An alias for combiplot."""
        return self.combiplot(**kwargs)

    @add_fig_kwargs
    def gridplot(self, units="eV", xlims=None, ylims=None, fontsize=8, **kwargs):
        """
        Plot multiple DOSes on a grid with matplotlib.

        Args:
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Axis_label and subtitle fontsize.

        Returns: |matplotlib-Figure|
        """
        titles = list(self._phdoses_dict.keys())
        phdos_list = list(self._phdoses_dict.values())

        nrows, ncols = 1, 1
        numeb = len(phdos_list)
        if numeb > 1:
            ncols = 2
            nrows = numeb // ncols + numeb % ncols

        # Build Grid
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()

        # don't show the last ax if numeb is odd.
        if numeb % ncols != 0: ax_list[-1].axis("off")

        for i, (label, phdos) in enumerate(self._phdoses_dict.items()):
            ax = ax_list[i]
            phdos.plot_dos_idos(ax, units=units)

            ax.set_xlabel('Energy %s' % abu.phunit_tag(units), fontsize=fontsize)
            ax.set_ylabel("DOS %s" % abu.phdos_label_from_units(units), fontsize=fontsize)
            ax.set_title(label, fontsize=fontsize)
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")
            if i % ncols != 0:
                ax.set_ylabel("")

        return fig

    @add_plotly_fig_kwargs
    def gridplotly(self, units="eV", xlims=None, ylims=None, fontsize=12, **kwargs):
        """
        Plot multiple DOSes on a grid with plotly.

        Args:
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
            fontsize: Axis_label and subtitle fontsize.

        Returns: |plotly.graph_objects.Figure|
        """
        titles = list(self._phdoses_dict.keys())
        phdos_list = list(self._phdoses_dict.values())

        nrows, ncols = 1, 1
        numeb = len(phdos_list)
        if numeb > 1:
            ncols = 2
            nrows = numeb // ncols + numeb % ncols

        # Build Grid Fig
        fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=titles, sharex=True, sharey=True)

        x_unit = abu.phunit_tag(units, unicode=True)
        y_unit = abu.phdos_label_from_units(units, unicode=True)
        for i, (label, phdos) in enumerate(self._phdoses_dict.items()):
            row, col = divmod(i, ncols)
            rcd = PlotlyRowColDesc(row, col, nrows, ncols)
            phdos.plotly_dos_idos(fig, rcd=rcd, units=units, trace_name=label, showlegend=False)
            fig.layout['xaxis'+str(rcd.iax)].title = {'text': 'Energy %s' % x_unit, "font": {"size" : fontsize}}
            if col%ncols==0:
                fig.layout['yaxis'+str(rcd.iax)].title = {"text": 'DOS %s' % y_unit, "font": {"size" : fontsize}}
            fig.layout.annotations[rcd.iax-1].font.size = fontsize
            plotly_set_lims(fig, xlims, "x")
            plotly_set_lims(fig, ylims, "y")

        return fig

    @add_fig_kwargs
    def plot_harmonic_thermo(self, tstart=5, tstop=300, num=50, units="eV", formula_units=1,
                             quantities="all", fontsize=8, **kwargs):
        """
        Plot thermodynamic properties from the phonon DOS within the harmonic approximation
        for all the files in the plotter with matplotlib.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 50.
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            formula_units: the number of formula units per unit cell. If unspecified, the
                thermodynamic quantities will be given on a per-unit-cell basis.
            quantities: List of strings specifying the thermodynamic quantities to plot.
                Possible values in ["internal_energy", "free_energy", "entropy", "c_v"].
            fontsize: Legend, axis_label and subtitles fontsize.

        Returns: |matplotlib-Figure|
        """
        quantities = list_strings(quantities) if quantities != "all" else \
            ["internal_energy", "free_energy", "entropy", "cv"]

        # Build grid of plots.
        ncols, nrows = 1, 1
        num_plots = len(quantities)
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_mat[-1, -1].axis("off")

        for iax, (qname, ax) in enumerate(zip(quantities, ax_mat.flat)):
            for i, (label, phdos) in enumerate(self._phdoses_dict.items()):
                # Compute thermodynamic quantity associated to qname.
                f1d = getattr(phdos, "get_" + qname)(tstart=tstart, tstop=tstop, num=num)
                ys = f1d.values
                if formula_units != 1: ys /= formula_units
                if units == "Jmol": ys = ys * abu.e_Cb * abu.Avogadro
                ax.plot(f1d.mesh, ys, label=label)

            ax.set_title(qname, fontsize=fontsize)
            ax.grid(True)
            ax.set_ylabel(_THERMO_YLABELS[qname][units], fontsize=fontsize)
            ax.set_xlabel("T (K)", fontsize=fontsize)
            if iax == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_plotly_fig_kwargs
    def plotly_harmonic_thermo(self, tstart=5, tstop=300, num=50, units="eV", formula_units=1,
                             quantities="all", fontsize=12, **kwargs):
        """
        Plot thermodynamic properties from the phonon DOS within the harmonic approximation
        for all the files in the plotter with plotly.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 50.
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            formula_units: the number of formula units per unit cell. If unspecified, the
                thermodynamic quantities will be given on a per-unit-cell basis.
            quantities: List of strings specifying the thermodynamic quantities to plot.
                Possible values in ["internal_energy", "free_energy", "entropy", "c_v"].
            fontsize: Legend, axis_label and subtitle fontsize.

        Returns: |plotly.graph_objects.Figure|
        """
        quantities = list_strings(quantities) if quantities != "all" else \
            ["internal_energy", "free_energy", "entropy", "cv"]

        # Build grid of plots.
        ncols, nrows = 1, 1
        num_plots = len(quantities)
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=quantities, sharex=True, sharey=False)

        import plotly.colors as pcolors
        l2color = pcolors.DEFAULT_PLOTLY_COLORS

        for iq, qname in enumerate(quantities):
            irow, icol = divmod(iq, ncols)
            for i, (label, phdos) in enumerate(self._phdoses_dict.items()):
                opt = {"color": l2color[i]}
                # Compute thermodynamic quantity associated to qname.
                f1d = getattr(phdos, "get_" + qname)(tstart=tstart, tstop=tstop, num=num)
                ys = f1d.values
                if formula_units != 1: ys /= formula_units
                if units == "Jmol": ys = ys * abu.e_Cb * abu.Avogadro
                if iq == 0:
                    fig.add_scatter(x=f1d.mesh, y=ys, mode="lines", name=label, legendgroup=label, showlegend=True,
                                    line=opt, row=irow + 1, col=icol + 1)
                else:
                    fig.add_scatter(x=f1d.mesh, y=ys, mode="lines", name=label, legendgroup=label, showlegend=False,
                                    line=opt, row=irow + 1, col=icol + 1)

            fig.layout.annotations[iq].font.size = fontsize
            fig.layout.legend.font.size = fontsize

            iax = iq + 1
            fig.layout['yaxis%u' % iax].title = {'text': _PLOTLY_THERMO_YLABELS[qname][units], 'font_size': fontsize}

            if irow == nrows - 1:
                fig.layout['xaxis%u' % iax].title = {'text': 'T (K)', 'font_size': fontsize}

        return fig

    def ipw_select_plot(self): # pragma: no cover
        """
        Return an ipython widget with controllers to select the plot.
        """
        def plot_callback(plot_type, units):
            getattr(self, plot_type)(units=units, show=True)

        import ipywidgets as ipw
        return ipw.interact_manual(
                plot_callback,
                plot_type=["combiplot", "gridplot"],
                units=["eV", "meV", "cm-1", "Thz", "Ha"],
            )

    def ipw_harmonic_thermo(self): # pragma: no cover
        """
        Return an ipython widget with controllers to plot thermodynamic properties
        from the phonon DOS within the harmonic approximation.
        """
        def plot_callback(tstart, tstop, num, units, formula_units):
            self.plot_harmonic_thermo(tstart=tstart, tstop=tstop, num=num,
                                      units=units, formula_units=formula_units, show=True)

        import ipywidgets as ipw
        return ipw.interact_manual(
                plot_callback,
                tstart=5, tstop=300, num=50, units=["eV", "Jmol"], formula_units=1)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.gridplot(show=False)
        yield self.plot_harmonic_thermo(show=False)
        #if self.has_same_formula():
        yield self.combiplot(show=False)

    def yield_plotly_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        #yield self.gridplotply(show=False)
        #yield self.plotly_harmonic_thermo(show=False)
        #if self.has_same_formula():
        yield self.combiplotly(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write an jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        # Use pickle files for data persistence.
        tmpfile = self.pickle_dump()

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("plotter = abilab.ElectronDosPlotter.pickle_load('%s')" % tmpfile),
            nbv.new_code_cell("print(plotter)"),
            nbv.new_code_cell("plotter.ipw_select_plot()"),
            nbv.new_code_cell("plotter.ipw_harmonic_thermo()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class RobotWithPhbands(object):
    """
    Mixin class for robots associated to files with |PhononBands|.
    """
    def combiplot_phbands(self, **kwargs):
        """Wraps combiplot method of |PhononBandsPlotter|. kwargs passed to combiplot."""
        return self.get_phbands_plotter().combiplot(**kwargs)

    def gridplot_phbands(self, **kwargs):
        """Wraps gridplot method of |PhononBandsPlotter|. kwargs passed to gridplot."""
        return self.get_phbands_plotter().gridplot(**kwargs)

    def boxplot_phbands(self, **kwargs):
        """Wraps boxplot method of |PhononBandsPlotter|. kwargs passed to boxplot."""
        return self.get_phbands_plotter().boxplot(**kwargs)

    def combiboxplot_phbands(self, **kwargs):
        """Wraps combiboxplot method of |PhononBandsPlotter|. kwargs passed to combiboxplot."""
        return self.get_phbands_plotter().combiboxplot(**kwargs)

    #def combiplot_phdos(self, **kwargs):
    #    """Wraps combiplot method of |ElectronDosPlotter|. kwargs passed to combiplot."""
    #    return self.get_phdos_plotter().combiplot(**kwargs)
    #
    #def gridplot_phdos(self, **kwargs):
    #    """Wraps gridplot method of |ElectronDosPlotter|. kwargs passed to gridplot."""
    #    return self.get_phdos_plotter().gridplot(**kwargs)

    def get_phbands_plotter(self, filter_abifile=None, cls=None) -> PhononBandsPlotter:
        """
        Build and return an instance of |PhononBandsPlotter| or a subclass is cls is not None.

        Args:
            filter_abifile: Function that receives an ``abifile`` object and returns
                True if the file should be added to the plotter.
            cls: subclass of |PhononBandsPlotter|
        """
        plotter = PhononBandsPlotter() if cls is None else cls()

        for label, abifile in self.items():
            if filter_abifile is not None and not filter_abifile(abifile): continue
            plotter.add_phbands(label, abifile.phbands)

        return plotter

    def get_phbands_dataframe(self, with_spglib=True) -> pd.DataFrame:
        """
        Build a |pandas-dataframe| with the most important results available in the band structures.
        """
        return dataframe_from_phbands([nc.phbands for nc in self.abifiles],
                                      index=self.labels, with_spglib=with_spglib)

    @add_fig_kwargs
    def plot_phdispl(self, qpoint, **kwargs):
        """
        Plot vertical bars with the contribution of the different atomic types to the phonon displacements
        at a given q-point. One panel for all phbands stored in the plotter.

        Args:
            qpoint: integer, vector of reduced coordinates or |Kpoint| object.
            kwargs: keyword arguments passed to phbands.plot_phdispl

        Returns: |matplotlib-Figure|
        """
        return self.get_phbands_plotter().plot_phdispl(qpoint, show=False, **kwargs)

    def get_phbands_code_cells(self, title=None):
        """Return list of notebook cells."""
        # Try not pollute namespace with lots of variables.
        nbformat, nbv = self.get_nbformat_nbv()
        title = "## Code to compare multiple PhononBands objects" if title is None else str(title)
        return [
            nbv.new_markdown_cell(title),
            nbv.new_code_cell("robot.get_phbands_plotter().ipw_select_plot();"),
            nbv.new_code_cell("#robot.plot_phdispl(qpoint=(0, 0, 0));"),
        ]


# TODO: PhdosRobot
class PhbstRobot(Robot, RobotWithPhbands):
    """
    This robot analyzes the results contained in multiple PHBST.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PhbstRobot
    """
    EXT = "PHBST"

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        plotter = self.get_phbands_plotter()
        for fig in plotter.yield_figs(): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If ``nbpath`` is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.PhbstRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
        ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_phbands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)
