# coding: utf-8
from __future__ import print_function, division, absolute_import # unicode_literals,

import sys
import functools
import numpy as np
import itertools
import pickle
import os
import six
import json
import warnings
import abipy.core.abinit_units as abu

from collections import OrderedDict
from monty.string import is_string, list_strings, marquee
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
from monty.termcolor import cprint
from pymatgen.core.units import eV_to_Ha, Energy
from pymatgen.core.periodic_table import Element
from abipy.core.func1d import Function1D
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_PhononBands, NotebookWriter
from abipy.core.kpoints import Kpoint, Kpath
from abipy.abio.robots import Robot
from abipy.iotools import ETSF_Reader
from abipy.tools import gaussian, duck
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims, get_axarray_fig_plt, set_visible, set_ax_xylabels
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos as PmgCompletePhononDos, PhononDos as PmgPhononDos


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
class PhononMode(object):
    """
    A phonon mode has a q-point, a frequency, a cartesian displacement and a |Structure|.
    """

    __slots__ = [
        "qpoint",
        "freq",
        "displ_cart", # Cartesian displacement.
        "structure"
    ]

    def __init__(self, qpoint, freq, displ_cart, structure):
        """
        Args:
            qpoint: qpoint in reduced coordinates.
            freq: Phonon frequency in eV.
            displ: Displacement (Cartesian coordinates in Angstrom)
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

    def to_string(self, with_displ=True, verbose=0):
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


class PhononBands(object):
    """
    Container object storing the phonon band structure.

    .. note::

        Frequencies are in eV. Cartesian displacements are in Angstrom.
    """
    @classmethod
    def from_file(cls, filepath):
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
            # TODO: Reading NonAnalyticalPh here is not safe because
            # it may happen that the netcdf file does not contain all the directions
            # required by AbiPy. For the time being we read NonAnalyticalPh only
            # if we know that calculation has been driven by AbiPy --> all directions are available.
            #if "non_analytical_directions" in r.rootgrp.variables:
            #    print("Found nonanal")
            #    non_anal_ph = NonAnalyticalPh.from_file(filepath)

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
    def as_phbands(cls, obj):
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
    def phfactor_ev2units(units):
        """
        Return conversion factor eV --> units (case-insensitive)
        """
        return abu.phfactor_ev2units(units)

    def read_non_anal_from_file(self, filepath):
        """
        Reads the non analytical directions, frequencies and displacements from the anaddb.nc file
        specified and adds them to the object.
        """
        self.non_anal_ph = NonAnalyticalPh.from_file(filepath)

    def __init__(self, structure, qpoints, phfreqs, phdispl_cart, non_anal_ph=None, amu=None,
                 epsinf=None, zcart=None, linewidths=None):
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

        # Dictionary with metadata e.g. nkpt, tsmear ...
        self.params = OrderedDict()

    # TODO: Replace num_qpoints with nqpt, deprecate num_qpoints
    @property
    def nqpt(self):
        """An alias for num_qpoints."""
        return  self.num_qpoints

    def __repr__(self):
        """String representation (short version)"""
        return "<%s, nk=%d, %s, id=%s>" % (self.__class__.__name__, self.num_qpoints, self.structure.formula, id(self))

    def __str__(self):
        return self.to_string()

    def to_string(self, title=None, with_structure=True, with_qpoints=False, verbose=0):
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

    def __add__(self, other):
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
    def minfreq(self):
        """Minimum phonon frequency."""
        return self.get_minfreq_mode()

    @property
    def maxfreq(self):
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
    def has_linewidths(self):
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

    def to_xmgrace(self, filepath, units="meV"):
        """
        Write xmgrace_ file with phonon band structure energies and labels for high-symmetry q-points.

        Args:
            filepath: String with filename or stream.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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
        w('@world xmax %d' % (self.num_qpoints -1))
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

    def qindex(self, qpoint):
        """
	Returns the index of the qpoint. Accepts integer or reduced coordinates.
	"""
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

    def get_dict4pandas(self, with_spglib=True):
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

    def get_phdos(self, method="gaussian", step=1.e-4, width=4.e-4):
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
            iqpt: index of qpoint in self
            filename: name of the XYZ file that will be created
            pre_factor: Multiplication factor of the displacements
            do_real: True if we want only real part of the displacement, False means imaginary part
            scale_matrix: Scaling matrix of the supercell
            max_supercell: Maximum size of the supercell with respect to primitive cell
        """
        if scale_matrix is None:
            if max_supercell is None:
                raise ValueError("If scale_matrix is not provided, please provide max_supercell !")

            scale_matrix = self.structure.get_smallest_supercell(self.qpoints[iqpt].frac_coords, max_supercell=max_supercell)

        natoms = int(np.round(len(self.structure)*np.linalg.det(scale_matrix)))
        with open(filename, "wt") as xyz_file:
            for imode in np.arange(self.num_branches):
                xyz_file.write(str(natoms) + "\n")
                xyz_file.write("Mode " + str(imode) + " : " + str(self.phfreqs[iqpt, imode]) + "\n")
                self.structure.write_vib_file(
                    xyz_file, self.qpoints[iqpt].frac_coords, pre_factor * np.reshape(self.phdispl_cart[iqpt, imode,:],(-1,3)),
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

            displ_list = np.zeros((self.num_branches, self.num_atoms, 3), dtype=np.complex)
            for i in range(self.num_atoms):
                displ_list[:,i,:] = self.phdispl_cart[iqpt,:,3*i:3*(i+1)] * \
                    np.exp(-2*np.pi*1j*np.dot(structure[i].frac_coords, self.qpoints[iqpt].frac_coords))

            displ_list = np.dot(np.dot(displ_list, structure.lattice.inv_matrix), ascii_basis) * pre_factor

            for imode in np.arange(self.num_branches):
                lines.append("#metaData: qpt=[{:.6f};{:.6f};{:.6f};{:.6f} \\".format(
                    q[0], q[1], q[2], self.phfreqs[iqpt, imode]))

                for displ in displ_list[imode]:
                    line = "#; "+ "; ".join("{:.6f}".format(i) for i in displ.real) + "; " \
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

        # http://henriquemiranda.github.io/phononwebsite/index.html
        data = {}
        data["name"] = name or self.structure.composition.reduced_formula
        data["natoms"] = self.num_atoms
        data["lattice"] = self.structure.lattice.matrix.tolist()
        data["atom_types"] = [e.name for e in self.structure.species]
        data["atom_numbers"] = self.structure.atomic_numbers
        data["formula"] = self.structure.formula.replace(" ", "")
        data["repetitions"] = repetitions or (3, 3, 3)
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
                data["highsym_qpts"] = list(six.moves.zip(*self._make_ticks_and_labels(None)))
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
            v /= np.linalg.norm(v[0,0,0])
            v = np.stack([v.real, v.imag], axis=-1)

            vectors.extend(v.tolist())

        data["qpoints"] = qpoints
        data["distances"] = distances
        data["eigenvalues"] = eigenvalues
        data["vectors"] = vectors
        #print("name", data["name"], "\nhighsym_qpts:", data["highsym_qpts"])

        with open(filename, 'wt') as json_file:
            json.dump(data, json_file, indent=indent)

    def decorate_ax(self, ax, units='eV', **kwargs):
        """
        Add q-labels, title and unit name to axis ax.
        Use units = "" to add k-labels without adding unit name.

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

    @add_fig_kwargs
    def plot(self, ax=None, units="eV", qlabels=None, branch_range=None, match_bands=False, temp=None,
             fontsize=12, **kwargs):
        r"""
        Plot the phonon band structure.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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
        if branch_range is None:
            branch_range = range(self.num_branches)
        else:
            branch_range = range(branch_range[0], branch_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, units=units, qlabels=qlabels)

        if "color" not in kwargs: kwargs["color"] = "black"
        if "linewidth" not in kwargs: kwargs["linewidth"] = 2.0

        # Plot the phonon branches.
        self.plot_ax(ax, branch_range, units=units, match_bands=match_bands, **kwargs)

        if temp is not None:
            # Scatter plot with Bose-Einstein occupation factor for T = temp
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

    @add_fig_kwargs
    def plot_colored_matched(self, ax=None, units="eV", qlabels=None, branch_range=None,
                             colormap="rainbow", max_colors=None, **kwargs):
        r"""
        Plot the phonon band structure with different color for each line.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            branch_range: Tuple specifying the minimum and maximum branch_i index to plot
                (default: all branches are plotted).
            colormap: matplotlib colormap to determine the colors available. The colors will be chosen not in a
                sequential order to avoid difficulties in distinguishing the lines.
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            max_colors: maximum number of colors to be used. If max_colors < num_braches the colors will be reapeated.
                It useful to better distinguish close bands when the number of branch is large.

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
    def plot_lt_character(self, units="eV", qlabels=None, ax=None, xlims=None, ylims=None,
                          colormap="jet", fontsize=12, **kwargs):
        r"""
        Plot the phonon band structure with colored lines. The color of the lines indicates
        the degree to which the mode is longitudinal:
        Red corresponds to longitudinal modes and black to purely transverse modes.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. ``qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: y-axis limits.
            colormap: Matplotlib colormap.
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        if self.zcart is None:
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
                # q x Z[atom] x disp[q, nu, atom]
                for nu in range(self.num_branches):
                    v = sum(np.dot(qcart, np.dot(self.zcart[iatom], dis[nu, iatom])) for iatom in range(self.num_atoms))
                    scatt_x.append(xx[iq])
                    scatt_y.append(ws[nu])
                    scatt_s.append(v * inv_qepsq)

            p_freqs = p_freqs * factor
            ax.plot(xx, p_freqs, **kwargs)
            first_xx = xx[-1]

        scatt_y = np.array(scatt_y) * factor
        scatt_s = np.abs(np.array(scatt_s))
        scatt_s /= scatt_s.max()
        scatt_s *= 50
        print("scatt_s", scatt_s, "min", scatt_s.min(), "max", scatt_s.max())

        ax.scatter(scatt_x, scatt_y, s=scatt_s,
            #c=None, marker=None, cmap=None, norm=None, vmin=None, vmax=None, alpha=None,
            #linewidths=None, verts=None, edgecolors=None, *, data=None
        )
        self.decorate_ax(ax, units=units, qlabels=None)
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
            #     ind_block = np.zeros((len(displ), self.num_branches), dtype=np.int)
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
                ind_block = np.zeros((len(displ), self.num_branches), dtype=np.int)
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
        #TODO should be modified in order to handle the "split" list of qpoints
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
        Plot phonon fatbands and, optionally, atom-projected phonon DOSes.
        The width of the band is given by ||v_{type}||
        where v is the (complex) phonon displacement (eigenvector) in cartesian coordinates and
        v_{type} selects only the terms associated to the atomic type.

        Args:
            use_eigvec: True if the width of the phonon branch should be computed from the eigenvectors.
                False to use phonon displacements. Note that the PHDOS is always decomposed in
                terms of (orthonormal) eigenvectors.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            colormap: Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            phdos_file: Used to activate fatbands + PJDOS plot.
                Accept string with path of PHDOS.nc file or :class:`PhdosFile` object.
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

    @add_fig_kwargs
    def plot_with_phdos(self, phdos, units="eV", qlabels=None, ax_list=None, width_ratios=(2, 1), **kwargs):
        r"""
        Plot the phonon band structure with the phonon DOS.

        Args:
            phdos: An instance of |PhononDos| or a netcdf file providing a PhononDos object.
            units: Units for plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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
        ax1.yaxis.set_view_interval(emin, emax)

        # Plot Phonon DOS
        phdos.plot_dos_idos(ax2, what="d", units=units, exchange_xy=True, **kwargs)

        ax2.grid(True)
        ax2.yaxis.set_ticks_position("right")
        #ax2.yaxis.set_label_position("right")

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
            branches: list of indices for the modes that shoul be represented. If None all the modes will be shown.
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
        for nu in branches:
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
                           label=symbol if nu == 0 else None, edgecolor='black',
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
                           label=symbol if nu == 0 else None, edgecolor='black',
                           hatch=hatches[igroup % len(hatches)] if hatches else None,
                           )
                    bottom += height

            xticks.append(x)
            xticklabels.append(format_w % w_qnu)
            x += (width + pad) / 2

        ax.set_xticks(xticks)
        ax.set_xticklabels((xticklabels))
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
            branches: list of indices for the modes that shoul be represented. If None all the modes will be shown.
            format_w: string used to format the values of the frequency. Default "%.3f".

        See plot_phdispl for the meaning of the other arguments.
        """
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=len(cart_dirs), ncols=1,
                                                sharex=True, sharey=True, squeeze=False)

        for i, (cart_dir, ax) in enumerate(zip(cart_dirs, ax_list.ravel())):
            self.plot_phdispl(qpoint, cart_dir=cart_dir, ax=ax, units=units, colormap=colormap,
                              is_non_analytical_direction=is_non_analytical_direction, use_eigvec=use_eigvec,
                              fontsize=fontsize, hatches=hatches, atoms_index=atoms_index, labels_groups=labels_groups,
                              normalize=normalize, use_sqrt=use_sqrt, branches=branches, show=False)
            # Disable artists.
            if i != 0:
                set_visible(ax, False, "legend", "title")
            #if len(cart_dirs) == 3 and i != 1:
            #    set_visible(ax, False, "ylabel")
            if i != len(cart_dirs) - 1:
                set_visible(ax, False, "xlabel")

        return fig

    def get_dataframe(self):
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
        """
        import pandas as pd
        rows = []
        for iq, qpoint in enumerate(self.qpoints):
            for nu in self.branches:
                rows.append(OrderedDict([
                           ("qidx", iq),
                           ("mode", nu),
                           ("freq", self.phfreqs[iq, nu]),
                           ("qpoint", self.qpoints[iq]),
                        ]))

        return pd.DataFrame(rows, columns=list(rows[0].keys()))

    @add_fig_kwargs
    def boxplot(self, ax=None, units="eV", mode_range=None, swarm=False, **kwargs):
        """
        Use seaborn_ to draw a box plot to show distributions of eigenvalues with respect to the mode index.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            mode_range: Only modes such as `mode_range[0] <= mode_index < mode_range[1]` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            kwargs: Keyword arguments passed to seaborn boxplot.
        """
        # Get the dataframe and select bands
        frame = self.get_dataframe()
        if mode_range is not None:
            frame = frame[(frame["mode"] >= mode_range[0]) & (frame["mode"] < mode_range[1])]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        factor = abu.phfactor_ev2units(units)
        yname = "freq %s" % abu.phunit_tag(units)
        frame[yname] = factor * frame["freq"]

        import seaborn as sns
        hue = None
        ax = sns.boxplot(x="mode", y=yname, data=frame, hue=hue, ax=ax, **kwargs)
        if swarm:
            sns.swarmplot(x="mode", y=yname, data=frame, hue=hue, color=".25", ax=ax)

        return fig

    def to_pymatgen(self, qlabels=None):
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
            # for q, phf in zip(split_q, split_phf)[1:-1]:
            for i, (q, phf, d) in enumerate(zip(split_q, split_phf, split_phdispl)):
                ph_freqs.append(phf)
                qpts.append(q)
                d = d.reshape(self.num_branches, self.num_atoms, 3)
                displ.append(d)
                # if the qpoint has a label it nees to be repeated. If it is one of the extrama either it should
                # not be repeated (if they are the real first or last point) or they will be already reapeated due
                # to the split.
                if any(np.allclose(q, labelled_q) for labelled_q in labelled_q_list):
                    if 0 < i <len(split_q) - 1:
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

    def acoustic_indices(self, qpoint, threshold=0.95, raise_on_no_indices=True):
        """
        Extract the indices of the three acoustic modes for a qpoint.
        Acoustic modes could be reasonably identified for Gamma and points close to Gamma.

        Args:
            qpoint: the qpoint. Accepts integer or reduced coordinates
            threshold: fractional value allowed for the matching of the displacements to identify acoustic modes.
            raise_on_no_indices: if True a RuntimeError will be raised if the acoustic mode will not be
                correctly identified
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

        if len(indices) != 3 and raise_on_no_indices:
            raise RuntimeError('wrong number of indices: {}'.format(indices))
        else:
            indices = [0, 1, 2]

        return indices

    def asr_breaking(self, units='eV', threshold=0.95, raise_on_no_indices=True):
        """
        Calculates the breaking of the acoustic sum rule.
        Requires the presence of Gamma.

        Args:
            units: Units for the output. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            threshold: fractional value allowed for the matching of the displacements to identify acoustic modes.
            raise_on_no_indices: if True a RuntimeError will be raised if the acoustic mode will not be
                correctly identified

        Returns:
            A namedtuple with:
                the three breaking of the acoustic modes
                the maximum breaking with sign
                the absolute value of the maximum breaking
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
        displ = self.phdispl_cart[qind, nmode].reshape((-1,3))

        return self.structure.frozen_phonon(qpoint=self.qpoints[qind].frac_coords, displ=displ, eta=eta,
                                            frac_coords=False, scale_matrix=scale_matrix, max_supercell=max_supercell)


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

    def __init__(self, filepath):
        """
        Args:
            path: path to the file
        """
        super(PhbstFile, self).__init__(filepath)
        self.reader = PHBST_Reader(filepath)

        # Initialize Phonon bands and add metadata from ncfile
        self._phbands = PhononBands.from_file(filepath)

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

        app(self.phbands.to_string(title=None, with_structure=True, with_qpoints=False, verbose=verbose))

        return "\n".join(lines)

    @property
    def structure(self):
        """|Structure| object"""
        return self.phbands.structure

    @property
    def qpoints(self):
        """List of q-point objects."""
        return self.phbands.qpoints

    @property
    def phbands(self):
        """|PhononBands| object"""
        return self._phbands

    def close(self):
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self):
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
        import pandas as pd
        frame = pd.DataFrame(d, columns=list(d.keys()))
        frame.qpoint = qpoint

        return frame

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


class PhononDos(Function1D):
    """
    This object stores the phonon density of states.
    An instance of ``PhononDos`` has a ``mesh`` (numpy array with the points of the mesh)
    and another numpy array, ``values``, with the DOS on the mesh.

    .. note::

        mesh is given in eV, values are in states/eV.
    """
    #def __init__(self, mesh, values, qmesh):
    #    super(PhononDos, self).__init__(mesh, values)
    #    self.qmesh = qmesh

    @classmethod
    def as_phdos(cls, obj, phdos_kwargs=None):
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

        raise TypeError("Don't know how to create `PhononDos` from %s" % type(obj))

    @lazy_property
    def iw0(self):
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
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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

    # TODO: This should be called plot_dos_idos!
    @add_fig_kwargs
    def plot(self, units="eV", **kwargs):
        """
        Plot Phonon DOS and IDOS on two distict plots.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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

    def get_internal_energy(self, tstart=5, tstop=300, num=50):
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

    def get_entropy(self, tstart=5, tstop=300, num=50):
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

    def get_free_energy(self, tstart=5, tstop=300, num=50):
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

    def get_cv(self, tstart=5, tstop=300, num=50):
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
                             quantities=None, fontsize=8, **kwargs):
        """
        Plot thermodynamic properties from the phonon DOSes within the harmonic approximation.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 50.
            quantities: List of strings specifying the thermodynamic quantities to plot.
                Possible values: ["internal_energy", "free_energy", "entropy", "c_v"].
                None means all.
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            formula_units: the number of formula units per unit cell. If unspecified, the
                thermodynamic quantities will be given on a per-unit-cell basis.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        quantities = list_strings(quantities) if quantities is not None else \
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
            # Compute thermodynamic quantity associated to qname.
            f1d = getattr(self, "get_" + qname)(tstart=tstart, tstop=tstop, num=num)
            ys = f1d.values
            if formula_units is not None: ys /= formula_units
            if units == "Jmol": ys = ys * abu.e_Cb * abu.Avogadro
            ax.plot(f1d.mesh, ys)

            ax.set_title(qname)
            ax.grid(True)
            ax.set_xlabel("Temperature (K)", fontsize=fontsize)
            ax.set_ylabel(_THERMO_YLABELS[qname][units], fontsize=fontsize)
            #ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def to_pymatgen(self):
        """
        Creates a pymatgen :class:`PmgPhononDos` object
        """
        factor = abu.phfactor_ev2units("thz")

        return PmgPhononDos(self.mesh*factor, self.values/factor)

    @property
    def debye_temp(self):
        """
        Debye temperature in K.
        """
        integrals = (self * self.mesh ** 2).spline_integral() / self.spline_integral()
        t_d = np.sqrt(5/3*integrals)/abu.kb_eVK

        return t_d

    def get_acoustic_debye_temp(self, nsites):
        """
        Acoustic Debye temperature in K, i.e. the Debye temperature divided by nsites**(1/3).

        Args:
            nsites: the number of sites in the cell.
        """
        return self.debye_temp/nsites**(1/3)


class PhdosReader(ETSF_Reader):
    """
    This object reads data from the PHDOS.nc file produced by anaddb.

    .. note::

            Frequencies are in eV, DOSes are in states/eV.
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

    def read_phdos(self):
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
        # [ntypat, nomega] array with PH-DOS projected over atom types."""
        values = self.read_pjdos_type()

        od = OrderedDict()
        for symbol in self.chemical_symbols:
            type_idx = self.typeidx_from_symbol(symbol)
            od[symbol] = PhononDos(self.wmesh, values[type_idx])

        return od

    # TODO
    #double msqd_dos_atom(number_of_atoms, three, three, number_of_frequencies)


class PhdosFile(AbinitNcFile, Has_Structure, NotebookWriter):
    """
    Container object storing the different DOSes stored in the
    PHDOS.nc file produced by anaddb.
    Provides helper function to visualize/extract data.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PhdosFile
    """

    def __init__(self, filepath):
        # Open the file, read data and create objects.
        super(PhdosFile, self).__init__(filepath)

        self.reader = r = PhdosReader(filepath)
        self.wmesh = r.wmesh

    def close(self):
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self):
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

    def to_string(self, verbose=0):
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
    def structure(self):
        """|Structure| object."""
        return self.reader.structure

    @lazy_property
    def phdos(self):
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

    @add_fig_kwargs
    def plot_pjdos_type(self, units="eV", stacked=True, colormap="jet", alpha=0.7, exchange_xy=False,
                        ax=None, xlims=None, ylims=None, fontsize=12, **kwargs):
        """
        Plot type-projected phonon DOS.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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

    @add_fig_kwargs
    def plot_pjdos_cartdirs_type(self, units="eV", stacked=True, colormap="jet", alpha=0.7,
                                 xlims=None, ylims=None, ax_list=None, fontsize=8, **kwargs):
        """
        Plot type-projected phonon DOS decomposed along the three cartesian directions.
        Three rows for each cartesian direction. Each row shows the contribution of each atomic type + Total Phonon DOS.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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
                                 xlims=None, ylims=None, ax_list=None, fontsize=8, **kwargs):
        """
        Plot phonon PJDOS for each atom in the unit cell. By default, only "inequivalent" atoms are shown.

        Args:
            view: "inequivalent", "all"
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            colormap: matplotlib colormap.
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        # Define num_plots and ax2atom depending on view.
        factor = abu.phfactor_ev2units(units)
        natom, ntypat = len(self.structure), self.structure.ntypesp
        lw = kwargs.pop("lw", 2)

        if view == "all" or natom == 1:
            iatom_list = np.arange(natom)

        elif view == "inequivalent":
            print("Calling spglib to find inequivalent sites.")
            print("Note that `symafm` magnetic symmetries (if any) are not taken into account.")
            ea = self.structure.spget_equivalent_atoms(printout=True)
            iatom_list = ea.irred_pos
        else:
            raise ValueError("Wrong value for view: %s" % str(view))

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
            for iatom in iatom_list:
                site = self.structure[iatom]
                symbol = str(site)
                color = cmap(float(iatom) / (len(iatom_list) - 1))
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
        ])

        return self._write_nb_nbpath(nb, nbpath)

    def to_pymatgen(self):
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
        units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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


def dataframe_from_phbands(phbands_objects, index=None, with_spglib=True):
    """
    Build a pandas dataframe with the most important results available in a list of band structures.

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

    import pandas as pd
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
    _LINE_COLORS = ["b", "r", "g", "m", "y", "k"]
    _LINE_STYLES = ["-", ":", "--", "-.",]
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

    def __repr__(self):
        """Invoked by repr"""
        return self.to_string(func=repr)

    def __str__(self):
        """Invoked by str"""
        return self.to_string(func=str)

    def add_plotter(self, other):
        """Merge two plotters, return new plotter."""
        if not isinstance(other, self.__class__):
            raise TypeError("Don't know to to add %s to %s" % (other.__class__, self.__class__))

        key_phbands = list(self._bands_dict.items()) + list(other._bands_dict.items())
        key_phdos = list(self._phdoses_dict.items()) + list(other._phdoses_dict.items())

        return self.__class__(key_phbands=key_phbands, key_phdos=key_phdos)

    def to_string(self, func=str, verbose=0):
        """String representation."""
        lines = []
        app = lines.append
        for i, (label, phbands) in enumerate(self.phbands_dict.items()):
            app("[%d] %s --> %s" % (i, label, func(phbands)))

        if self.phdoses_dict:
            for i, (label, phdos) in enumerate(self.phdoses_dict.items()):
                app("[%d] %s --> %s" % (i, label, func(phdos)))

        return "\n".join(lines)

    def has_same_formula(self):
        """
        True of plotter contains structures with same chemical formula.
        """
        structures = [phbands.structure for phbands in self.phbands_dict.values()]
        if structures and any(s.formula != structures[0].formula for s in structures): return False
        return True

    def get_phbands_frame(self, with_spglib=True):
        """
        Build a |pandas-DataFrame| with the most important results available in the band structures.
        """
        return dataframe_from_phbands(list(self.phbands_dict.values()),
                                      index=list(self.phbands_dict.keys()), with_spglib=with_spglib)

    @property
    def phbands_dict(self):
        """Dictionary with the mapping label --> phbands."""
        return self._bands_dict

    # TODO: Just an alias. To be removed in 0.4
    bands_dict = phbands_dict

    @property
    def phdoses_dict(self):
        """Dictionary with the mapping label --> phdos."""
        return self._phdoses_dict

    @property
    def phbands_list(self):
        """"List of |PhononBands| objects."""
        return list(self._bands_dict.values())

    @property
    def phdoses_list(self):
        """"List of |PhononDos|."""
        return list(self._phdoses_dict.values())

    def iter_lineopt(self):
        """Generates matplotlib linestyles."""
        for o in itertools.product(self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def add_phbands(self, label, bands, phdos=None, dos=None, phdos_kwargs=None):
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

    #def bands_statdiff(self, ref=0):
    #    """
    #    Compare the reference bands with index ref with the other bands stored in the plotter.
    #    """
    #    for i, label in enumerate(self._bands_dict.keys()):
    #        if i == ref:
    #            ref_label = label
    #            break
    #    else:
    #        raise ValueError("ref index %s is > number of bands" % ref)

    #    ref_bands = self._bands_dict[ref_label]

    #    text = []
    #    for label, bands in self._bands_dict.items():
    #        if label == ref_label: continue
    #        stat = ref_bands.statdiff(bands)
    #        text.append(str(stat))

    #    return "\n\n".join(text)

    @add_fig_kwargs
    def combiplot(self, qlabels=None, units='eV', ylims=None, width_ratios=(2, 1), fontsize=8,
                  linestyle_dict=None, **kwargs):
        r"""
        Plot the band structure and the DOS on the same figure.
        Use ``gridplot`` to plot band structures on different figures.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the k-points.
                The values are the labels e.g. ``klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            width_ratios: Ratio between the width of the phonon bands plots and the DOS plots.
                Used if plotter has DOSes.
            fontsize: fontsize for titles and legend.
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
            cprint("WARNING: Bands have different number of k-points:\n%s" % str(nqpt_list), "yellow")

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

    @add_fig_kwargs
    def gridplot(self, with_dos=True, units="eV", fontsize=8, **kwargs):
        """
        Plot multiple phonon bandstructures and optionally DOSes on a grid.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            with_dos: True to plot phonon DOS (if available).
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        titles = list(self._bands_dict.keys())
        phb_objects = list(self._bands_dict.values())
        phdos_objects = None
        if self.phdoses_dict and with_dos:
            phdos_objects = list(self.phdoses_dict.values())

        return phbands_gridplot(phb_objects, titles=titles, phdos_objects=phdos_objects, units=units, fontsize=fontsize, show=False)

    @add_fig_kwargs
    def gridplot_with_hue(self, hue, with_dos=False, units="eV", width_ratios=(2, 1), ylims=None, fontsize=8, **kwargs):
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
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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

        from abipy.tools import sort_and_groupby, getattrd, hasattrd

        # Need index to handle all_phdos_objects if DOSes are wanted.
        if callable(hue):
            items = [(hue(phb), phb, i, label) for i, (phb, label) in enumerate(zip(all_phb_objects, all_labels))]
        else:
            # Assume string. Either phbands.hue or phbands.params[hue].
            if hasattrd(all_phb_objects[0], hue):
                items = [(getattrd(phb, hue), phb, i, label) for i, (phb, label) in enumerate(zip(all_phb_objects, all_labels))]
            else:
                items = [(phb.params[hue], phb, i, label) for i, (phb, label) in enumerate(zip(all_phb_objects, all_labels))]

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
        Use seaborn_ to draw a box plot to show distributions of eigenvalues with respect to the band index.
        Band structures are drawn on different subplots.

        Args:
            mode_range: Only bands such as ``mode_range[0] <= nu_index < mode_range[1]`` are included in the plot.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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
        Use seaborn_ to draw a box plot comparing the distributions of the frequencies.
        Phonon Band structures are drawn on the same plot.

        Args:
            mode_range: Only bands such as ``mode_range[0] <= nu_index < mode_range[1]`` are included in the plot.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            swarm: True to show the datapoints on top of the boxes
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: Keyword arguments passed to seaborn_ boxplot.
        """
        frames = []
        for label, phbands in self.phbands_dict.items():
            # Get the dataframe, select bands and add column with label
            frame = phbands.get_dataframe()
            if mode_range is not None:
                frame = frame[(frame["mode"] >= mode_range[0]) & (frame["mode"] < mode_range[1])]
            frame["label"] = label
            frames.append(frame)

        # Merge frames ignoring index (not meaningful here)
        import pandas as pd
        data = pd.concat(frames, ignore_index=True)

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
        at a given q-point. One panel for all phbands stored in the plotter.

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
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            width_ratios: Ratio between the band structure plot and the dos plot.
                Used when there are DOS stored in the plotter.
            show: True if the animation should be shown immediately

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
        fig = plotter.gridplot()
    """
    def __init__(self, key_phdos=None, phdos_kwargs=None):
        self._phdoses_dict = OrderedDict()
        if key_phdos is None: key_phdos = []
        for label, phdos in key_phdos:
            self.add_phdos(label, phdos, phdos_kwargs=phdos_kwargs)

    @property
    def phdos_list(self):
        """List of phonon DOSes"""
        return list(self._phdoses_dict.values())

    def add_phdos(self, label, phdos, phdos_kwargs=None):
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
    #    True of plotter contains structures with same chemical formula.
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
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
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

    def plot(self, **kwargs):
        """An alias for combiplot."""
        return self.combiplot(**kwargs)

    @add_fig_kwargs
    def gridplot(self, units="eV", xlims=None, ylims=None, fontsize=8, **kwargs):
        """
        Plot multiple DOSes on a grid.

        Args:
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        titles = list(self._phdoses_dict.keys())
        phdos_list = list(self._phdoses_dict.values())

        import matplotlib.pyplot as plt
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

    @add_fig_kwargs
    def plot_harmonic_thermo(self, tstart=5, tstop=300, num=50, units="eV", formula_units=1,
                             quantities="all", fontsize=8, **kwargs):
        """
        Plot thermodynamic properties from the phonon DOS within the harmonic approximation.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 50.
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            formula_units: the number of formula units per unit cell. If unspecified, the
                thermodynamic quantities will be given on a per-unit-cell basis.
            quantities: List of strings specifying the thermodynamic quantities to plot.
                Possible values in ["internal_energy", "free_energy", "entropy", "c_v"].
            fontsize: Legend and title fontsize.

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
            ax.set_xlabel("Temperature (K)", fontsize=fontsize)
            if iax == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        #fig.tight_layout()
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


class NonAnalyticalPh(Has_Structure):
    """
    Phonon data at gamma including non analytical contributions
    Read from anaddb.nc
    """

    def __init__(self, structure, directions, phfreqs, phdispl_cart, amu=None):
        """
        Args:
            structure: |Structure| object.
            directions: Cartesian directions along which the non analytical frequencies have been calculated
            phfreqs: Phonon frequencies with non analytical contribution in eV along directions
            phdispl_cart: Displacement in Angstrom in Cartesian coordinates with non analytical contribution
                along directions
            amu: dictionary that associates the atomic species present in the structure to the values of the atomic
                mass units used for the calculation
        """
        self._structure = structure
        self.directions = directions
        self.phfreqs = phfreqs
        self.phdispl_cart = phdispl_cart
        self.amu = amu
        self.amu_symbol = None
        if amu is not None:
            self.amu_symbol = {}
            for z, m in amu.items():
                el = Element.from_Z(int(z))
                self.amu_symbol[el.symbol] = m

    @classmethod
    def from_file(cls, filepath):
        """
        Reads the non analytical directions, frequencies and displacements from the anaddb.nc file specified.
        Non existence of displacements is accepted for compatibility with abinit 8.0.6
        Raises an error if the other values are not present in anaddb.nc.
        """
        with ETSF_Reader(filepath) as r:
            directions = r.read_value("non_analytical_directions")
            phfreq = r.read_value("non_analytical_phonon_modes")

            # need a default as the first abinit version including IFCs in the netcdf doesn't have this attribute
            phdispl_cart = r.read_value("non_analytical_phdispl_cart", cmode="c", default=None)

            structure = r.read_structure()

            amu_list = r.read_value("atomic_mass_units", default=None)
            if amu_list is not None:
                # ntypat arrays
                atomic_numbers = r.read_value("atomic_numbers")
                amu = {at: a for at, a in zip(atomic_numbers, amu_list)}
            else:
                amu = None

            return cls(structure=structure, directions=directions, phfreqs=phfreq, phdispl_cart=phdispl_cart, amu=amu)

    @lazy_property
    def dyn_mat_eigenvect(self):
        """
        [ndirection, 3*natom, 3*natom] array with the orthonormal eigenvectors of the dynamical matrix.
        in Cartesian coordinates.
        """
        return get_dyn_mat_eigenvec(self.phdispl_cart, self.structure, amu=self.amu)

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def index_direction(self, direction, cartesian=False):
        """
        Returns: the index of direction. Raises: `ValueError` if not found.

        Args:
            direction: a 3 element list indicating the direction. Can be a generic vector
            cartesian: if True the direction are already in cartesian coordinates, if False it
                will be converted to match the internal description of the directions.
        """
        if not cartesian:
            direction = self.structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(direction)
        else:
            direction = np.array(direction)
        direction = direction / np.linalg.norm(direction)

        for i, d in enumerate(self.directions):
            d = d / np.linalg.norm(d)
            if np.allclose(d, direction):
                return i

        raise ValueError("Cannot find direction: `%s` with cartesian: `%s` in non_analytical cartesian directions:\n%s" %
                (str(direction), cartesian, str(self.directions)))

    def has_direction(self, direction, cartesian=False):
        """
        Checks if the input direction is among those available.

        Args:
            direction: a 3 element list indicating the direction. Can be a generic vector
            cartesian: if True the direction are already in cartesian coordinates, if False it
                will be converted to match the internal description of the directions.
        """
        try:
            self.index_direction(direction, cartesian=cartesian)
            return True
        except ValueError:
            return False


class InteratomicForceConstants(Has_Structure):
    """
    The interatomic force constants calculated by anaddb.
    Read from anaddb.nc
    """

    def __init__(self, structure, atoms_indices, neighbours_indices, ifc_cart_coord,
                 ifc_cart_coord_short_range, local_vectors, distances):
        """
        Args:
            structure: |Structure| object.
            atoms_index: List of integers representing the indices in the structure of the analyzed atoms.
            neighbours_index: List of integers representing the indices in the structure of the neighbour atoms.
            ifc_cart_coord: ifc in Cartesian coordinates
            ifc_cart_coord_short_range: short range part of the ifc in Cartesian coordinates
            local_vectors: local basis used to determine the ifc_local_coord
        """
        self._structure = structure
        self.atoms_indices = atoms_indices
        self.neighbours_indices = neighbours_indices
        self.ifc_cart_coord = ifc_cart_coord
        self.ifc_cart_coord_short_range = ifc_cart_coord_short_range
        self.local_vectors = local_vectors
        self.distances = distances

    @property
    def number_of_atoms(self):
        """Number of atoms is structure."""
        return len(self.structure)

    @classmethod
    def from_file(cls, filepath):
        """Create the object from a netcdf_ file."""
        with ETSF_Reader(filepath) as r:
            try:
                structure = r.read_structure()
                atoms_indices = r.read_value("ifc_atoms_indices") - 1
                neighbours_indices = r.read_value("ifc_neighbours_indices") - 1
                distances = r.read_value("ifc_distances")
                ifc_cart_coord = r.read_value("ifc_matrix_cart_coord")
                ifc_cart_coord_short_range = r.read_value("ifc_matrix_cart_coord_short_range", default=None)
                local_vectors = r.read_value("ifc_local_vectors", default=None)
            except:
                import traceback
                msg = traceback.format_exc()
                msg += ("Error while trying to read IFCs from file.\n"
                       "Verify that the required variables are used in anaddb: ifcflag, natifc, atifc, ifcout\n")
                raise ValueError(msg)

            return cls(structure=structure, atoms_indices=atoms_indices, neighbours_indices=neighbours_indices,
                       ifc_cart_coord=ifc_cart_coord,
                       ifc_cart_coord_short_range=ifc_cart_coord_short_range, local_vectors=local_vectors,
                       distances=distances)

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    @property
    def number_of_neighbours(self):
        """Number of neighbouring atoms for which the ifc are present. ifcout in anaddb."""
        return np.shape(self.neighbours_indices)[1]

    @lazy_property
    def ifc_cart_coord_ewald(self):
        """Ewald part of the ifcs in cartesian coordinates"""
        if self.ifc_cart_coord_short_range is None:
            return None
        else:
            return self.ifc_cart_coord-self.ifc_cart_coord_short_range

    @lazy_property
    def ifc_local_coord(self):
        """Ifcs in local coordinates"""
        if self.local_vectors is None:
            return None
        else:
            return np.einsum("ktli,ktij,ktuj->ktlu", self.local_vectors, self.ifc_cart_coord, self.local_vectors)

    @lazy_property
    def ifc_local_coord_short_range(self):
        """Short range part of the ifcs in cartesian coordinates"""
        if self.local_vectors is None:
            return None
        else:
            return np.einsum("ktli,ktij,ktuj->ktlu", self.local_vectors, self.ifc_cart_coord_short_range, self.local_vectors)

    @lazy_property
    def ifc_local_coord_ewald(self):
        """Ewald part of the ifcs in local coordinates"""
        return np.einsum("ktli,ktij,ktuj->ktlu", self.local_vectors, self.ifc_cart_coord_ewald, self.local_vectors)

    def _filter_ifc_indices(self, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None, max_dist=None):
        """
        Internal method that provides the indices of the neighouring atoms in self.neighbours_indices that satisfy
        the required conditions. All the arguments are optional. If None the filter will not be applied.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
        """

        if atom_indices is not None and atom_element is not None:
            raise ValueError("atom_index and atom_element cannot be specified simultaneously")

        if atom_indices is not None and not isinstance(atom_indices, (list, tuple)):
            atom_indices = [atom_indices]

        if atom_element:
            atom_indices = self.structure.indices_from_symbol(atom_element)

        if atom_indices is None:
            atom_indices = range(len(self.structure))

        # apply the filter: construct matrices of num_atoms*num_neighbours size, all conditions should be satisfied.
        ind = np.where(
            (np.tile(np.in1d(self.atoms_indices, atom_indices), [self.number_of_neighbours, 1])).T &
            (self.distances > min_dist if min_dist is not None else True) &
            (self.distances < max_dist if max_dist is not None else True) &
            (np.in1d(self.neighbours_indices, self.structure.indices_from_symbol(neighbour_element))
             .reshape(self.number_of_atoms, self.number_of_neighbours) if neighbour_element is not None else True)
        )

        return ind

    def get_ifc_cartesian(self, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None, max_dist=None):
        """
        Filters the IFCs in cartesian coordinates
        All the arguments are optional. If None the filter will not be applied.
        Returns two arrays containing the distances and the corresponding filtered ifcs.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
        """
        ind = self._filter_ifc_indices(atom_indices=atom_indices, atom_element=atom_element,
                                       neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist)

        return self.distances[ind], self.ifc_cart_coord[ind]

    def get_ifc_local(self, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None, max_dist=None):
        """
        Filters the IFCs in local coordinates
        All the arguments are optional. If None the filter will not be applied.
        Returns two arrays containing the distances and the corresponding filtered ifcs.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        ind = self._filter_ifc_indices(atom_indices=atom_indices, atom_element=atom_element,
                                       neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist)

        return self.distances[ind], self.ifc_local_coord[ind]

    def get_plot_ifc(self, ifc, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None,
                     max_dist=None, ax=None, **kwargs):
        """
        Plots the specified ifcs, filtered according to the optional arguments.
        An array with shape number_of_atoms*number_of_neighbours, so only one of the components of the ifc matrix can
        be plotted at a time.

        Args:
            ifc: an array with shape number_of_atoms * number_of_neighbours of the ifc that should be plotted
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ind = self._filter_ifc_indices(atom_indices=atom_indices, atom_element=atom_element,
                                       neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist)

        dist, filtered_ifc =  self.distances[ind], ifc[ind]

        if 'color' not in kwargs:
            kwargs['color'] = 'blue'

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['lw'] = 0

        ax.set_xlabel('Distance (Bohr)')
        ax.set_ylabel(r'IFC (Ha/Bohr$^2$)')
        ax.grid(True)

        ax.plot(dist, filtered_ifc, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_longitudinal_ifc(self, atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None,
                              max_dist=None, ax=None, **kwargs):
        """
        Plots the total longitudinal ifcs in local coordinates, filtered according to the optional arguments.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns: |matplotlib-Figure|
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        return self.get_plot_ifc(self.ifc_local_coord[:, :, 0, 0], atom_indices=atom_indices, atom_element=atom_element,
                                 neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist, ax=ax, **kwargs)

    @add_fig_kwargs
    def plot_longitudinal_ifc_short_range(self, atom_indices=None, atom_element=None, neighbour_element=None,
                                          min_dist=None, max_dist=None, ax=None, **kwargs):
        """
        Plots the short range longitudinal ifcs in local coordinates, filtered according to the optional arguments.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns: |matplotlib-Figure|
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        if self.ifc_local_coord_short_range is None:
            raise ValueError("Ewald contribution is missing, Run anaddb with dipdip=1")

        return self.get_plot_ifc(self.ifc_local_coord_short_range[:, :, 0, 0], atom_indices=atom_indices,
                                 atom_element=atom_element, neighbour_element=neighbour_element, min_dist=min_dist,
                                 max_dist=max_dist, ax=ax, **kwargs)

    @add_fig_kwargs
    def plot_longitudinal_ifc_ewald(self, atom_indices=None, atom_element=None, neighbour_element=None,
                                          min_dist=None, max_dist=None, ax=None, **kwargs):
        """
        Plots the Ewald part of the ifcs in local coordinates, filtered according to the optional arguments.

        Args:
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns: |matplotlib-Figure|
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        if self.ifc_local_coord_ewald is None:
            raise ValueError("Ewald contribution is missing, Run anaddb with dipdip=1")

        return self.get_plot_ifc(self.ifc_local_coord_ewald[:, :, 0, 0], atom_indices=atom_indices,
                                 atom_element=atom_element, neighbour_element=neighbour_element, min_dist=min_dist,
                                 max_dist=max_dist, ax=ax, **kwargs)


# TODO: amu should become mandatory.
def get_dyn_mat_eigenvec(phdispl, structure, amu=None, amu_symbol=None):
    """
    Converts the phonon displacements to the orthonormal eigenvectors of the dynamical matrix.
    Small discrepancies with the original values may be expected due to the different values of the atomic masses in
    abinit and pymatgen.

    .. note::

        These eigenvectors are orthonormalized and should be very close to the ones computed by Abinit in a.u.
        Note, however, that the output vectors are given in atomic units so dividing then by the sqrt(Mass)
        won't give the dipl_cart used in PhononBands that are in Angstrom.

    Args:
        phdispl: a numpy array containing the displacements in cartesian coordinates. The last index should have
            size 3*(num atoms), but the rest of the shape is arbitrary. If qpts is not None the first dimension
            should match the q points.
        structure: |Structure| object.
        amu: dictionary that associates the atomic numbers present in the structure to the values of the atomic
            mass units used for the calculation. Incompatible with amu_sumbol. If None and amu_symbol is None, values
            from pymatgen will be used.  Note that this will almost always lead to inaccuracies in the conversion.
        amu_symbol: dictionary that associates the symbol present in the structure to the values of the atomic
            mass units used for the calculation. Incompatible with amu. If None and amu_symbol is None, values from
            pymatgen will be used. that this will almost always lead to inaccuracies in the conversion.

    Returns:
        A |numpy-array| of the same shape as phdispl containing the eigenvectors of the dynamical matrix
    """
    eigvec = np.zeros(np.shape(phdispl), dtype=np.complex)

    if amu is not None and amu_symbol is not None:
        raise ValueError("Only one between amu and amu_symbol should be provided!")

    if amu is not None:
        amu_symbol = {Element.from_Z(n).symbol: v for n, v in amu.items()}

    if amu_symbol is None:
        warnings.warn("get_dyn_mat_eigenvec has been called with amu=None. Eigenvectors may not be orthonormal.")
        amu_symbol = {e.symbol: e.atomic_mass for e in structure.composition.elements}

    for j, a in enumerate(structure):
        eigvec[...,3*j:3*(j+1)] = phdispl[...,3*j:3*(j+1)] * np.sqrt(amu_symbol[a.specie.symbol]*abu.amu_emass) / abu.Bohr_Ang

    return eigvec


def match_eigenvectors(v1, v2):
    """
    Given two list of vectors, returns the pair matching based on the complex scalar product.
    Returns the indices of the second list that match the vectors of the first list in ascending order.
    """
    prod = np.absolute(np.dot(v1, v2.transpose().conjugate()))

    indices = np.zeros(len(v1), dtype=np.int)
    missing_v1 = [True] * len(v1)
    missing_v2 = [True] * len(v1)
    for m in reversed(np.argsort(prod, axis=None)):
        i, j = np.unravel_index(m, prod.shape)
        if missing_v1[i] and missing_v2[j]:
            indices[i] = j
            missing_v1[i] = missing_v2[j] = False
            if not any(missing_v1):
                if any(missing_v2):
                    raise RuntimeError('Something went wrong in matching vectors: {} {}'.format(v1, v2))
                break

    return indices


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
        """Wraps combiboxplot method of |ElectronDosPlotter|. kwargs passed to combiboxplot."""
        return self.get_phbands_plotter().combiboxplot(**kwargs)

    #def combiplot_phdos(self, **kwargs):
    #    """Wraps combiplot method of |ElectronDosPlotter|. kwargs passed to combiplot."""
    #    return self.get_phdos_plotter().combiplot(**kwargs)
    #
    #def gridplot_phdos(self, **kwargs):
    #    """Wraps gridplot method of |ElectronDosPlotter|. kwargs passed to gridplot."""
    #    return self.get_phdos_plotter().gridplot(**kwargs)

    def get_phbands_plotter(self, filter_abifile=None, cls=None):
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

    def get_phbands_dataframe(self, with_spglib=True):
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


def open_file_phononwebsite(filename, port=8000,
                            website="http://henriquemiranda.github.io/phononwebsite",
                            host="localhost", browser=None): # pragma: no cover
    """
    Take a file, detect the type and open it on the phonon website
    Based on a similar function in <https://github.com/henriquemiranda/phononwebsite/phononweb.py>

    Args:
        filename: file with phonon data in phononwebsite format.
        port: Initial port.
        website: Website URL
        host: localhost name.
        browser: Open webpage in ``browser``. Use default if $BROWSER if None.
    """
    if filename.endswith(".json"):
        filetype = "json"
    elif filename.endswith(".yaml"):
        filetype = "yaml"
    else:
        filetype = "rest"

    try:
        from http.server import HTTPServer, SimpleHTTPRequestHandler
    except ImportError:
        from BaseHTTPServer import HTTPServer
        # python 2 requires internal implementation
        from abipy.tools.SimpleHTTPServer import SimpleHTTPRequestHandler

    # Add CORS header to the website
    class CORSRequestHandler (SimpleHTTPRequestHandler):
        def end_headers (self):
            #self.send_header('Access-Control-Allow-Origin', website)
            self.send_header('Access-Control-Allow-Origin', "http://henriquemiranda.github.io")
            SimpleHTTPRequestHandler.end_headers(self)
        def log_message(self, format, *args):
            return

    # Initialize http server thread
    print('Starting HTTP server at port %d ...' % port, end=" ")
    trial, max_ntrial = 0, 50
    while trial < max_ntrial:
        try:
            server = HTTPServer(('', port), CORSRequestHandler)
            #print("got port:", port)
            break
        except OSError:
            trial += 1
            port += 1
            print(port, end=", ")
    else:
        raise RuntimeError("Cannot find available port after %s attempts" % max_ntrial)

    # Create threads python
    server.url = 'http://{}:{}'.format(host, server.server_port)
    from threading import Thread
    t = Thread(target=server.serve_forever)
    t.daemon = True
    t.start()

    # Open website with the file
    try:
        from urllib.parse import quote
    except ImportError:
        from urllib import quote

    url_filename = 'http://{}:{}/{}'.format(host, server.server_port, quote(filename))
    url = '%s/phonon.html?%s=%s' % (website, filetype, url_filename)
    print("\nOpening URL:", url)
    print("Using default browser, if the webpage is not displayed correctly",
          "\ntry to change browser either via command line options or directly in the shell with e.g:\n\n"
          "     export BROWSER=firefox\n")
    print('Press Ctrl+C to terminate HTTP server')
    import webbrowser
    webbrowser.get(browser).open_new_tab(url)

    # Quit application when SIGINT is received
    def signal_handler(signal, frame):
        sys.exit(0)

    import signal
    signal.signal(signal.SIGINT, signal_handler)
    signal.pause()