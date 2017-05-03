# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import functools
import numpy as np
import itertools
import pickle
import os
import six
import json
import abipy.core.abinit_units as abu

from collections import OrderedDict
from monty.string import is_string, list_strings, marquee
from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.termcolor import cprint
from monty.dev import deprecated
from pymatgen.core.units import eV_to_Ha, Energy
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from abipy.core.func1d import Function1D
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_PhononBands, NotebookWriter
from abipy.core.kpoints import Kpoint, KpointList
from abipy.iotools import ETSF_Reader
from abipy.tools import gaussian, duck
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims



__all__ = [
    "frame_from_phbands",
    "PhononBands",
    "PhononBandsPlotter",
    "PhbstFile",
    "PhononDos",
    "PhononDosPlotter",
    "PhdosReader",
    "PhdosFile",
]

def _factor_ev2units(units):
    """
    Return conversion factor eV --> units (case-insensitive)
    """
    eV_to_cm1 = 8065.5440044136285
    d = {"ev": 1, "mev": 1000, "ha": eV_to_Ha,
         "cm-1": eV_to_cm1, 'cm^-1': eV_to_cm1, "thz": abu.eV_to_THz,
         }
    try:
        return d[units.lower()]
    except KeyError:
        raise KeyError('Value for units `{}` unknown\nPossible values are:\n {}'.format(units, list(d.keys())))


def _unit_tag(units):
    d = {"ev": "[eV]", "mev": "[meV]", "ha": '[Ha]',
         "cm-1": "[cm$^{-1}$]", 'cm^-1': "[cm$^{-1}$]", "thz": '[Thz]',
         }
    try:
        return d[units.lower()]
    except KeyError:
        raise KeyError('Value for units `{}` unknown\nPossible values are:\n {}'.format(units, list(d.keys())))


def _dos_label_from_units(units):
    d = {"ev": "[states/eV]", "mev": "[states/meV]", "ha": '[states/Ha]',
         "cm-1": "[states/cm$^{-1}$]", 'cm^-1': "[states/cm$^{-1}$]",
         "thz": '[states/Thz]',
         }
    try:
        return d[units.lower()]
    except KeyError:
        raise KeyError('Value for units `{}` unknown\nPossible values are:\n {}'.format(units, list(d.keys())))


@functools.total_ordering
class PhononMode(object):
    """
    A phonon mode has a q-point, a frequency, a cartesian displacement and a structure.
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
            displ: Displacement (Cartesian coordinates, Angstrom)
            structure: :class:`Structure` object.
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

    def to_string(self, with_displ=True):
        """
        String representation

        Args:
           with_displ: True to print phonon displacement.
	"""
        lines = ["%s: q-point %s, frequency %.5f [eV]" % (self.__class__.__name__, self.qpoint, self.freq)]
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
        """Create the object from a netCDF file."""
        with PHBST_Reader(filepath) as r:
            structure = r.read_structure()

            # Build the list of q-points
            qpoints = KpointList(structure.reciprocal_lattice,
                                 frac_coords=r.read_qredcoords(),
                                 weights=r.read_qweights(),
                                 names=None)

            # Read amu
            amu_list = r.read_amu()
            if amu_list is not None:
                atom_species = r.read_value("atomic_numbers")
                amu = {at: a for at, a in zip(atom_species, amu_list)}
            else:
                cprint("Warning: file %s does not contain atomic_numbers.\nParticular methods need them!" %
                        filepath, "red")
                amu = None

            return cls(structure=structure,
                       qpoints=qpoints,
                       phfreqs=r.read_phfreqs(),
                       phdispl_cart=r.read_phdispl_cart(),
                       amu=amu)

    @classmethod
    def as_phbands(cls, obj):
        """
        Return an instance of :class:`PhononBands` from a generic obj.
        Supports:

            - instances of cls
            - files (string) that can be open with abiopen and that provide a `phbands` attribute.
            - objects providing a `phbands` attribute.
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
    def factor_ev2units(units):
        """
        Return conversion factor eV --> units (case-insensitive)
        """
        return _factor_ev2units(units)

    def read_non_anal_from_file(self, filepath):
        """
        Reads the non analytical directions, frequencies and eigendisplacements from the anaddb.nc file specified and
        adds them to the object.
        """
        self.non_anal_ph = NonAnalyticalPh.from_file(filepath)

    def __init__(self, structure, qpoints, phfreqs, phdispl_cart, non_anal_ph=None, amu=None):
        """
        Args:
            structure: :class:`Structure` object.
            qpoints: :class:`KpointList` instance.
            phfreqs: Phonon frequencies in eV.
            phdispl_cart: [nqpt, 3*natom, 3*natom] array with displacement in Cartesian coordinates in Angstrom.
                The last dimension stores the cartesian components.
            non_anal_ph: :class:`NonAnalyticalPh` with information of the non analytical contribution
                None if contribution is not present
            amu: dictionary that associates the atomic species present in the structure to the values of the atomic
                mass units used for the calculation
        """
        self.structure = structure

        # `KpointList` with the q-points
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

    def __str__(self):
        return self.to_string()

    def to_string(self, title=None, with_structure=True, with_qpoints=False, verbose=0):
        """
        Human-readable string with useful info such as band gaps, position of HOMO, LOMO...

        Args:
            with_structure: False if structural info shoud not be displayed.
            with_qpoints: False if q-point info shoud not be displayed.
            verbose: Verbosity level.
        """
        tos = str if verbose else repr
        lines = []; app = lines.append
        if title is not None:
            app(marquee(title, mark="="))

        if with_structure:
            app(tos(self.structure))
            app("")

        app("Number of q-points: %d" % self.num_qpoints)
        app("Atomic mass units: %s" % str(self.amu))
        has_dipdip = self.non_anal_ph is not None
        app("Has non-analytical contribution for q --> 0: %s" % has_dipdip)
        if verbose and has_dipdip:
            app(str(self.non_anal_ph))

        if with_qpoints:
            app(marquee("Q-points", mark="="))
            app(tos(self.qpoints))
            app("")

        return "\n".join(lines)

    #def displ_of_specie(self, specie):
    #    """Returns the displacement vectors for the given specie."""
    #    # TODO recheck the ordering
    #    # (nqpt, 3*natom, natom, 2) the last dimension stores the cartesian components.
    #    #raise NotImplementedError("")
    #    displ_specie = []
    #    for i, site in enumerate(self.structure):
    #        if site.specie == specie:
    #            displ_specie.append(self.phdispl_cart[:, :, i, :])

    #    return displ_specie

    @lazy_property
    def _auto_qlabels(self):
        # Find the q-point names in the pymatgen database.
        # We'll use _auto_klabels to label the point in the matplotlib plot
        # if qlabels are not specified by the user.
        _auto_qlabels = OrderedDict()
        for idx, qpoint in enumerate(self.qpoints):
            name = self.structure.findname_in_hsym_stars(qpoint)
            if name is not None:
                _auto_qlabels[idx] = name

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
        """Maximum phonon frequency."""
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

    @lazy_property
    def dyn_mat_eigenvect(self):
        """Eigenvalues of the dynamical matrix."""
        return get_dyn_mat_eigenvec(self.phdispl_cart, self.structure, self.amu)

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
        Write xmgrace file with phonon band structure energies and labels for high-symmetry q-points.

        Args:
            filepath: Filename
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
        """
        f = open(filepath, "wt")
        def w(s):
            f.write(s)
            f.write("\n")

        factor = _factor_ev2units(units)
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
        w('@yaxis  label "Phonon %s"' % _unit_tag(units))
        w('@yaxis  label char size 1.500000')
        w('@yaxis  ticklabel char size 1.500000')
        for nu in self.branches:
            w('@    s%d line color %d' % (nu, 1))

        # TODO: LO-TO splitting
        for nu in self.branches:
            w('@target G0.S%d' % nu)
            w('@type xy')
            for iq in range(self.num_qpoints):
                w('%d %.8E' % (iq, wqnu_units[iq, nu]))
            w('&')

        f.close()

    def qindex(self, qpoint):
        """
	Returns the index of the qpoint. Accepts integer or reduced coordinates.
	"""
        if duck.is_intlike(qpoint):
            return int(qpoint)
        else:
            return self.qpoints.index(qpoint)

    def qindex_qpoint(self, qpoint):
        """Returns (qindex, qpoint) from an integer or a qpoint."""
        qindex = self.qindex(qpoint)
        qpoint = self.qpoints(qindex)
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
                if freq < below_mev * 1000:
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

    def get_dict4frame(self, with_spglib=True):
        """
        Return a :class:`OrderedDict` with the most important parameters:

            - Chemical formula and number of atoms.
            - Lattice lengths, angles and volume.
            - The spacegroup number computed by Abinit (set to None if not available).
            - The spacegroup number and symbol computed by spglib (set to None not `with_spglib`).

        Useful to construct pandas DataFrames

        Args:
            with_spglib: If True, spglib is invoked to get the spacegroup symbol and number
        """
        odict = OrderedDict([
            ("nqpt", self.num_qpoints), ("nmodes", self.num_branches),
            ("min_freq", self.minfreq), ("max_freq", self.maxfreq),
            ("mean_freq", self.phfreqs.mean()), ("std_freq", self.phfreqs.std())

        ])
        odict.update(self.structure.get_dict4frame(with_spglib=with_spglib))

        return odict

    def get_phdos(self, method="gaussian", step=1.e-4, width=4.e-4):
        """
        Compute the phonon DOS on a linear mesh.

        Args:
            method: String defining the method
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Returns:
            :class:`PhononDos` object.

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
            raise ValueError("Method %s is not supported" % method)

        return PhononDos(mesh, values)

    def create_xyz_vib(self, iqpt, filename, pre_factor=200, do_real=True, scale_matrix=None, max_supercell=None):
        """
        Create vibration XYZ file for visualization of phonons.

        Args:
            iqpt: index of qpoint in self
            filename: name of the XYZ file that will be created
            pre_factor: Multiplication factor of the eigendisplacements
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
        This format can be read with V_sim or ascii-phonons.

        Args:
            iqpts: an index or a list of indices of the qpoints in self. Note that at present only V_sim supports
                an ascii file with multiple qpoints.
            filename: name of the ascii file that will be created.
            pre_factor: Multiplication factor of the eigendisplacements.
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

    def create_phononwebsite_json(self, filename, name=None, repetitions=None, highsym_qpts=None, match_bands=True,
                                  highsym_qpts_mode="split"):
        """
        Writes a json file that can be parsed from phononwebsite https://github.com/henriquemiranda/phononwebsite

        Args:
            filename: name of the json file that will be created
            name: name associated with the data.
            repetitions: number of repetitions of the cell. List of three integers. Defaults to [3,3,3].
            highsym_qpts: list of tuples. The first element of each tuple should be a list with the coordinates
                of a high symmetry point, the second element of the tuple should be its label.
            match_bands: if True tries to follow the band along the path based on the scalar product of the eigenvectors.
            highsym_qpts_mode: if highsym_qpts is None high symmetry q-points can be automatically determined.
                Accepts the following values:
                    'split' will split the path based on points where the path changes direction in the Brillouin zone.
                        Similar to the what is done in phononwebsite. Only Gamma will be labeled.
                    'std' uses the standard generation procedure for points and labels used in PhononBands.
                    None does not set any point.
        """

        def split_non_collinear(qpts):
            r"""
            function that splits the list of qpoints at repetitions (only the first point will be considered as
            high symm) and where the direction changes. Also sets $\Gamma$ for [0,0,0].
            Similar to what is done in phononwebsite.
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

        data = dict()
        data["name"] = name or self.structure.composition.reduced_formula
        data["natoms"] = self.num_atoms
        data["lattice"] = self.structure.lattice.matrix.tolist()
        data["atom_types"] = [e.name for e in self.structure.species]
        data["atom_numbers"] = self.structure.atomic_numbers
        data["chemical_symbols"] = self.structure.symbol_set
        data["atomic_numbers"] = list(set(self.structure.atomic_numbers))
        data["formula"] = self.structure.formula.replace(" ", "")
        data["repetitions"] = repetitions or (3,3,3)
        data["atom_pos_car"] = self.structure.cart_coords.tolist()
        data["atom_pos_red"] = self.structure.frac_coords.tolist()

        qpoints = []
        for q_sublist in self.split_qpoints:
            qpoints.extend(q_sublist.tolist())

        if highsym_qpts is None:
            if highsym_qpts_mode is None:
                data["highsym_qpts"] = []
            elif highsym_qpts_mode == 'split':
                data["highsym_qpts"] = split_non_collinear(qpoints)
            elif highsym_qpts_mode == 'std':
                data["highsym_qpts"] = list(six.zip(self._make_ticks_and_labels(None)))
        else:
            data["highsym_qpts"] = highsym_qpts

        distances = [0]
        for i in range(1,len(qpoints)):
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
            # the eigenvectors are needed (the mass factor is removed)
            vect = get_dyn_mat_eigenvec(phdispl_sublist, self.structure, amu=self.amu)
            # since phononwebsite will multiply again by exp(2*pi*q.r) this factor should be removed,
            # because in abinit convention it is included in the eigenvectors.
            for iqpt in range(len(qpts)):
                q = self.qpoints[iqpt].frac_coords
                for ai in range(self.num_atoms):
                    vect[iqpt, :, 3*ai:3*(ai + 1)] = vect[iqpt, :, 3*ai:3*(ai + 1)] * \
                        np.exp(-2*np.pi*1j*np.dot(self.structure[ai].frac_coords, q))

            if match_bands:
                vect = vect[np.arange(vect.shape[0])[:, None, None],
                            self.split_matched_indices[i][...,None],
                            np.arange(vect.shape[2])[None, None,:]]
            v = vect.reshape((len(vect), self.num_branches,self.num_atoms, 3))
            v = np.stack([v.real, v.imag], axis=-1)

            vectors.extend(v.tolist())

        data["qpoints"] = qpoints
        data["distances"] = distances
        data["eigenvalues"] = eigenvalues
        data["vectors"] = vectors

        with open(filename, 'wt') as json_file:
            json.dump(data, json_file, indent=2)

    def decorate_ax(self, ax, units='eV', **kwargs):
        title = kwargs.pop("title", None)
        if title is not None: ax.set_title(title)
        ax.grid(True)

        # Handle conversion factor.
        # TODO: Encapsulate this part.
        units = units.lower()
        if units == 'ev':
            ax.set_ylabel('Energy [eV]')
        elif units == 'mev':
            ax.set_ylabel('Energy [meV]')
        elif units == 'ha':
            ax.set_ylabel('Energy [Ha]')
        elif units in ('cm-1', 'cm^-1'):
            ax.set_ylabel(r'Frequency [cm$^{-1}$]')
        elif units == 'thz':
            ax.set_ylabel(r'Frequency [Thz]')
        else:
            raise ValueError('Value for units `{}` unknown'.format(units))

        # Set ticks and labels.
        ticks, labels = self._make_ticks_and_labels(kwargs.pop("qlabels", None))
        if ticks:
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False)

    @add_fig_kwargs
    def plot(self, ax=None, units="eV", qlabels=None, branch_range=None, match_bands=False, **kwargs):
        r"""
        Plot the phonon band structure.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.
            branch_range: Tuple specifying the minimum and maximum branch index to plot (default: all branches are plotted).
            match_bands: if True the bands will be matched based on the scalar product between the eigenvectors.

        Returns:
            `matplotlib` figure.
        """
        # Select the band range.
        if branch_range is None:
            branch_range = range(self.num_branches)
        else:
            branch_range = range(branch_range[0], branch_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, units=units, qlabels=qlabels)

        if "color" not in kwargs and not match_bands:
            kwargs["color"] = "black"

        if "linewidth" not in kwargs:
            kwargs["linewidth"] = 2.0

        # Plot the phonon branches.
        self.plot_ax(ax, branch_range, units=units, match_bands=match_bands, **kwargs)

        return fig

    def plot_ax(self, ax, branch, units='eV', match_bands=False, **kwargs):
        """
        Plots the frequencies for the given branches indices as a function of the q index on axis ax.
        If branch is None, all phonon branches are plotted.

        Return:
            The list of 'matplotlib' lines added.
        """
        if branch is None:
            branch_range = range(self.num_branches)
        elif isinstance(branch, (list, tuple, np.ndarray)):
            branch_range = branch
        else:
            branch_range = [branch]

        first_xx = 0
        lines = []

        factor = _factor_ev2units(units)

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
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}
            branch_range: Tuple specifying the minimum and maximum branch_i index to plot (default: all branches are plotted).
            colormap: matplotlib colormap to determine the colors available. The colors will be chosen not in a
                sequential order to avoid difficulties in distinguishing the lines.
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            max_colors: maximum number of colors to be used. If max_colors < num_braches the colors will be reapeated.
                It useful to better distinguish close bands when the number of branch is large.

        Returns:
            `matplotlib` figure.
        """
        # Select the band range.
        if branch_range is None:
            branch_range = range(self.num_branches)
        else:
            branch_range = range(branch_range[0], branch_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, units=units, qlabels=qlabels)

        first_xx = 0
        lines = []
        factor = _factor_ev2units(units)

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
            #     eigenvectors = get_dyn_mat_eigenvec(displ, self.structure, self.amu)
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
                eigenvectors = get_dyn_mat_eigenvec(displ, self.structure, self.amu)
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

    def _get_non_anal_freqs(self, direction):
        # directions for the qph2l in anaddb are given in cartesian coordinates
        direction = self.structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(direction)
        direction = direction / np.linalg.norm(direction)

        for i, d in enumerate(self.non_anal_directions):
            d = d / np.linalg.norm(d)
            if np.allclose(direction, d):
                return self.non_anal_phfreqs[i]

        raise ValueError("Non analytical contribution has not been calcolated for direction {0} ".format(direction))

    def _get_non_anal_phdispl(self, direction):
        # directions for the qph2l in anaddb are given in cartesian coordinates
        direction = self.structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(direction)
        direction = direction / np.linalg.norm(direction)

        for i, d in enumerate(self.non_anal_directions):
            d = d / np.linalg.norm(d)
            if np.allclose(direction, d):
                return self.non_anal_phdispl_cart[i]

        raise ValueError("Non analytical contribution has not been calcolated for direction {0} ".format(direction))

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

    @add_fig_kwargs
    def plot_fatbands(self, units="eV", colormap="jet", phdos_file=None,
                      alpha=0.7, max_stripe_width_mev=3.0, width_ratios=(2, 1),
                      qlabels=None, ylims=None,
                      **kwargs):
                      #cart_dir=None
        r"""
        Plot phonon fatbands and, optionally, atom-projected DOS.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            colormap: Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            phdos_file: Used to activate fatbands + PJDOS plot.
                Accept string with path of PHDOS.nc file or :class:`PhdosFile` object.
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            max_stripe_width_mev: The maximum width of the stripe in meV. Will be rescaled according to `units`.
            width_ratios: Ratios between the width of the fatbands plots and the DOS plots.
                Used if `phdos_file` is not None
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels. e.g. qlabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.

        Returns:
            `matplotlib` figure.
        """
        lw = kwargs.pop("lw", 2)
        factor = _factor_ev2units(units)
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
        from matplotlib.gridspec import GridSpec #, GridSpecFromSubplotSpec

        fig = plt.figure()
        nrows, ncols = (ntypat, 1) if phdos_file is None else (ntypat, 2)
        gspec = GridSpec(nrows=nrows, ncols=ncols, width_ratios=width_ratios if ncols == 2 else None,
                         wspace=0.05, hspace=0.1)

        cmap = plt.get_cmap(colormap)
        qq = list(range(self.num_qpoints))

        # phonon_displacements are in cartesian coordinates and stored in an array with shape
        # (nqpt, 3*natom, 3*natom) where the last dimension stores the cartesian components.
        # Precompute normalization factor:
        #   d2[q, nu] = \sum_{i=0}^{3*Nat-1) |d^{q\nu}_i|**2

        # FIXME there's a bug in anaddb since we should orthogonalize
        # wrt the phonon displacement as done (correctly) here
        d2_qnu = np.empty((self.num_qpoints, self.num_branches))
        for iq in range(self.num_qpoints):
            for nu in self.branches:
                cvect = self.phdispl_cart[iq, nu, :]
                d2_qnu[iq, nu] = np.vdot(cvect, cvect).real

        # Plot fatbands: one plot per atom type.
        ax00 = None
        for ax_row, symbol in enumerate(self.structure.symbol_set):
            ax = plt.subplot(gspec[ax_row, 0], sharex=ax00, sharey=ax00)
            if ax_row == 0: ax00 = ax
            self.decorate_ax(ax, units=units, qlabels=qlabels)
            color = cmap(float(ax_row) / (ntypat - 1))

            # dir_indices lists the coordinate indices for the atoms of the same type.
            atom_indices = self.structure.indices_from_symbol(symbol)
            dir_indices = []
            for aindx in atom_indices:
                start = 3 * aindx
                dir_indices.extend([start, start + 1, start + 2])
            dir_indices = np.array(dir_indices)

            for nu in self.branches:
                yy_qq = self.phfreqs[:, nu] * factor

                # Exctract the sub-vector associated to this atom type.
                displ_type = self.phdispl_cart[:, nu, dir_indices]
                d2_type = np.empty(self.num_qpoints)
                for iq in range(self.num_qpoints):
                    d2_type[iq] = np.vdot(displ_type[iq], displ_type[iq]).real

                # Normalize and scale by max_stripe_width_mev taking into account units.
                # The stripe is centered on the phonon branch hence the factor 2
                d2_type = factor * max_stripe_width_mev * 1.e-3 * d2_type / (2. * d2_qnu[:, nu])

                # Plot the phonon branch and the stripe.
                if nu == 0:
                    ax.plot(qq, yy_qq, lw=lw, label=symbol, color=color)
                else:
                    ax.plot(qq, yy_qq, lw=lw, color=color)

                ax.fill_between(qq, yy_qq + d2_type, yy_qq - d2_type, facecolor=color, alpha=alpha, linewidth=0)

            set_axlims(ax, ylims, "y")
            ax.legend(loc="best")

        # Type projected DOSes.
        if phdos_file is not None:
            ax01 = None
            for ax_row, symbol in enumerate(self.structure.symbol_set):
                color = cmap(float(ax_row) / (ntypat - 1))
                ax = plt.subplot(gspec[ax_row, 1], sharex=ax01, sharey=ax00)
                if ax_row == 0: ax01 = ax

                # Get PJDOS
                # Dictionary symbol --> partial PhononDos
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
    def plot_with_phdos(self, phdos, units="eV", qlabels=None, axlist=None, **kwargs):
        r"""
        Plot the phonon band structure with the phonon DOS.

        Args:
            phdos: An instance of :class:`PhononDos` or netcdf file providing a PhononDos object.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels e.g. qlabels = {(0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L"}.
            axlist: The axes for the bandstructure plot and the DOS plot. If axlist is None, a new figure
                is created and the two axes are automatically generated.

        Returns:
            `matplotlib` figure.
        """
        phdos = PhononDos.as_phdos(phdos, phdos_kwargs=None)

        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        if axlist is None:
            # Build axes and align bands and DOS.
            fig = plt.figure()
            gspec = GridSpec(1, 2, width_ratios=[2, 1])
            gspec.update(wspace=0.05)
            ax1 = plt.subplot(gspec[0])
            ax2 = plt.subplot(gspec[1], sharey=ax1)
        else:
            # Take them from axlist.
            ax1, ax2 = axlist
            fig = plt.gcf()

        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the phonon band structure.
        self.plot_ax(ax1, branch=None, units=units, **kwargs)
        self.decorate_ax(ax1, units=units, qlabels=qlabels)

        factor = _factor_ev2units(units)
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

    def to_dataframe(self):
        """
        Return a pandas DataFrame with the following columns:

            ['qidx', 'mode', 'freq', 'qpoint']

        where:

        ==============  ==========================
        Column          Meaning
        ==============  ==========================
        qidx            q-point index.
        mode            phonon branch index.
        freq            Phonon frequency in eV.
        qpoint          :class:`Kpoint` object
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
        Use seaborn to draw a box plot to show distributions of eigenvalues with respect to the mode index.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            mode_range: Only modes such as `mode_range[0] <= mode_index < mode_range[1]` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            kwargs: Keyword arguments passed to seaborn boxplot.
        """
        # Get the dataframe and select bands
        frame = self.to_dataframe()
        if mode_range is not None:
            frame = frame[(frame["mode"] >= mode_range[0]) & (frame["mode"] < mode_range[1])]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        factor = _factor_ev2units(units)
        yname = "freq %s" % _unit_tag(units)
        frame[yname] = factor * frame["freq"]

        import seaborn.apionly as sns
        hue = None
        ax = sns.boxplot(x="mode", y=yname, data=frame, hue=hue, ax=ax, **kwargs)
        if swarm:
            sns.swarmplot(x="mode", y=yname, data=frame, hue=hue, color=".25", ax=ax)

        return fig

    def to_pymatgen(self, qlabels=None):
        r"""
        Creates a pymatgen PhononBandStructureSymmLine object.

        Args:
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the q-points.
                The values are the labels e.g. qlabels = {(0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L"}.
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

        ph_freqs = np.transpose(ph_freqs)
        qpts = np.array(qpts)
        displ = np.transpose(displ, (1, 0, 2, 3))

        return PhononBandStructureSymmLine(qpoints=qpts, frequencies=ph_freqs,
                                           lattice=self.structure.reciprocal_lattice,
                                           has_nac=self.non_anal_ph is not None, eigendisplacements=displ,
                                           labels_dict=labels_dict, structure=self.structure)


class PHBST_Reader(ETSF_Reader):
    """This object reads data from PHBST.nc file produced by anaddb."""

    def read_qredcoords(self):
        """Array with the reduced coordinates of the q-points."""
        return self.read_value("qpoints")

    def read_qweights(self):
        """The weights of the q-points"""
        return self.read_value("qweights")

    def read_phfreqs(self):
        """Array with the phonon frequencies in eV."""
        return self.read_value("phfreqs")

    def read_phdispl_cart(self):
        """
        Complex array with the Cartesian displacements in **Angstrom**
        shape is (num_qpoints,  mu_mode,  cart_direction).
        """
        return self.read_value("phdispl_cart", cmode="c")

    def read_amu(self):
        """The atomic mass units"""
        return self.read_value("atomic_mass_units", default=None)


class PhbstFile(AbinitNcFile, Has_Structure, Has_PhononBands, NotebookWriter):

    def __init__(self, filepath):
        """
        Object used to access data stored in the PHBST file produced by ABINIT.

        Args:
            path: path to the file
        """
        super(PhbstFile, self).__init__(filepath)
        self.reader = PHBST_Reader(filepath)

        # Initialize Phonon bands
        self._phbands = PhononBands.from_file(filepath)

    @property
    def structure(self):
        """:class:`Structure` object"""
        return self.phbands.structure

    @property
    def qpoints(self):
        """List of q-point objects."""
        return self.phbands.qpoints

    @property
    def phbands(self):
        """:class:`PhononBands` object"""
        return self._phbands

    def close(self):
        """Close the file."""
        self.reader.close()

    def qindex(self, qpoint):
        """
        Returns the index of the qpoint. Accepts integer or reduced coordinates.
        """
        if duck.is_intlike(qpoint):
            return int(qpoint)
        else:
            return self.qpoints.index(qpoint)

    def qindex_qpoint(self, qpoint):
        """
        Returns (qindex, qpoint) from an integer or a qpoint.
        """
        qindex = self.qindex(qpoint)
        return qindex, self.qpoints[qindex]

    def get_phframe(self, qpoint, with_structure=True):
        """
        Return a pandas :class:`DataFrame` with the phonon frequencies at the given q-point and
        information on the crystal structure (used for convergence studies).

        Args:
            qpoint: integer, vector of reduced coordinates or :class:`Kpoint` object.
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
            d.update(self.structure.get_dict4frame(with_spglib=True))

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

    def write_notebook(self, nbpath=None):
        """
        Write an jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("fig = ncfile.phbands.plot()"),
            nbv.new_code_cell("fig = ncfile.phbands.qpoints.plot()"),
            #nbv.new_code_cell("fig = ncfile.phbands.get_phdos().plot()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


_THERMO_YLABELS = {  # [name][units] --> latex string
    "internal_energy": {"eV": "$U(T)$ [eV/cell]", "Jmol": "$U(T)$ [J/mole]"},
    "free_energy": {"eV": "$F(T) + ZPE$ [eV/cell]", "Jmol": "$F(T) + ZPE$ [J/mole]"},
    "entropy": {"eV": "$S(T)$ [eV/cell]", "Jmol": "$S(T)$ [J/mole]"},
    "cv": {"eV": "$C_V(T)$ [eV/cell]", "Jmol": "$C_V(T)$ [J/mole]"},
}

class PhononDos(Function1D):
    """
    This object stores the phonon density of states.
    An instance of `PhononDos` has a `mesh` (numpy array with the points of the mesh)
    and another numpy array, `values`, with the DOS on the mesh..

    .. note::

        mesh is given in eV, values are in states/eV.
    """
    #def __init__(self, mesh, values, qmesh):
    #    super(PhononDos, self).__init__(mesh, values)
    #    self.qmesh = qmesh

    @classmethod
    def as_phdos(cls, obj, phdos_kwargs=None):
        """
        Return an instance of :class:`PhononDOS` from a generic obj.
        Supports:

            - instances of cls
            - files (string) that can be open with abiopen and that provide one of the following attributes:
                [`phdos`, `phbands`]
            - instances of `PhononBands`
            - objects providing a `phbands`attribute.

        Args:
            phdos_kwargs: optional dictionary with the options passed to `get_phdos` to compute the phonon DOS.
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
        index of the first point in the mesh whose value is >= 0
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
        Helper function to plot DOS/IDOS on the axis ax.

        Args:
            ax: matplotlib axis
            what: string selecting the quantity to plot:
                "d" for DOS, "i" for IDOS. chars can be concatenated
                hence what="id" plots both IDOS and DOS. (default "d").
            exchange_xy: True to exchange axis
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            kwargs: Options passed to matplotlib plot method.

        Return:
            list of lines added to the plot
        """
        opts = [c.lower() for c in what]
        lines = []

        for c in opts:
            f = {"d": self, "i": self.idos}[c]
            xfactor = _factor_ev2units(units)
            # Don't rescale IDOS
            yfactor = 1 / xfactor if c == "d" else 1

            ls = f.plot_ax(ax, exchange_xy=exchange_xy, xfactor=xfactor, yfactor=yfactor, **kwargs)
            lines.extend(ls)

        return lines

    @add_fig_kwargs
    def plot(self, units="eV", **kwargs):
        """
        Plot Phonon DOS and IDOS on two distict plots.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            kwargs: Keyword arguments passed to :mod:`matplotlib`.

        Returns:
            `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        fig = plt.figure()
        gspec = GridSpec(2, 1, height_ratios=[1, 2], wspace=0.05)
        ax1 = plt.subplot(gspec[0])
        ax2 = plt.subplot(gspec[1])

        for ax in (ax1, ax2):
            ax.grid(True)

        ax2.set_xlabel('Energy %s' % _unit_tag(units))
        ax1.set_ylabel("IDOS [states]")
        ax2.set_ylabel("DOS %s" % _dos_label_from_units(units))

        self.plot_dos_idos(ax1, what="i", units=units, **kwargs)
        self.plot_dos_idos(ax2, what="d", units=units, **kwargs)

        return fig

    def get_internal_energy(self, tstart=5, tstop=300, num=50):
        """
        Returns the internal energy, in eV, in the harmonic approximation for different temperatures
        Zero point energy is included.

        tstart: The starting value (in Kelvin) of the temperature mesh.
        tstop: The end value (in Kelvin) of the mesh.
        num: int, optional Number of samples to generate. Default is 50.

        Return:
            :class:`Function1D` with U(T) + ZPE
        """
        tmesh = np.linspace(tstart, tstop, num=num)
        w, gw = self.mesh[self.iw0:], self.values[self.iw0:]
        coth = lambda x: 1.0 / np.tanh(x)

        vals = np.empty(len(tmesh))
        for it, temp in enumerate(tmesh):
            wd2kt = w / (2 * abu.kb_eVK * temp)
            vals[it] = np.trapz(w * coth(wd2kt) * gw, x=w)

        return Function1D(tmesh, 0.5 * vals + self.zero_point_energy)

    def get_entropy(self, tstart=5, tstop=300, num=50):
        """
        Returns the entropy, in eV/K, in the harmonic approximation for different temperatures

        tstart: The starting value (in Kelvin) of the temperature mesh.
        tstop: The end value (in Kelvin) of the mesh.
        num: int, optional Number of samples to generate. Default is 50.

        Return:
            :class:`Function1D` with S(T).
        """
        tmesh = np.linspace(tstart, tstop, num=num)
        w, gw = self.mesh[self.iw0:], self.values[self.iw0:]
        coth = lambda x: 1.0 / np.tanh(x)

        vals = np.empty(len(tmesh))
        for it, temp in enumerate(tmesh):
            wd2kt = w / (2 * abu.kb_eVK * temp)
            vals[it] = np.trapz((wd2kt * coth(wd2kt) - np.log(2 * np.sinh(wd2kt))) * gw, x=w)

        return Function1D(tmesh, abu.kb_eVK * vals)

    def get_free_energy(self, tstart=5, tstop=300, num=50):
        """
        Returns the free energy, in eV, in the harmonic approximation for different temperatures
        Zero point energy is included.

        tstart: The starting value (in Kelvin) of the temperature mesh.
        tstop: The end value (in Kelvin) of the mesh.
        num: int, optional Number of samples to generate. Default is 50.

        Return:
            :class:`Function1D` with F(T) = U(T) + ZPE - T x S(T)
        """
        uz = self.get_internal_energy(tstart=tstart, tstop=tstop, num=num)
        s = self.get_entropy(tstart=tstart, tstop=tstop, num=num)

        return Function1D(uz.mesh, uz.values - s.mesh * s.values)

    def get_cv(self, tstart=5, tstop=300, num=50):
        """
        Returns the constant-volume specific heat, in eV/K, in the harmonic approximation for different temperatures

        tstart: The starting value (in Kelvin) of the temperature mesh.
        tstop: The end value (in Kelvin) of the mesh.
        num: int, optional Number of samples to generate. Default is 50.

        Return:
            :class:`Function1D` with C_v(T).
        """
        tmesh = np.linspace(tstart, tstop, num=num)
        w, gw = self.mesh[self.iw0:], self.values[self.iw0:]
        csch2 = lambda x: 1.0 / (np.sinh(x) ** 2)

        vals = np.empty(len(tmesh))
        for it, temp in enumerate(tmesh):
            wd2kt = w / (2 * abu.kb_eVK * temp)
            vals[it] = np.trapz(wd2kt ** 2 * csch2(wd2kt) * gw, x=w)

        return Function1D(tmesh, abu.kb_eVK * vals)

    @add_fig_kwargs
    def plot_harmonic_thermo(self, tstart=5, tstop=300, num=50, units="eV", formula_units=None,
                             quantities=None, **kwargs):
        """
        Plot thermodinamic properties from the phonon DOSes within the harmonic approximation.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 50.
            quantities: List of strings specifying the thermodinamic quantities to plot.
                Possible values: ["internal_energy", "free_energy", "entropy", "c_v"].
                None means all.
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            formula_units:
                the number of formula units per unit cell. If unspecified, the
                thermodynamic quantities will be given on a per-unit-cell basis.

        Returns:
            matplotlib figure.
        """
        quantities = list_strings(quantities) if quantities is not None else \
            ["internal_energy", "free_energy", "entropy", "cv"]

        # Build grid of plots.
        ncols, nrows = 1, 1
        num_plots = len(quantities)
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        import matplotlib.pyplot as plt
        fig, axmat = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=False, squeeze=False)
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: axmat[-1, -1].axis("off")

        for iax, (qname, ax) in enumerate(zip(quantities, axmat.flat)):
            # Compute thermodinamic quantity associated to qname.
            f1d = getattr(self, "get_" + qname)(tstart=tstart, tstop=tstop, num=num)
            ys = f1d.values
            if formula_units is not None: ys /= formula_units
            if units == "Jmol": ys = ys * abu.e_Cb * abu.Avogadro
            ax.plot(f1d.mesh, ys)

            ax.set_title(qname)
            ax.grid(True)
            ax.set_xlabel("Temperature [K]")
            ax.set_ylabel(_THERMO_YLABELS[qname][units])
            #ax.legend(loc="best")

        fig.tight_layout()
        return fig


class PhdosReader(ETSF_Reader):
    """
    This object reads data from the PHDOS.nc file produced by anaddb.

    .. note::

            Frequencies are in eV, DOSes are in states/eV.
    """
    @lazy_property
    def structure(self):
        """The crystalline structure."""
        return self.read_structure()

    @lazy_property
    def wmesh(self):
        """The frequency mesh for the PH-DOS in eV."""
        return self.read_value("wmesh")

    def read_pjdos_type(self):
        """[ntypat, nomega] array with PH-DOS projected over atom types."""
        return self.read_value("pjdos_type")

    def read_pjdos_atdir(self):
        """
        Return [natom, three, nomega] array with Phonon DOS projected over atoms and reduced directions.
        """
        return self.read_value("pjdos")

    def read_phdos(self):
        """Return the :class:`PhononDOS`. with the total phonon DOS"""
        return PhononDos(self.wmesh, self.read_value("phdos"))

    def read_pjdos_symbol_rc_dict(self):
        """
        Return `OrderedDict` mapping element symbol --> [3, nomega] array
        with the the phonon DOSes summed over atom-types and decomposed along
        the three reduced directions.
        """
        # phdos_rc_type[ntypat, 3, nomega]
        values = self.read_value("pjdos_rc_type")

        od = OrderedDict()
        for symbol in self.chemical_symbols:
           type_idx = self.typeidx_from_symbol(symbol)
           od[symbol] = values[type_idx]

        return od

    def read_pjdos_symbol_dict(self):
        """
        Ordered dictionary mapping element symbol --> `PhononDos`
        where PhononDos is the contribution to the total DOS summed over atoms
        with chemical symbol `symbol`.
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
    PHDOS.nc file produced by anaddb. Provides helper function
    to visualize/extract data.
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
    def structure(self):
        """Returns the :class:`Structure` object."""
        return self.reader.structure

    @lazy_property
    def phdos(self):
        """:class:`PhononDos` object."""
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
    def plot_pjdos_type(self, units="eV", stacked=True, colormap="jet", alpha=0.7,
                        ax=None, xlims=None, ylims=None, **kwargs):
        """
        Plot type-projected phonon DOS.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            colormap: Have a look at the colormaps
                `here <http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html>`_
                and decide which one you'd like:
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            ylims: y-axis limits.

        Returns:
            matplotlib figure.
        """
        lw = kwargs.pop("lw", 2)
        factor = _factor_ev2units(units)

        ax, fig, plt = get_ax_fig_plt(ax)
        cmap = plt.get_cmap(colormap)

        ax.grid(True)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        ax.set_xlabel('Frequency %s' % _unit_tag(units))
        ax.set_ylabel('PJDOS %s' % _dos_label_from_units(units))

        # Type projected DOSes.
        num_plots = len(self.pjdos_symbol)
        cumulative = np.zeros(len(self.wmesh))

        for i, (symbol, pjdos) in enumerate(self.pjdos_symbol.items()):
            x, y = pjdos.mesh * factor, pjdos.values / factor
            color = cmap(float(i) / (num_plots - 1))
            if not stacked:
                ax.plot(x, y, lw=lw, label=symbol, color=color)
            else:
                ax.plot(x, cumulative + y, lw=lw, label=symbol, color=color)
                ax.fill_between(x, cumulative, cumulative + y, facecolor=color, alpha=alpha)
                cumulative += y

        # Total PHDOS
        x, y = self.phdos.mesh * factor, self.phdos.values / factor
        ax.plot(x, y, lw=lw, label="Total PHDOS", color='black')
        ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_pjdos_redirs_type(self, units="eV", stacked=True, colormap="jet", alpha=0.7,
                               xlims=None, ylims=None, axlist=None, **kwargs):
        """
        Plot type-projected phonon DOS decomposed along the three reduced directions.
        Three rows for each reduced direction. Each row shows the contribution of each atomic type + Total PH DOS.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            colormap: Have a look at the colormaps
                `here <http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html>`_
                and decide which one you'd like:
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            ylims: y-axis limits.
            axlist: List of matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            matplotlib figure.
        """
        lw = kwargs.pop("lw", 2)
        ntypat = self.structure.ntypesp
        factor = _factor_ev2units(units)

        # Three rows for each reduced direction.
        # Each row shows the contribution of each atomic type + Total PH DOS.
        import matplotlib.pyplot as plt
        nrows, ncols = 3, 1
        if axlist is None:
            fig, axlist = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, squeeze=True)
        else:
            axlist = np.reshape(axlist, (nrows, ncols)).ravel()

        cmap = plt.get_cmap(colormap)

        # symbol --> [three, number_of_frequencies]
        pjdos_symbol_rc = self.reader.read_pjdos_symbol_rc_dict()

        xx = self.phdos.mesh * factor
        for idir, ax in enumerate(axlist):
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")

            if idir in (0, 2):
                ax.set_ylabel(r'PJDOS along $L_{%d}$' % idir)
                if idir == 2:
                    ax.set_xlabel('Frequency %s' % _unit_tag(units))

            # Plot Type projected DOSes along reduced direction idir
            cumulative = np.zeros(len(self.wmesh))
            for itype, symbol in enumerate(self.reader.chemical_symbols):
                color = cmap(float(itype) / (ntypat - 1))
                yy = pjdos_symbol_rc[symbol][idir] / factor

                if not stacked:
                    ax.plot(xx, yy, label=symbol, color=color)
                else:
                    ax.plot(xx, cumulative + yy, lw=lw, label=symbol, color=color)
                    ax.fill_between(xx, cumulative, cumulative + yy, facecolor=color, alpha=alpha)
                    cumulative += yy

            # Add Total PHDOS
            ax.plot(xx, self.phdos.values / factor, lw=lw, label="Total PHDOS", color='black')
            ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_pjdos_redirs_site(self, view="inequivalent", units="eV", stacked=True, colormap="jet", alpha=0.7,
                               xlims=None, ylims=None, axlist=None, **kwargs):
        """
        Plot phonon PJDOS for each atom in the unit cell. By default, only "inequivalent" atoms are shown.

        Args:
            view: "inequivalent", "all"
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            stacked: True if DOS partial contributions should be stacked on top of each other.
            colormap: matplotlib colormap.
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            axlist: List of matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        # Define num_plots and ax2atom depending on view.
        factor = _factor_ev2units(units)
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

        # Three rows for each reduced direction.
        # Each row shows the contribution of each site + Total PH DOS.
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(colormap)
        nrows, ncols = 3, 1
        if axlist is None:
            fig, axlist = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, squeeze=True)
        else:
            axlist = np.reshape(axlist, (nrows, ncols)).ravel()

        # [natom, three, nomega] array with PH-DOS projected over atoms and reduced directions"""
        pjdos_atdir = self.reader.read_pjdos_atdir()

        xx = self.phdos.mesh * factor
        for idir, ax in enumerate(axlist):
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")

            if idir in (0, 2):
                ax.set_ylabel(r'PJDOS along $L_{%d}$' % idir)
                if idir == 2:
                    ax.set_xlabel('Frequency %s' % _unit_tag(units))

            # Plot Type projected DOSes along reduced direction idir
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
            ax.legend(loc="best")

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("fig = ncfile.phdos.plot()"),
            nbv.new_code_cell("fig = ncfile.plot_pjdos_type()"),
            nbv.new_code_cell("fig = ncfile.plot_pjdos_redirs_type()"),
            #nbv.new_code_cell("fig = ncfile.plot_pjdos_redirs_site()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


# FIXME: Remove. Use PhononBandsPlotter API.
@add_fig_kwargs
def phbands_gridplot(phb_objects, titles=None, phdos_objects=None, phdos_kwargs=None, units="eV", **kwargs):
    """
    Plot multiple phonon bandstructures and optionally DOSes on a grid.

    Args:
        phb_objects: List of objects from which the phonon band structures are extracted.
            Each item in phb_objects is either a string with the path of the netcdf file,
            or one of the abipy object with an `phbands` attribute or a :class:`PhononBands` object.
        phdos_objects:
            List of objects from which the phonon DOSes are extracted.
            Accept filepaths or :class:`PhononDos` objects. If phdos_objects is not None,
            each subplot in the grid contains a band structure with DOS else a simple bandstructure plot.
        titles:
            List of strings with the titles to be added to the subplots.
        phdos_kwargs: optional dictionary with the options passed to `get_phdos` to compute the phonon DOS.
            Used only if `phdos_objects` is not None.
        units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.

    Returns:
        matplotlib figure.
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
        # Plot grid with bands only.
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, squeeze=False)
        axes = axes.ravel()
        # don't show the last ax if numeb is odd.
        if numeb % ncols != 0: axes[-1].axis("off")

        for i, (phbands, ax) in enumerate(zip(phbands_list, axes)):
            phbands.plot(ax=ax, units=units, show=False)
            if titles is not None: ax.set_title(titles[i])
            if i % ncols != 0:
                ax.set_ylabel("")

    else:
        # Plot grid with bands + DOS
        # see http://matplotlib.org/users/gridspec.html
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        fig = plt.figure()
        gspec = GridSpec(nrows, ncols)

        for i, (phbands, phdos) in enumerate(zip(phbands_list, phdos_list)):
            subgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=gspec[i], width_ratios=[2, 1], wspace=0.05)
            # Get axes and align bands and DOS.
            ax1 = plt.subplot(subgrid[0])
            ax2 = plt.subplot(subgrid[1], sharey=ax1)
            phbands.plot_with_phdos(phdos, axlist=(ax1, ax2), units=units, show=False)

            if titles is not None: ax1.set_title(titles[i])
            if i % ncols != 0:
                for ax in (ax1, ax2):
                    ax.set_ylabel("")

    return fig


def frame_from_phbands(phbands_objects, index=None, with_spglib=True):
    """
    Build a pandas dataframe with the most important results available in a list of band structures.

    Args:
        phbands_objects: List of objects that can be converted to structure.
            Support netcdf filenames or :class:`PhononBands` objects
            See `PhononBands.as_ebands` for the complete list.
        index: Index of the dataframe.
        with_spglib: If True, spglib is invoked to get the spacegroup symbol and number.

    Return:
        pandas :class:`DataFrame`
    """
    phbands_list = [PhononBands.as_phbands(obj) for obj in phbands_objects]
    # Use OrderedDict to have columns ordered nicely.
    odict_list = [(phbands.get_dict4frame(with_spglib=with_spglib)) for phbands in phbands_list]

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
    _LINE_COLORS = ["b", "r", "k", "g"]
    _LINE_STYLES = ["-", ":", "--", "-.",]
    _LINE_WIDTHS = [2, ]

    def __init__(self, key_phbands=None, key_phdos=None, phdos_kwargs=None):
        """
        Args:
            key_phbands: List of (label, phbands) tuples.
                phbands is any object that can be converted into :class:`PhononBands` e.g. ncfile, path.
            key_phdos: List of (label, phdos) tuples.
                phdos is any object that can be converted into :class:`PhononDos`
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

    def get_phbands_frame(self, with_spglib=True):
        """
        Build a pandas dataframe with the most important results available in the band structures.
        """
        return frame_from_phbands(list(self.phbands_dict.values()),
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
        """"List of `:class:PhononBands`."""
        return list(self._bands_dict.values())

    @property
    def phdoses_list(self):
        """"List of :class:`PhononDos`."""
        return list(self._phdoses_dict.values())

    def iter_lineopt(self):
        """Generates style options for lines."""
        for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    @deprecated(message="add_phbands_from_file method of PhononBandsPlotter has been replaced by add_phbands. It will be removed in 0.4")
    def add_phbands_from_file(self, filepath, label=None):
        """
        Adds a band structure for plotting. Reads data from a Netcdfile
        """
        from abipy.abilab import abiopen
        with abiopen(filepath) as ncfile:
            if label is None: label = ncfile.filepath
            self.add_phbands(label, ncfile.phbands)

    def add_phbands(self, label, bands, phdos=None, dos=None, phdos_kwargs=None):
        """
        Adds a band structure for plotting.

        Args:
            label: label for the bands. Must be unique.
            bands: :class:`PhononBands` object.
            phdos: :class:`PhononDos` object.
            phdos_kwargs: optional dictionary with the options passed to `get_phdos` to compute the phonon DOS.
              Used only if `phdos` is not None.
        """
        if dos is not None:
            import warnings
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
    def combiplot(self, qlabels=None, units='eV', ylims=None, **kwargs):
        r"""
        Plot the band structure and the DOS on the same figure.
        Use `gridplot` to plot band structures on different figures.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            qlabels: dictionary whose keys are tuples with the reduced coordinates of the k-points.
                The values are the labels e.g. klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used

        Returns:
            matplotlib figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        # Build grid of plots.
        fig = plt.figure()
        if self.phdoses_dict:
            gspec = GridSpec(1, 2, width_ratios=[2, 1])
            gspec.update(wspace=0.05)
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

        # Plot bands.
        lines, legends = [], []
        my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        for (label, bands), lineopt in zip(self._bands_dict.items(), self.iter_lineopt()):
            i += 1
            my_kwargs.update(lineopt)
            opts_label[label] = my_kwargs.copy()

            l = bands.plot_ax(ax1, branch=None, units=units, **my_kwargs)
            lines.append(l[0])

            # Use relative paths if label is a file.
            if os.path.isfile(label):
                legends.append("%s" % os.path.relpath(label))
            else:
                legends.append("%s" % label)

            # Set ticks and labels, legends.
            if i == 0:
                bands.decorate_ax(ax1, qlabels=qlabels, units=units)

        ax1.legend(lines, legends, loc='best', shadow=True)

        # Add DOSes
        if self.phdoses_dict:
            ax = ax_list[1]
            for label, dos in self.phdoses_dict.items():
                dos.plot_dos_idos(ax, exchange_xy=True, units=units, **opts_label[label])

        return fig

    @deprecated(message="plot method of PhononBandsPlotter has been replaced by combiplot. It will be removed in 0.4")
    def plot(self, *args, **kwargs):
        return self.combiplot(*args, **kwargs)

    @add_fig_kwargs
    def gridplot(self, with_dos=True, units="eV", **kwargs):
        """
        Plot multiple electron bandstructures and optionally DOSes on a grid.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            with_dos: True if DOS should be printed.

        Returns:
            matplotlib figure.
        """
        titles = list(self._bands_dict.keys())
        phb_objects = list(self._bands_dict.values())
        phdos_objects = None
        if self.phdoses_dict and with_dos:
            phdos_objects = list(self.phdoses_dict.values())

        return phbands_gridplot(phb_objects, titles=titles, phdos_objects=phdos_objects, units=units, show=False)

    @add_fig_kwargs
    def boxplot(self, mode_range=None, units="eV", swarm=False, **kwargs):
        """
        Use seaborn to draw a box plot to show distributions of eigenvalues with respect to the band index.
        Band structures are drawn on different subplots.

        Args:
            mode_range: Only bands such as `mode_range[0] <= nu_index < mode_range[1]` are included in the plot.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            swarm: True to show the datapoints on top of the boxes
            kwargs: Keywork arguments passed to seaborn boxplot.
        """
        # Build grid of plots.
        num_plots, ncols, nrows = len(self.phbands_dict), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, squeeze=False)
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
        Use seaborn to draw a box plot comparing the distributions of the frequencies.
        Phonon Band structures are drawn on the same plot.

        Args:
            mode_range: Only bands such as `mode_range[0] <= nu_index < mode_range[1]` are included in the plot.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            swarm: True to show the datapoints on top of the boxes
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            kwargs: Keyword arguments passed to seaborn boxplot.
        """
        frames = []
        for label, phbands in self.phbands_dict.items():
            # Get the dataframe, select bands and add column with label
            frame = phbands.to_dataframe()
            if mode_range is not None:
                frame = frame[(frame["mode"] >= mode_range[0]) & (frame["mode"] < mode_range[1])]
            frame["label"] = label
            frames.append(frame)

        # Merge frames ignoring index (not meaningful here)
        import pandas as pd
        data = pd.concat(frames, ignore_index=True)

        import seaborn.apionly as sns
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        # Create column with frequencies in `units`.
        factor = _factor_ev2units(units)
        yname = "freq %s" % _unit_tag(units)
        data[yname] = factor * data["freq"]

        sns.boxplot(x="mode", y=yname, data=data, hue="label", ax=ax, **kwargs)
        if swarm:
            sns.swarmplot(x="mode", y=yname, data=data, hue="label", color=".25", ax=ax)

        return fig

    def animate(self, interval=250, savefile=None, units="eV", width_ratios=(2, 1), show=True):
        """
        Use matplotlib to animate a list of band structure plots (with or without DOS).

        Args:
            interval: draws a new frame every interval milliseconds.
            savefile: Use e.g. 'myanimation.mp4' to save the animation in mp4 format.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            width_ratios: Defines the ratio between the band structure plot and the dos plot.
                Used when there are DOS stored in the plotter.
            show: True if the animation should be shown immediately

        Returns:
            Animation object.

        See also:
            http://matplotlib.org/api/animation_api.html
            http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

        Note:
            It would be nice to animate the title of the plot, unfortunately
            this feature is not available in the present version of matplotlib.
            See: http://stackoverflow.com/questions/17558096/animated-title-in-matplotlib
        """
        phbands_list, phdos_list = self.phbands_list, self.phdoses_list
        if phdos_list and len(phdos_list) != len(phbands_list):
            raise ValueError("The number of objects for DOS must be equal to the number of bands")

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
            gspec = GridSpec(1, 2, width_ratios=width_ratios)
            gspec.update(wspace=0.05)
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

    def ipw_select_plot(self):
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
            phdos: :class:`PhononDos` object.
            phdos_kwargs: optional dictionary with the options passed to `get_phdos` to compute the phonon DOS.
                Used when phdos is not already an instance of `cls` or when we have to compute the DOS from obj.
        """
        if label in self._phdoses_dict:
            raise ValueError("label %s is already in %s" % (label, list(self._phdoses_dict.keys())))

        self._phdoses_dict[label] = PhononDos.as_phdos(phdos, phdos_kwargs)

    @add_fig_kwargs
    def combiplot(self, ax=None, units="eV", xlims=None, ylims=None, **kwargs):
        """
        Plot DOSes on the same figure. Use `gridplot` to plot DOSes on different figures.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            ylims: y-axis limits.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        ax.grid(True)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        ax.set_xlabel('Energy %s' % _unit_tag(units))
        ax.set_ylabel('DOS %s' % _dos_label_from_units(units))

        lines, legends = [], []
        for label, dos in self._phdoses_dict.items():
            l = dos.plot_dos_idos(ax, units=units, **kwargs)[0]
            lines.append(l)
            legends.append("DOS: %s" % label)

        # Set legends.
        ax.legend(lines, legends, loc='best', shadow=True)

        return fig

    @deprecated(message="plot method of PhononDos has been replaced by combiplot. It will be removed in 0.4")
    def plot(self, **kwargs):
        return self.combiplot(**kwargs)

    @add_fig_kwargs
    def gridplot(self, units="eV", xlims=None, ylims=None, **kwargs):
        """
        Plot multiple DOSes on a grid.

        Args:
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used

        Returns:
            matplotlib figure.
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
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, squeeze=False)
        axes = axes.ravel()
        # don't show the last ax if numeb is odd.
        if numeb % ncols != 0: axes[-1].axis("off")

        for i, (label, phdos) in enumerate(self._phdoses_dict.items()):
            ax = axes[i]
            phdos.plot_dos_idos(ax, units=units)

            ax.set_xlabel('Energy %s' % _unit_tag(units))
            ax.set_ylabel("DOS %s" % _dos_label_from_units(units))
            ax.set_title(label)
            ax.grid(True)
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")
            if i % ncols != 0:
                ax.set_ylabel("")

        return fig

    @add_fig_kwargs
    def plot_harmonic_thermo(self, tstart=5, tstop=300, num=50, units="eV", formula_units=1,
                             quantities="all", **kwargs):
        """
        Plot thermodinamic properties from the phonon DOS within the harmonic approximation.

        Args:
            tstart: The starting value (in Kelvin) of the temperature mesh.
            tstop: The end value (in Kelvin) of the mesh.
            num: int, optional Number of samples to generate. Default is 50.
            units: eV for energies in ev/unit_cell, Jmol for results in J/mole.
            formula_units:
                the number of formula units per unit cell. If unspecified, the
                thermodynamic quantities will be given on a per-unit-cell basis.
            quantities: List of strings specifying the thermodinamic quantities to plot.
                Possible values: ["internal_energy", "free_energy", "entropy", "c_v"].

        Returns:
            matplotlib figure.
        """
        quantities = list_strings(quantities) if quantities != "all" else \
            ["internal_energy", "free_energy", "entropy", "cv"]

        # Build grid of plots.
        ncols, nrows = 1, 1
        num_plots = len(quantities)
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        import matplotlib.pyplot as plt
        fig, axmat = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=False, squeeze=False)
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: axmat[-1, -1].axis("off")

        for iax, (qname, ax) in enumerate(zip(quantities, axmat.flat)):

            for i, (label, phdos) in enumerate(self._phdoses_dict.items()):
                # Compute thermodinamic quantity associated to qname.
                f1d = getattr(phdos, "get_" + qname)(tstart=tstart, tstop=tstop, num=num)
                ys = f1d.values
                if formula_units != 1: ys /= formula_units
                if units == "Jmol": ys = ys * abu.e_Cb * abu.Avogadro
                ax.plot(f1d.mesh, ys, label=label)

            ax.set_title(qname)
            ax.grid(True)
            ax.set_xlabel("Temperature [K]")
            ax.set_ylabel(_THERMO_YLABELS[qname][units])
            ax.legend(loc="best")

        fig.tight_layout()
        return fig

    def ipw_select_plot(self):
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

    def ipw_harmonic_thermo(self):
        """
        Return an ipython widget with controllers to plot thermodinamic properties
        from the phonon DOS within the harmonic approximation.
        """
        def plot_callback(tstart, tstop, num, units, formula_units):
            self.plot_harmonic_thermo(tstart=tstart, tstop=tstop, num=num,
                                      units=units, formula_units=formula_units, show=True)

        import ipywidgets as ipw
        return ipw.interact_manual(
                plot_callback,
                tstart=5, tstop=300, num=50, units=["eV", "Jmol"], formula_units=1)

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
            structure: :class:`Structure` object.
            directions: Cartesian directions along which the non analytical frequencies have been calculated
            phfreqs: Phonon frequencies with non analytical contribution in eV along directions
            phdispl_cart: Displacement in Angstromg in Cartesian coordinates with non analytical contribution
                along directions
            amu: dictionary that associates the atomic species present in the structure to the values of the atomic
                mass units used for the calculation
        """
        self._structure = structure
        self.directions = directions
        self.phfreqs = phfreqs
        self.phdispl_cart = phdispl_cart
        self.amu = amu

    @classmethod
    def from_file(cls, filepath):
        """
        Reads the non analytical directions, frequencies and eigendisplacements from the anaddb.nc file specified.
        Non existence of eigendisplacements is accepted for compatibility with abinit 8.0.6
        Raises an error if the other values are not present in anaddb.nc.
        """

        with ETSF_Reader(filepath) as r:
            directions = r.read_value("non_analytical_directions")
            phfreq = r.read_value("non_analytical_phonon_modes")

            #needs a default as the first abinit version including IFCs in the netcdf doesn't have this attribute
            phdispl_cart = r.read_value("non_analytical_phdispl_cart", cmode="c", default=None)

            structure = r.read_structure()

            amu_list = r.read_value("atomic_mass_units", default=None)
            if amu_list is not None:
                atom_species = r.read_value("atomic_numbers")
                amu = {at: a for at, a in zip(atom_species, amu_list)}
            else:
                amu = None

            return cls(structure=structure, directions=directions, phfreqs=phfreq, phdispl_cart=phdispl_cart, amu=amu)

    @lazy_property
    def dyn_mat_eigenvect(self):
        return get_dyn_mat_eigenvec(self.phdispl_cart, self.structure, amu=self.amu)

    @property
    def structure(self):
        """:class:`Structure` object"""
        return self._structure

    def has_direction(self, direction, cartesian=False):
        """
        Checks if the input direction is among those available.

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

        for d in self.directions:
            d = d / np.linalg.norm(d)
            if np.allclose(d, direction):
                return True

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
            structure: :class:`Structure` object.
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
        return len(self.structure)

    @classmethod
    def from_file(cls, filepath):
        """Create the object from a netCDF file."""
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
        """:class:`Structure` object"""
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
            return np.einsum("ktli,ktij,ktuj->ktlu",self.local_vectors,self.ifc_cart_coord,self.local_vectors)

    @lazy_property
    def ifc_local_coord_short_range(self):
        """Short range part of the ifcs in cartesian coordinates"""
        if self.local_vectors is None:
            return None
        else:
            return np.einsum("ktli,ktij,ktuj->ktlu",self.local_vectors,self.ifc_cart_coord_short_range,self.local_vectors)

    @lazy_property
    def ifc_local_coord_ewald(self):
        """Ewald part of the ifcs in local coordinates"""
        return np.einsum("ktli,ktij,ktuj->ktlu",self.local_vectors,self.ifc_cart_coord_ewald,self.local_vectors)

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
            atom_indices =self.structure.indices_from_symbol(atom_element)

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
        be plot at a time.

        Args:
            ifc: an array with shape number_of_atoms*number_of_neighbours of the ifc that should be plotted
            atom_indices: a list of atom indices in the structure. Only neighbours of these atoms will be considered.
            atom_element: symbol of an element in the structure. Only neighbours of these atoms will be considered.
            neighbour_element: symbol of an element in the structure. Only neighbours of this specie will be considered.
            min_dist: minimum distance between atoms and neighbours.
            max_dist: maximum distance between atoms and neighbours.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns:
            matplotlib figure
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        ind = self._filter_ifc_indices(atom_indices=atom_indices, atom_element=atom_element,
                                       neighbour_element=neighbour_element, min_dist=min_dist, max_dist=max_dist)

        dist, filtered_ifc =  self.distances[ind], ifc[ind]

        if 'color' not in kwargs:
            kwargs['color'] = 'blue'

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['lw'] = 0

        ax.set_xlabel('Distance [Bohr]')
        ax.set_ylabel(r'IFC [Ha/Bohr$^2$]')
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
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns:
            matplotlib figure
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
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns:
            matplotlib figure
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
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to blue, symbol to 'o' and lw to 0

        Returns:
            matplotlib figure
        """
        if self.local_vectors is None:
            raise ValueError("Local coordinates are missing. Run anaddb with ifcana=1")

        if self.ifc_local_coord_ewald is None:
            raise ValueError("Ewald contribution is missing, Run anaddb with dipdip=1")

        return self.get_plot_ifc(self.ifc_local_coord_ewald[:, :, 0, 0], atom_indices=atom_indices,
                                 atom_element=atom_element, neighbour_element=neighbour_element, min_dist=min_dist,
                                 max_dist=max_dist, ax=ax, **kwargs)


def get_dyn_mat_eigenvec(phdispl, structure, amu=None):
    """
    Converts the phonon eigendisplacements to the orthonormal eigenvectors of the dynamical matrix.
    Small discrepancies with the original values may be expected due to the different values of the atomic masses in
    abinit and pymatgen.

    Args:
        phdispl: a numpy array containing the eigendisplacements in cartesian coordinates. The last index should have
            size 3*(num atoms), but the rest of the shape is arbitrary. If qpts is not None the first dimension
            should match the q points.
        structure: :class:`Structure` object.
        amu: dictionary that associates the atomic species present in the structure to the values of the atomic
            mass units used for the calculation. If None, values from pymatgen will be used. Note that this will
            almost always lead to inaccuracies in the conversion.

    Returns:
        A numpy array of the same shape as phdispl containing the eigenvectors of the dynamical matrix
    """
    eigvec = np.zeros(np.shape(phdispl), dtype=np.complex)

    if amu is None:
        amu = {e.number: e.atomic_mass for e in structure.composition.elements}

    for j, a in enumerate(structure):
        eigvec[...,3*j:3*(j+1)] = phdispl[...,3*j:3*(j+1)]*np.sqrt(amu[a.specie.number]*abu.amu_emass)/abu.Bohr_Ang

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
