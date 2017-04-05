# coding: utf-8
"""Classes for the analysis of electronic structures."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import copy
import itertools
import json
import pickle
import numpy as np
import pymatgen.core.units as units

from collections import OrderedDict, namedtuple, Iterable
from monty.string import is_string, marquee
from monty.termcolor import cprint
from monty.json import MSONable, MontyEncoder
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
from monty.bisect import find_le, find_gt
from monty.dev import deprecated
from pymatgen.serializers.json_coders import pmg_serialize
from abipy.core.func1d import Function1D
from abipy.core.mixins import NotebookWriter
from abipy.core.kpoints import (Kpoint, KpointList, Kpath, IrredZone, KSamplingInfo, KpointsReaderMixin,
    kmesh_from_mpdivs, Ktables, has_timrev_from_kptopt, map_bz2ibz)
from abipy.core.structure import Structure
from abipy.iotools import ETSF_Reader, bxsf_write
from abipy.tools import gaussian
from abipy.tools.plotting import set_axlims, add_fig_kwargs, get_ax_fig_plt, Marker


import logging
logger = logging.getLogger(__name__)


__all__ = [
    "ElectronBands",
    "ElectronDos",
    "frame_from_ebands",
    "ElectronBandsPlotter",
    "ElectronDosPlotter",
]


class Electron(namedtuple("Electron", "spin kpoint band eig occ kidx")):
    """
    Single-particle state.

    .. Attributes:

        spin: spin index (C convention, i.e >= 0)
        kpoint: :class:`Kpoint` object.
        band: band index. (C convention, i.e >= 0).
        eig: KS eigenvalue.
        occ: Occupation factor.
        kidx: Index of the k-point in the initial array.

    .. note::
        Energies are in eV.
    """
    #def __eq__(self, other):
    #    return (self.spin = other.spin and
    #            self.kpoint == other.kpoint and
    #            self.band == other.band and
    #            self.eig == other.eig
                 # and self.occ == other.occ
    #            )

    #def __ne__(self, other):
    #    return not self == other

    def __str__(self):
        return "spin=%d, kpt=%s, band=%d, eig=%.3f, occ=%.3f" % (
            self.spin, self.kpoint, self.band, self.eig, self.occ)

    @property
    def skb(self):
        """Tuple with (spin, kpoint, band)."""
        return self.spin, self.kpoint, self.band

    def copy(self):
        return Electron(**{f: copy.copy(getattr(self, f)) for f in self._fields})

    @classmethod
    def get_fields(cls, exclude=()):
        fields = list(cls._fields)
        for e in exclude:
            fields.remove(e)

        return tuple(fields)

    def as_dict(self):
        """Convert self into a dict."""
        return super(Electron, self)._asdict()

    def to_strdict(self, fmt=None):
        """Ordered dictionary mapping fields --> strings."""
        d = self.as_dict()
        for k, v in d.items():
            if np.iscomplexobj(v):
                if abs(v.imag) < 1.e-3:
                    d[k] = "%.2f" % v.real
                else:
                    d[k] = "%.2f%+.2fj" % (v.real, v.imag)
            elif isinstance(v, int):
                d[k] = "%d" % v
            else:
                try:
                    d[k] = "%.2f" % v
                except TypeError as exc:
                    #print("k", k, str(exc))
                    d[k] = str(v)
        return d

    @property
    def tips(self):
        """Dictionary with the description of the fields."""
        return self.__class__.TIPS()

    @classmethod
    def TIPS(cls):
        """
        Class method that returns a dictionary with the description of the fields.
        The string are extracted from the class doc string.
        """
        try:
            return cls._TIPS

        except AttributeError:
            # Parse the doc string.
            cls._TIPS = _TIPS = {}
            lines = cls.__doc__.splitlines()

            for i, line in enumerate(lines):
                if line.strip().startswith(".. Attributes"):
                    lines = lines[i+1:]
                    break

            def num_leadblanks(string):
                """Returns the number of the leading whitespaces in a string"""
                return len(string) - len(string.lstrip())

            for field in cls._fields:
                for i, line in enumerate(lines):

                    if line.strip().startswith(field + ":"):
                        nblanks = num_leadblanks(line)
                        desc = []
                        for s in lines[i+1:]:
                            if nblanks == num_leadblanks(s) or not s.strip():
                                break
                            desc.append(s.lstrip())

                        _TIPS[field] = "\n".join(desc)

            diffset = set(cls._fields) - set(_TIPS.keys())
            if diffset:
                raise RuntimeError("The following fields are not documented: %s" % str(diffset))

            return _TIPS


class ElectronTransition(object):
    """
    This object describes an electronic transition between two single-particle states.
    """
    def __init__(self, in_state, out_state):
        """
        Args:
            in_state, out_state: Initial and finale state (:class:`Electron` instances).
        """
        self.in_state = in_state
        self.out_state = out_state

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        app("Energy: %.3f [eV]" % self.energy)
        app("Initial state: %s" % str(self.in_state))
        app("Final state:   %s" % str(self.out_state))

        return "\n".join(lines)

    #def __eq__(self, other):
    #    if other is None: return False
    #    if not isinstance(other, self.__class__): return False
    #    return self.in_state == other.in_state and
    #           self.out_state == other.out_state

    #def __ne__(self, other):
    #    return not self == other

    #def __ge__(self, other):
    #    return self.energy >= other.energy

    @property
    def energy(self):
        """Transition energy in eV."""
        return self.out_state.eig - self.in_state.eig

    @property
    def qpoint(self):
        """k_final - k_initial"""
        return self.out_state.kpoint - self.in_state.kpoint

    @property
    def is_direct(self):
        """True if direct transition."""
        return self.in_state.kpoint == self.out_state.kpoint


class Smearing(AttrDict):
    """
    Stores data and information about the smearing technique.
    """
    _MANDATORY_KEYS = [
        "scheme",
        "occopt",
        "tsmear_ev",
    ]

    @classmethod
    def from_dict(cls, d):
        """
        Makes Smearing obey the general json interface used in pymatgen for easier serialization.
        """
        return cls(**{k: d[k] for k in cls._MANDATORY_KEYS})

    @pmg_serialize
    def as_dict(self):
        """
        Makes Smearing obey the general json interface used in pymatgen for easier serialization.
        """
        return self

    def to_json(self):
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    def __init__(self, *args, **kwargs):
        super(Smearing, self).__init__(*args, **kwargs)
        for mkey in self._MANDATORY_KEYS:
            if mkey not in self:
                raise ValueError("Mandatory key %s must be provided" % str(mkey))

    def __str__(self):
        return "smearing scheme: %s, tsmear_eV: %.3f, occopt: %d" % (self.scheme, self.tsmear_ev, self.occopt)

    @property
    def has_metallic_scheme(self):
        """True if we are using a metallic scheme for occupancies."""
        return self.occopt in [3, 4, 5, 6, 7, 8]


class StatParams(namedtuple("StatParams", "mean stdev min max")):
    """Named tuple with statistical parameters."""
    def __str__(self):
        return "mean = %.3f, stdev = %.3f, min = %.3f, max = %.3f [eV]" % (
            self.mean, self.stdev, self.min, self.max)


class ElectronBandsError(Exception):
    """Exceptions raised by ElectronBands."""


class ElectronBands(object):
    """
    This object stores the electronic band structure.

    .. attribute:: fermie

            Fermi level in eV. Note that, if the band structure has been computed
            with a NSCF run, fermie corresponds to the fermi level obtained
            in the SCF run that produced the density used for the band structure calculation.
    """
    Error = ElectronBandsError

    # FIXME
    # Increase a bit the value of fermie used in bisection routines to solve the problem mentioned below
    pad_fermie = 1e-6
    # One should check whether fermie is recomputed at the end of the SCF cyle
    # I have problems in finding homos/lumos in semiconductors (e.g. Si)
    # because fermie is slightly smaller than the CBM:

    # fermie 5.59845327874
    # homos [Electron(spin=0, kpoint=[0.000, 0.000, 0.000], band=1, eig=5.5984532787385985, occ=2.0)]
    # lumos [Electron(spin=0, kpoint=[0.000, 0.000, 0.000], band=2, eig=5.5984532788661543, occ=2.0)]
    #
    # There's also another possible problem if the DEN is computed on a grid that does not contain the CBM (e.g. Gamma)
    # because the CBM obtained with the NSCF band structure run will be likely above the Ef computed previously.

    @classmethod
    def from_file(cls, filepath):
        """
        Initialize an instance of :class:`ElectronBands` from the netCDF file `filepath`.
        """
        if filepath.endswith(".nc"):
            with ElectronsReader(filepath) as r:
                new = r.read_ebands()
        else:
            raise NotImplementedError("ElectronBands can only be initialized from nc files")

        assert new.__class__ == cls
        return new

    @classmethod
    def from_dict(cls, d):
        d = d.copy()
        kd = d["kpoints"].copy()
        kd.pop("@module")

        kpoints_cls = KpointList.subclass_from_name(kd.pop("@class"))
        #kpoints = kpoints_cls(**kd)
        kpoints = kpoints_cls.from_dict(kd)

        # Needed to support old dictionaries
        if "nspden" not in d: d["nspden"] = 1
        if "nspinor" not in d: d["nspinor"] = 1
        return cls(Structure.from_dict(d["structure"]), kpoints,
                   d["eigens"], d["fermie"], d["occfacts"], d["nelect"], d["nspinor"], d["nspden"],
                   nband_sk=d["nband_sk"], smearing=d["smearing"],
                   #markers=None, widths=None
                   )

    @pmg_serialize
    def as_dict(self):
        """Return dictionary with JSON serialization."""
        return dict(
            structure=self.structure.as_dict(),
            kpoints=self.kpoints.as_dict(),
            eigens=self.eigens.tolist(),
            fermie=float(self.fermie),
            occfacts=self.occfacts.tolist(),
            nelect=float(self.nelect),
            nspinor=self.nspinor,
            nspden=self.nspden,
            nband_sk=self.nband_sk.tolist(),
            smearing=self.smearing.as_dict(),
        )

    @classmethod
    def as_ebands(cls, obj):
        """
        Return an instance of :class:`ElectronBands` from a generic `obj`.
        Supports:

            - instances of cls
            - files (string) that can be open with abiopen and that provide an `ebands` attribute.
            - objects providing an `ebands` attribute
        """
        if isinstance(obj, cls):
            return obj
        elif is_string(obj):
            # path?
            if obj.endswith(".pickle"):
                with open(obj, "rb") as fh:
                    return cls.as_ebands(pickle.load(fh))

            from abipy.abilab import abiopen
            with abiopen(obj) as abifile:
                return abifile.ebands
        elif hasattr(obj, "ebands"):
            # object with ebands
            return obj.ebands

        raise TypeError("Don't know how to extract ebands from object %s" % type(obj))

    def to_json(self):
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    def __init__(self, structure, kpoints, eigens, fermie, occfacts, nelect, nspinor, nspden,
                 nband_sk=None, smearing=None, markers=None, widths=None):
        """
        Args:
            structure: pymatgen structure.
            kpoints: :class:`KpointList` instance.
            eigens: Array-like object with the eigenvalues (eV) stored as [s,k,b]
                    where s: spin , k: kpoint, b: band index
            fermie: Fermi level in eV.
            occfacts: Occupation factors (same shape as eigens)
            nelect: Number of valence electrons in the unit cell.
            nspinor: Number of spinorial components
            nspden: Number of indipendent density components.
            smearing: :class:`Smearing` object storing information on the smearing technique.
            nband_sk: Array-like object with the number of bands treated at each [spin,kpoint]
                      If not given, nband_sk is initialized from eigens.
            markers: Optional dictionary containing markers labelled by a string.
                     Each marker is a list of tuple(x, y, s) where x,and y are the position
                     in the graph and s is the size of the marker.
                     Used for plotting purpose e.g. QP data, energy derivatives...
            widths: Optional dictionary containing data used for the so-called fatbands
                    Each entry is an array of shape [nsppol, nkpt, mband] giving the width
                    of the band at that particular point.
                    Used for plotting purpose e.g. L-projections.
        """
        self.structure = structure

        # Eigenvalues and occupancies are stored in ndarrays ordered by [spin,kpt,band]
        self._eigens = np.atleast_3d(eigens)
        self._occfacts = np.atleast_3d(occfacts)
        assert self._eigens.shape == self._occfacts.shape
        self.nsppol, self.nkpt, self.mband = self.eigens.shape
        self.nspinor, self.nspden = nspinor, nspden

        if nband_sk is not None:
            self.nband_sk = np.array(nband_sk)
        else:
            self.nband_sk = np.array(self.nsppol * self.nkpt * [self.mband])
            self.nband_sk.shape = (self.nsppol, self.nkpt)
        #print(nband_sk)

        self.kpoints = kpoints
        assert self.nkpt == len(self.kpoints)

        self.smearing = {} if smearing is None else smearing
        self.nelect = float(nelect)
        self.fermie = float(fermie)

        # Recompute the Fermi level (in principle should do this only if
        # bands are computed on a BZ mesh with a NSCF run.
        #if self.kpoints.is_ibz: # and iscf < 0
        #    self.recalc_fermie()

        if markers is not None:
            for key, xys in markers.items():
                self.set_marker(key, xys)

        if widths is not None:
            for key, width in widths.items():
                self.set_width(key, width)

    @lazy_property
    def _auto_klabels(self):
        # Find the k-point names in the pymatgen database.
        # We'll use _auto_klabels to label the point in the matplotlib plot
        # if klabels are not specified by the user.
        _auto_klabels = OrderedDict()
        for idx, kpoint in enumerate(self.kpoints):
            name = self.structure.findname_in_hsym_stars(kpoint)
            if name is not None:
                _auto_klabels[idx] = name
                if kpoint.name is None:
                    kpoint.set_name(name)
        return _auto_klabels

    #def __repr__(self):
    #    """String representation (short version)"""

    def __str__(self):
        """
        String representation
        """
        return self.to_string()

    # Handy variables used to loop
    @property
    def spins(self):
        """Spin range"""
        return range(self.nsppol)

    @property
    def nband(self):
        try:
            return self._nband
        except AttributeError:
            assert np.all(self.nband_sk == self.nband_sk[0])
            self._nband = self.nband_sk[0,0]
            return self._nband

    @property
    def kidxs(self):
        """Range with the index of the k-points."""
        return range(self.nkpt)

    @property
    def eigens(self):
        """Eigenvalues in eV. ndarray with shape (nspin, nkpt, mband)."""
        return self._eigens

    @property
    def occfacts(self):
        """Occupation factors. ndarray with shape (nspin, nkpt, mband)."""
        return self._occfacts

    @property
    def reciprocal_lattice(self):
        """Reciprocal lattice vectors in Angstrom."""
        return self.structure.reciprocal_lattice

    @property
    def shape(self):
        """Shape of the array with the eigenvalues."""
        return self.nsppol, self.nkpt, self.mband

    @property
    def has_metallic_scheme(self):
        """True if we are using a metallic scheme for occupancies."""
        return self.smearing.has_metallic_scheme

    def recalc_fermie(self, nelect=None, method="gaussian", step=0.001, width=0.002):
        """
        Recompute the Fermi level.
        """
        if nelect is None: nelect = self.nelect
        edos = self.get_edos(method=method, step=step, width=width)
        ef = edos.find_mu(nelect)
        self.set_fermie(ef)
        return ef

    def set_fermie(self, fermie):
        self.fermie = fermie
        # TODO: Recalculate occupations.

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
            ("nsppol", self.nsppol), ("nkpt", self.nkpt), ("nband", self.nband_sk.min()),
            ("nelect", self.nelect), ("fermie", self.fermie),

        ])
        odict.update(self.structure.get_dict4frame(with_spglib=with_spglib))
        odict.update(self.smearing)

        bws = self.bandwidths
        for spin in self.spins:
            odict["bandwidth_spin%d" % spin] = bws[spin]

        enough_bands = (self.mband > self.nspinor * self.nelect // 2)
        if enough_bands:
            fundamental_gaps = self.fundamental_gaps
            for spin in self.spins:
                odict["fundgap_spin%d" % spin] = fundamental_gaps[spin].energy

            direct_gaps = self.direct_gaps
            for spin in self.spins:
                odict["dirgap_spin%d" % spin] = direct_gaps[spin].energy

        return odict

    @property
    def markers(self):
        try:
            return self._markers
        except AttributeError:
            return {}

    def del_marker(self, key):
        """
        Delete the entry in self.markers with the specied key. All markers are removed if key is None.
        """
        if key is not None:
            try:
                del self._markers[key]
            except AttributeError:
                pass
        else:
            try:
                del self._markers
            except AttributeError:
                pass

    def set_marker(self, key, xys, extend=False):
        """
        Set an entry in the markers dictionary.

        Args:
            key: string used to label the set of markers.
            xys: Three iterables x,y,s where x[i],y[i] gives the
                 positions of the i-th markers in the plot and s[i] is the size of the marker.
            extend: True if the values xys should be added to a pre-existing marker.
        """
        if not hasattr(self, "_markers"):
            self._markers = OrderedDict()

        if extend:
            if key not in self.markers:
                self._markers[key] = Marker(*xys)
            else:
                # Add xys to the previous marker set.
                self._markers[key].extend(*xys)

        else:
            if key in self.markers:
                raise ValueError("Cannot overwrite key %s in data" % key)

            self._markers[key] = Marker(*xys)

    @property
    def widths(self):
        """Widths dictionary"""
        try:
            return self._widths
        except AttributeError:
            return {}

    def del_width(self, key):
        """
        Delete the entry in self.widths with the specified key. All keys are removed if key is None.
        """
        if key is not None:
            try:
                del self._widths[key]
            except AttributeError:
                pass
        else:
            try:
                del self._widths
            except AttributeError:
                pass

    def set_width(self, key, width, overwrite=True):
        """
        Set an entry in the widths dictionary.

        Args:
            key: string used to label the set of markers.
            width: array-like of positive numbers, shape is [nsppol, nkpt, mband].
        """
        width = np.reshape(width, self.shape)

        if not hasattr(self, "_widths"):
            self._widths = OrderedDict()

        if not overwrite and key in self.widths:
            if not np.allclose(width, self.widths[key]):
                raise ValueError("Cannot overwrite key %s in data" % key)

        if np.any(np.iscomplex(width)):
            raise ValueError("Found ambiguous complex entry %s" % str(width))

        if np.any(width < 0.0):
            raise ValueError("Found negative entry in width array %s" % str(width))

        self._widths[key] = width

    @property
    def has_bzmesh(self):
        """True if the k-point sampling is homogeneous."""
        return isinstance(self.kpoints, IrredZone)

    @property
    def has_bzpath(self):
        """True if the bands are computed on a k-path."""
        return isinstance(self.kpoints, Kpath)

    @lazy_property
    def kptopt(self):
        """The value of the kptopt input variable."""
        if hasattr(self, "kpoints.ksampling.kptopt"):
            return self.kpoints.ksampling.kptopt
        else:
            cprint("ebands.kpoints.ksampling.kptopt is not defined, assuming kptopt = 1", "red")
            return 1

    @lazy_property
    def has_timrev(self):
        """True if time-reversal symmetry is used in the BZ sampling."""
        return has_timrev_from_kptopt(self.kptopt)

    def kindex(self, kpoint):
        """
        The index of the k-point in the internal list of k-points.
        Accepts: :class:`Kpoint` instance or integer.
        """
        if isinstance(kpoint, int):
            return kpoint
        else:
            return self.kpoints.index(kpoint)

    def sb_iter(self):
        """Iterator over (spin, band) indices."""
        for spin in self.spins:
            for band in self.nband_sk[spin,k]:
                yield spin, band

    def skb_iter(self):
        """Iterator over (spin, k, band) indices."""
        for spin in self.spins:
            for k in self.kidxs:
                for band in self.nband_sk[spin,k]:
                    yield spin, k, band

    def show_bz(self, **kwargs):
        """Call `matplotlib` to show the Brillouin zone."""
        return self.structure.show_bz(**kwargs)

    def copy(self):
        """Shallow copy of self."""
        return copy.copy(self)

    def deepcopy(self):
        """Deep copy of self."""
        return copy.deepcopy(self)

    def degeneracies(self, spin, kpoint, bands_range, tol_ediff=1.e-3):
        """
        Returns a list with the indices of the degenerate bands.

        Args:
            spin: Spin index.
            kpoint: K-point index or :class:`Kpoint` object
            bands_range: List of band indices to analyze.
            tol_ediff: Tolerance on the energy difference (in eV)

        Returns:
            List of tuples [(e0, bands_e0), (e1, bands_e1, ....]
            Each tuple stores the degenerate energy and a list with the band indices.
            The band structure of silicon at Gamma, for example, will produce something like:
            [(-6.3, [0]), (5.6, [1, 2, 3]), (8.1, [4, 5, 6])]
        """
        # Find the index of the k-point
        k = self.kindex(kpoint)

        # Extract the energies we are interested in.
        bands_list = list(bands_range)
        energies = [self.eigens[spin,k,band] for band in bands_list]

        # Group bands according to their degeneracy.
        bstart, deg_ebands = bands_list[0], []
        e0, bs = energies[0], [bstart]

        for band, e in enumerate(energies[1:]):
            band += (bstart + 1)
            new_deg = abs(e-e0) > tol_ediff

            if new_deg:
                ebs = (e0, bs)
                deg_ebands.append(ebs)
                e0, bs = e, [band]
            else:
                bs.append(band)

        deg_ebands.append((e0, bs))
        return deg_ebands

    def enemin(self, spin=None, band=None):
        """Compute the minimum of the eigenvalues."""
        spin_range = self.spins
        if spin is not None:
            assert isinstance(spin, int)
            spin_range = [spin]

        my_kidxs = self.kidxs

        if band is not None:
            assert isinstance(band, int)
            my_bands = [band]

        emin = np.inf
        for spin in spin_range:
            for k in my_kidxs:
                if band is None:
                    my_bands = range(self.nband_sk[spin,k])
                for band in my_bands:
                    e = self.eigens[spin,k,band]
                    emin = min(emin, e)
        return emin

    def enemax(self, spin=None, band=None):
        """Compute the maximum of the eigenvalues."""
        spin_range = self.spins
        if spin is not None:
            assert isinstance(spin, int)
            spin_range = [spin]

        my_kidxs = self.kidxs

        if band is not None:
            assert isinstance(band, int)
            my_bands = [band]

        emax = -np.inf
        for spin in spin_range:
            for k in my_kidxs:
                if band is None:
                    my_bands = range(self.nband_sk[spin,k])
                for band in my_bands:
                    e = self.eigens[spin,k,band]
                    emax = max(emax, e)
        return emax

    def dispersionless_states(self, erange=None, deltae=0.05, kfact=0.9):
        """
        This function detects dispersionless states.
        A state is dispersionless if there are more that (nkpt * kfact) energies
        in the energy intervale [e0-deltae, e0+deltae]

        Args:
            erange=Energy range to be analyzed in the form [emin, emax]
            deltae: Defines the energy interval in eV around the KS eigenvalue.
            kfact: Can be used to change the criterion used to detect dispersionless states.

        Returns:
            List of :class:`Electron` objects. Each item contains information on
            the energy, the occupation and the location of the dispersionless band.
        """
        if erange is None: erange = [-np.inf, np.inf]
        kref = 0
        dless_states = []
        for spin in self.spins:
            for band in range(self.nband_sk[spin,kref]):
                e0 = self.eigens[spin, kref, band]
                if not erange[1] > e0 > erange[0]: continue
                hrange = [e0 - deltae, e0 + deltae]
                hist, bin_hedges = np.histogram(self.eigens[spin,:,:],
                    bins=2, range=hrange, weights=None, density=False)
                #print("hist", hist, "hrange", hrange, "bin_hedges", bin_hedges)

                if hist.sum() > self.nkpt * kfact:
                    state = self._electron_state(spin, kref, band)
                    dless_states.append(state)

        return dless_states

    def raw_print(self, stream=sys.stdout):
        """Print k-points and energies on stream."""
        stream.write("# Band structure energies in Ev.\n")
        stream.write("# idx   kpt_red(1:3)  ene(b1) ene(b2) ...\n")

        fmt_k = lambda k: " %.6f" % k
        fmt_e = lambda e: " %.6f" % e
        for spin in self.spins:
            stream.write("# spin = " + str(spin) + "\n")
            for k, kpoint in enumerate(self.kpoints):
                nb = self.nband_sk[spin,k]
                ene_sk = self.eigens[spin,k,:nb]
                st = str(k+1)
                for c in kpoint: st += fmt_k(c)
                for e in ene_sk: st += fmt_e(e)
                stream.write(st+"\n")

        stream.flush()

    def to_dataframe(self, e0="fermie"):
        """
        Return a pandas DataFrame with the following columns:

          ['spin', 'kidx', 'band', 'eig', 'occ', 'kpoint']

        where:

        ==============  ==========================
        Column          Meaning
        ==============  ==========================
        spin            spin index
        kidx            k-point index
        band            band index
        eig             KS eigenvalue in eV.
        occ             Occupation of the state.
        kpoint          :class:`Kpoint` object
        ==============  ==========================

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
                The Fermi energy is saved in frame.fermie
        """
        import pandas as pd
        rows = []
        e0 = self.get_e0(e0)
        for spin in self.spins:
            for k, kpoint in enumerate(self.kpoints):
                for band in range(self.nband_sk[spin,k]):
                    eig = self.eigens[spin,k,band] - e0
                    rows.append(OrderedDict([
                               ("spin", spin),
                               ("kidx", k),
                               ("band", band),
                               ("eig", eig),
                               ("occ", self.occfacts[spin, k, band]),
                               ("kpoint", self.kpoints[k]),
                            ]))

        frame = pd.DataFrame(rows, columns=list(rows[0].keys()))
        frame.fermie = e0
        return frame

    # Alias to maintain compatibility
    # TODO: Remove it in 0.4
    to_pdframe = to_dataframe

    @add_fig_kwargs
    def boxplot(self, ax=None, e0="fermie", brange=None, swarm=False, **kwargs):
        """
        Use seaborn to draw a box plot to show distributions of eigenvalues with respect to the band index.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            brange: Only bands such as `brange[0] <= band_index < brange[1]` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            kwargs: Keyword arguments passed to seaborn boxplot.
        """
        # Get the dataframe and select bands
        frame = self.to_dataframe(e0=e0)
        if brange is not None: frame = frame[brange[0] <= frame["band"] < brange[1]]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        import seaborn.apionly as sns
        hue = None if self.nsppol == 1 else "spin"
        ax = sns.boxplot(x="band", y="eig", data=frame, hue=hue, ax=ax, **kwargs)
        if swarm:
            sns.swarmplot(x="band", y="eig", data=frame, hue=hue, color=".25", ax=ax)
        return fig

    def to_pymatgen(self):
        """
        Return a pymatgen bandstructure object.
        """
        from pymatgen.electronic_structure.core import Spin
        from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine

        assert np.all(self.nband_sk == self.nband_sk[0,0])

        # eigenvals is a dict of energies for spin up and spin down
        # {Spin.up:[][],Spin.down:[][]}, the first index of the array
        # [][] refers to the band and the second to the index of the
        # kpoint. The kpoints are ordered according to the order of the
        # kpoints array. If the band structure is not spin polarized, we
        # only store one data set under Spin.up
        eigenvals = {Spin.up: self.eigens[0,:,:].T.copy().tolist()}
        if self.nsppol == 2:
            eigenvals[Spin.down] = self.eigens[1,:,:].T.copy().tolist()

        if self.kpoints.is_path:
            labels_dict = {k.name: k.frac_coords for k in self.kpoints if k.name is not None}
            return BandStructureSymmLine(self.kpoints.frac_coords, eigenvals, self.reciprocal_lattice, self.fermie,
                                         labels_dict, coords_are_cartesian=False, structure=self.structure, projections=None)
        else:
            return BandStructure(self.kpoints.frac_coords, eigenvals, self.reciprocal_lattice, self.fermie,
                                 labels_dict=None, coords_are_cartesian=False, structure=self.structure, projections=None)

    def _electron_state(self, spin, kpoint, band):
        """
        Build an instance of :class:`Electron` from the spin, kpoint and band index"""
        kidx = self.kindex(kpoint)
        eig = self.eigens[spin, kidx, band]
        return Electron(spin=spin,
                        kpoint=self.kpoints[kidx],
                        band=band,
                        eig=eig,
                        occ=self.occfacts[spin, kidx, band],
                        kidx=kidx,
                        #fermie=self.fermie
                        )

    @property
    def lomos(self):
        """lomo states for each spin channel as a list of nsppol :class:`Electron`."""
        lomos = self.nsppol * [None]
        for spin in self.spins:
            lomo_kidx = self.eigens[spin,:,0].argmin()
            lomos[spin] = self._electron_state(spin, lomo_kidx, 0)

        return lomos

    def lomo_sk(self, spin, kpoint):
        """
        Returns the LOMO state for the given spin, kpoint.

        Args:
            spin: Spin index
            kpoint: Index of the kpoint or :class:`Kpoint` object.
        """
        return self._electron_state(spin, kpoint, 0)

    def homo_sk(self, spin, kpoint):
        """
        Returns the HOMO state for the given spin, kpoint.

        Args:
            spin: Spin index
            kpoint: Index of the kpoint or :class:`Kpoint` object.
        """
        k = self.kindex(kpoint)
        # Find rightmost value less than or equal to fermie.
        b = find_le(self.eigens[spin,k,:], self.fermie + self.pad_fermie)
        return self._electron_state(spin, k, b)

    def lumo_sk(self, spin, kpoint):
        """
        Returns the LUMO state for the given spin, kpoint.

        Args:
            spin: Spin index
            kpoint: Index of the kpoint or :class:`Kpoint` object.
        """
        k = self.kindex(kpoint)
        # Find leftmost value greater than fermie.
        b = find_gt(self.eigens[spin,k,:], self.fermie + self.pad_fermie)
        return self._electron_state(spin, k, b)

    @property
    def homos(self):
        """
        homo states for each spin channel as a list of nsppol :class:`Electron`.
        """
        homos = self.nsppol * [None]

        for spin in self.spins:
            blist, enes = [], []
            for k in self.kidxs:
                # Find rightmost value less than or equal to fermie.
                b = find_le(self.eigens[spin,k,:], self.fermie + self.pad_fermie)
                blist.append(b)
                enes.append(self.eigens[spin,k,b])

            homo_kidx = np.array(enes).argmax()
            homo_band = blist[homo_kidx]

            # Build Electron instance.
            homos[spin] = self._electron_state(spin, homo_kidx, homo_band)

        return homos

    @property
    def lumos(self):
        """
        lumo states for each spin channel as a list of nsppol :class:`Electron`.
        """
        lumos = self.nsppol * [None]

        for spin in self.spins:
            blist, enes = [], []
            for k in self.kidxs:
                # Find leftmost value greater than fermie.
                b = find_gt(self.eigens[spin,k,:], self.fermie + self.pad_fermie)
                blist.append(b)
                enes.append(self.eigens[spin,k,b])

            lumo_kidx = np.array(enes).argmin()
            lumo_band = blist[lumo_kidx]

            # Build Electron instance.
            lumos[spin] = self._electron_state(spin, lumo_kidx, lumo_band)

        return lumos

    #def is_metal(self, spin)
    #    """True if this spin channel is metallic."""
    #    if not self.has_metallic_scheme: return False
    #    for k in self.kidxs:
    #        # Find leftmost value greater than x.
    #        b = find_gt(self.eigens[spin,k,:], self.fermie)
    #        if self.eigens[spin,k,b] < self.fermie + 0.01:
    #            return True

    #def is_semimetal(self, spin)
    #    """True if this spin channel is semi-metal."""
    #    fun_gaps = self.fundamental_gaps
    #    for spin in self.spins:
    #       if abs(fun_gaps.ene) <  TOL_EGAP

    @lazy_property
    def bandwidths(self):
        """The bandwidth for each spin channel i.e. the energy difference (homo - lomo)."""
        return [self.homos[spin].eig - self.lomos[spin].eig for spin in self.spins]

    @lazy_property
    def fundamental_gaps(self):
        """List of :class:`ElectronTransition` with info on the fundamental gaps for each spin."""
        return [ElectronTransition(self.homos[spin], self.lumos[spin]) for spin in self.spins]

    @lazy_property
    def direct_gaps(self):
        """List of :class:`ElectronTransition` with info on the direct gaps for each spin."""
        dirgaps = self.nsppol * [None]
        for spin in self.spins:
            gaps = []
            for k in self.kidxs:
                homo_sk = self.homo_sk(spin, k)
                lumo_sk = self.lumo_sk(spin, k)
                gaps.append(lumo_sk.eig - homo_sk.eig)

            # Find the index of the k-point where the direct gap is located.
            kdir = np.array(gaps).argmin()
            dirgaps[spin] = ElectronTransition(self.homo_sk(spin, kdir), self.lumo_sk(spin, kdir))

        return dirgaps

    def to_string(self, title=None, with_structure=True, with_kpoints=True, **kwargs):
        """
        Human-readable string with useful info such as band gaps, position of HOMO, LOMO...

        Args:
            with_structure: False if structural info shoud not be displayed.
            with_kpoints: False if k-point info shoud not be displayed.
        """
        lines = []; app = lines.append
        if title is not None:
            app(marquee(title, mark="="))

        if with_structure:
            app(str(self.structure))
            app("")

        app("Number of electrons: %s, Fermi level: %.3f [eV]" % (self.nelect, self.fermie))
        app("nsppol: %d, nkpt: %d, mband: %d, nspinor: %s, nspden: %s" % (
           self.nsppol, self.nkpt, self.mband, self.nspinor, self.nspden))
        app(str(self.smearing))

        def indent(s):
            return "    " + s.replace("\n", "\n    ")

        if not self.has_metallic_scheme:
            enough_bands = (self.mband > self.nspinor * self.nelect // 2)
            for spin in self.spins:
                if self.nsppol == 2:
                    app(">>> For spin %s" % spin)
                if enough_bands:
                    app("Direct gap:\n%s" % indent(str(self.direct_gaps[spin])))
                    app("Fundamental gap:\n%s" % indent(str(self.fundamental_gaps[spin])))
                app("Bandwidth: %.3f [eV]" % self.bandwidths[spin])
                app("Valence minimum located at:\n%s" % indent(str(self.lomos[spin])))
                app("Valence max located at:\n%s" % indent(str(self.homos[spin])))
                app("")

        if with_kpoints:
            app(marquee("K-points", mark="="))
            app(str(self.kpoints))
            app("")

        return "\n".join(lines)

    def spacing(self, axis=None):
        """
        Compute the statistical parameters of the energy spacing, i.e. e[b+1] - e[b]

        Returns:
            `namedtuple` with the statistical parameters in eV
        """
        ediff = self.eigens[1:,:,:] - self.eigens[:self.mband-1,:,:]

        return StatParams(
            mean=ediff.mean(axis=axis),
            stdev=ediff.std(axis=axis),
            min=ediff.min(axis=axis),
            max=ediff.max(axis=axis))

    def statdiff(self, other, axis=None, numpy_op=np.abs):
        """
        Compare the eigenenergies of two bands and compute the
        statistical parameters: mean, standard deviation, min and max
        The bands are aligned wrt to their fermi level.

        Args:
            other: :class:`ElectronBands` object.
            axis:  Axis along which the statistical parameters are computed.
                   The default is to compute the parameters of the flattened array.
            numpy_op: Numpy function to apply to the difference of the eigenvalues. The
                      default computes `|self.eigens - other.eigens|`.

        Returns:
            `namedtuple` with the statistical parameters in eV
        """
        ediff = numpy_op(self.eigens - self.fermie - other.eigens + other.fermie)
        return StatParams(mean=ediff.mean(axis=axis),
                          stdev=ediff.std(axis=axis),
                          min=ediff.min(axis=axis),
                          max=ediff.max(axis=axis)
                          )

    def ipw_dos(self):
        """
        Return an ipython widget with controllers to compute the electron DOS.
        """
        import ipywidgets as ipw

        def plot_dos(method, step, width):
            edos = self.get_edos(method=method, step=step, width=width)
            edos.plot()

        return ipw.interactive(
                plot_dos,
                method=["gaussian", "tetra"],
                step=ipw.FloatSlider(value=0.1, min=1e-6, max=1, step=0.05, description="Step of linear mesh [eV]"),
                width=ipw.FloatSlider(value=0.2, min=1e-6, max=1, step=0.05, description="Gaussian broadening [eV]"),
                __manual=True,
            )

    def get_edos(self, method="gaussian", step=0.1, width=0.2):
        """
        Compute the electronic DOS on a linear mesh.

        Args:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Returns:
            :class:`ElectronDos` object.
        """
        # Weights must be normalized to one.
        wsum = self.kpoints.sum_weights()
        if abs(wsum - 1) > 1.e-6:
            err_msg = "Kpoint weights should sum up to one while sum_weights is %.3f\n" % wsum
            err_msg += "The list of kpoints does not represent a homogeneous sampling of the BZ\n"
            err_msg += str(type(self.kpoints)) # + "\n" + str(self.kpoints)
            raise ValueError(err_msg)

        # Compute the linear mesh.
        epad = 3.0 * width
        e_min = self.enemin() - epad
        e_max = self.enemax() + epad

        nw = int(1 + (e_max - e_min) / step)
        mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)
        dos = np.zeros((self.nsppol, nw))

        # TODO: Write cython version.
        if method == "gaussian":
            for spin in self.spins:
                for k, kpoint in enumerate(self.kpoints):
                    weight = kpoint.weight
                    for band in range(self.nband_sk[spin,k]):
                        e = self.eigens[spin,k,band]
                        dos[spin] += weight * gaussian(mesh, width, center=e)

        else:
            raise ValueError("Method %s is not supported" % method)

        # Use fermie from Abinit if we are not using metallic scheme for occopt.
        fermie = None
        #if self.smearing["occopt"] == 1:
        #    print("using fermie from GSR")
        #    fermie = self.fermie
        edos = ElectronDos(mesh, dos, self.nelect, fermie=fermie)
        #print("ebands.fermie", self.fermie, "edos.fermie", edos.fermie)
        return edos

    def get_ejdos(self, spin, valence, conduction, method="gaussian", step=0.1, width=0.2, mesh=None):
        r"""
        Compute the join density of states at q == 0.
            :math:`\sum_{kbv} f_{vk} (1 - f_{ck}) \delta(\omega - E_{ck} + E_{vk})`

        Args:
            spin: Spin index.
            valence: Int or iterable with the valence indices.
            conduction: Int or iterable with the conduction indices.
            method: String defining the method.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            mesh: Frequency mesh to use. If None, the mesh is computed automatically from the eigenvalues.

        Returns:
            :class:`Function1D` object.
        """
        wsum = self.kpoints.sum_weights()
        if abs(wsum - 1) > 1.e-6:
            err_msg =  "Kpoint weights should sum up to one while sum_weights is %.3f\n" % wsum
            err_msg += "The list of kpoints does not represent a homogeneous sampling of the BZ\n"
            err_msg += str(type(self.kpoints)) + "\n" + str(self.kpoints)
            raise ValueError(err_msg)

        if not isinstance(valence, Iterable): valence = [valence]
        if not isinstance(conduction, Iterable): conduction = [conduction]

        if mesh is None:
            # Compute the linear mesh.
            cmin, cmax = +np.inf, -np.inf
            vmin, vmax = +np.inf, -np.inf
            for c in conduction:
                cmin = min(cmin, self.eigens[spin,:,c].min())
                cmax = max(cmax, self.eigens[spin,:,c].max())
            for v in valence:
                vmin = min(vmin, self.eigens[spin,:,v].min())
                vmax = max(vmax, self.eigens[spin,:,v].max())

            e_min = cmin - vmax
            e_min -= 0.1 * abs(e_min)
            e_max = cmax - vmin
            e_max += 0.1 * abs(e_max)

            nw = int(1 + (e_max - e_min) / step)
            mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)
        else:
            nw = len(mesh)

        jdos = np.zeros(nw)

        # Normalize the occupation factors.
        full = 2.0 if self.nsppol == 1 else 1.0

        if method == "gaussian":
            for k, kpoint in enumerate(self.kpoints):
                weight = kpoint.weight
                for c in conduction:
                    ec = self.eigens[spin,k,c]
                    fc = 1.0 - self.occfacts[spin,k,c] / full
                    for v in valence:
                        ev = self.eigens[spin,k,v]
                        fv = self.occfacts[spin,k,v] / full
                        fact = weight * fv * fc
                        jdos += fact * gaussian(mesh, width, center=ec-ev)

        else:
            raise ValueError("Method %s is not supported" % method)

        return Function1D(mesh, jdos)

    @add_fig_kwargs
    def plot_ejdosvc(self, vrange, crange, method="gaussian", step=0.1, width=0.2,
                     cumulative=True, ax=None, alpha=0.7, **kwargs):
        """
        Plot the decomposition of the joint-density of States (JDOS).

        Args:
            vrange: Int or `Iterable` with the indices of the valence bands to consider.
            crange: Int or `Iterable` with the indices of the conduction bands to consider.
            method: String defining the method.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            cumulative: True for cumulative plots (default).
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        if not isinstance(crange, Iterable): crange = [crange]
        if not isinstance(vrange, Iterable): vrange = [vrange]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)
        ax.set_xlabel('Energy [eV]')
        cmap = plt.get_cmap("jet")
        lw = 1.0

        for s in self.spins:
            spin_sign = +1 if s == 0 else -1

            # Get total JDOS for this spin
            tot_jdos = spin_sign * self.get_ejdos(s, vrange, crange, method=method, step=step, width=width)

            # Decomposition in terms of v --> c transitions.
            jdos_vc = OrderedDict()
            for v in vrange:
                for c in crange:
                    jd = self.get_ejdos(s, v, c, method=method, step=step, width=width, mesh=tot_jdos.mesh)
                    jdos_vc[(v, c)] = spin_sign * jd

            # Plot data for this spin.
            if cumulative:
                cumulative = np.zeros(len(tot_jdos))
                num_plots, i = len(jdos_vc), 0
                for (v, c), jdos in jdos_vc.items():
                    label = r"$v=%s \rightarrow c=%s, \sigma=%s$" % (v, c, s)
                    color = cmap(float(i)/num_plots)
                    x, y = jdos.mesh, jdos.values
                    ax.plot(x, cumulative + y, lw=lw, label=label, color=color)
                    ax.fill_between(x, cumulative, cumulative + y, facecolor=color, alpha=alpha)
                    cumulative += jdos.values
                    i += 1
            else:
                num_plots, i = len(jdos_vc), 0
                for (v, c), jdos in jdos_vc.items():
                    color = cmap(float(i)/num_plots)
                    jdos.plot_ax(ax, color=color, lw=lw, label=r"$v=%s \rightarrow c=%s, \sigma=%s$" % (v, c, s))
                    i += 1

            tot_jdos.plot_ax(ax, color="k", lw=lw, label=r"Total JDOS, $\sigma=%s$" % s)

        ax.legend(loc="best")
        return fig

    def apply_scissors(self, scissors):
        """
        Modify the band structure with the scissors operator.

        Args:
            scissors: An instance of :class:`Scissors`.

        Returns:
            New instance of :class:`ElectronBands` with modified energies.
        """
        if self.nsppol == 1 and not isinstance(scissors, Iterable):
            scissors = [scissors]
        if self.nsppol == 2 and len(scissors) != 2:
            raise ValueError("Expecting two scissors operators for spin up and down")

        # Create new array with same shape as self.
        qp_energies = np.zeros(self.shape)

        # Calculate quasi-particle energies with the scissors operator.
        for spin in self.spins:
            sciss = scissors[spin]
            for k in self.kidxs:
                for band in range(self.nband_sk[spin,k]):
                    e0 = self.eigens[spin,k,band]
                    try:
                        qp_ene = e0 + sciss.apply(e0)
                    except sciss.Error:
                        raise

                    # Update the energy.
                    qp_energies[spin,k,band] = qp_ene

        # Apply the scissors to the Fermi level as well.
        # NB: This should be ok for semiconductors in which fermie == CBM (abinit convention)
        # and there's usually one CBM state whose QP correction is expected to be reproduced
        # almost exactly by the polyfit.
        # Not sure about metals. Besides occupations are not changed here!
        fermie = self.fermie + scissors[0].apply(self.fermie)
        #fermie = self.fermie
        print("KS fermie", self.fermie, "--> QP fermie", fermie, "Delta(QP-KS)=", fermie - self.fermie)

        return self.__class__(
            self.structure, self.kpoints, qp_energies, fermie, self.occfacts, self.nelect, self.nspinor, self.nspden,
            nband_sk=self.nband_sk, smearing=self.smearing, markers=self.markers)

    @add_fig_kwargs
    def plot(self, ax=None, klabels=None, band_range=None, e0="fermie",
             ylims=None, marker=None, width=None, **kwargs):
        r"""
        Plot the band structure.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            klabels: dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels. e.g.
                klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0):"L"}.
            band_range: Tuple specifying the minimum and maximum band to plot (default: all bands are plotted)
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            marker: String defining the marker to plot. Accepts the syntax `markername:fact` where
                fact is a float used to scale the marker size.
            width: String defining the width to plot. Accepts the syntax widthname:fact where
                fact is a float used to scale the stripe size.

        Returns:
            `matplotlib` figure
        """
        # Select the band range.
        if band_range is None:
            band_range = range(self.mband)
        else:
            band_range = range(band_range[0], band_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, klabels=klabels)
        set_axlims(ax, ylims, "y")

        # Plot the band energies.
        for spin in self.spins:
            if spin == 0:
                opts = {"color": "black", "linewidth": 2.0}
            else:
                opts = {"color": "red", "linewidth": 2.0}

            for band in band_range:
                self.plot_ax(ax, e0, spin=spin, band=band, **opts)

        # Add markers to the plot.
        if marker is not None:
            try:
                key, fact = marker.split(":")
            except ValueError:
                key = marker
                fact = 1
            fact = float(fact)

            self.plot_marker_ax(ax, key, e0, fact=fact)

        # Plot fatbands.
        if width is not None:
            try:
                key, fact = width.split(":")
            except ValueError:
                key = width
                fact = 1

            self.plot_width_ax(ax, key, e0, fact=fact)

        return fig

    def decorate_ax(self, ax, **kwargs):
        """
        Decoracte matplotlib Axis

        Accept:
            title:
            klabels
            klabel_size:
        """
        title = kwargs.pop("title", None)
        if title is not None: ax.set_title(title)

        ax.grid(True)
        ax.set_ylabel('Energy [eV]')

        # Set ticks and labels.
        klabels = kwargs.pop("klabels", None)
        ticks, labels = self._make_ticks_and_labels(klabels)
        if ticks:
            # Don't show label if previous k-point is the same.
            for il in range(1, len(labels)):
                if labels[il] == labels[il-1]: labels[il] = ""
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))
            ax.set_xlim(0, ticks[-1])

    def get_e0(self, e0):
        """
        e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
        """
        if e0 is None:
            return 0.0
        elif is_string(e0):
            if e0 == "fermie":
                return self.fermie
            elif e0 == "None":
                return 0.0
            else:
                raise ValueError("Wrong value for e0: %s" % e0)
        else:
            # Assume number
            return e0

    def plot_ax(self, ax, e0, spin=None, band=None, **kwargs):
        """
        Helper function to plot the energies for (spin, band) on the axis ax.

        Args:
            ax: Matplotlib axis.
            e0: Option used to define the zero of energy in the band structure plot.
            spin: Spin index. If None, all spins are plotted.
            band: Band index, If None, all bands are plotted.

        Return matplotlib lines
        """
        spin_range = range(self.nsppol) if spin is None else [spin]
        band_range = range(self.mband) if band is None else [band]

        # Disable labels.
        if "label" not in kwargs:
            kwargs["label"] = "_no_legend_" # Actively suppress.

        xx, lines = range(self.nkpt), []
        e0 = self.get_e0(e0)
        for spin in spin_range:
            for band in band_range:
                yy = self.eigens[spin,:,band] - e0
                lines.extend(ax.plot(xx, yy, **kwargs))

        return lines

    def plot_width_ax(self, ax, key, e0, spin=None, band=None, fact=1.0, **kwargs):
        """
        Helper function to plot fatbands for the given (spin,band) on the axis ax.

        Args:
            e0: Option used to define the zero of energy in the band structure plot.
        """
        spin_range = range(self.nsppol) if spin is None else [spin]
        band_range = range(self.mband) if band is None else [band]

        facecolor = kwargs.pop("facecolor", "blue")
        alpha = kwargs.pop("alpha", 0.7)

        x, width = range(self.nkpt), fact * self.widths[key]

        e0 = self.get_e0(e0)
        for spin in spin_range:
            for band in band_range:
                y, w = self.eigens[spin,:,band] - e0, width[spin,:,band] * fact
                ax.fill_between(x, y-w/2, y+w/2, facecolor=facecolor, alpha=alpha)

    def plot_marker_ax(self, ax, key, e0, fact=1.0):
        """
        Helper function to plot the markers on the axes ax.

        Args:
            e0: Option used to define the zero of energy in the band structure plot.
        """
        pos, neg = self.markers[key].posneg_marker()
        e0 = self.get_e0(e0)

        # Use different symbols depending on the value of s. Cannot use negative s.
        if pos:
            ax.scatter(pos.x, pos.y - e0, s=np.abs(pos.s)*fact, marker="^", label=key + " >0")

        if neg:
            ax.scatter(neg.x, neg.y - e0, s=np.abs(neg.s)*fact, marker="v", label=key + " <0")

    def _make_ticks_and_labels(self, klabels):
        """Return ticks and labels from the mapping qlabels."""
        if klabels is not None:
            d = OrderedDict()
            for kcoord, kname in klabels.items():
                # Build Kpoint instance.
                ktick = Kpoint(kcoord, self.reciprocal_lattice)
                for idx, kpt in enumerate(self.kpoints):
                    if ktick == kpt: d[idx] = kname

        else:
            d = self._auto_klabels

        # Return ticks, labels
        return list(d.keys()), list(d.values())

    @add_fig_kwargs
    def plot_with_edos(self, dos, klabels=None, axlist=None, e0="fermie", ylims=None, **kwargs):
        r"""
        Plot the band structure and the DOS.

        Args:
            dos: An instance of :class:`ElectronDos`.
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. `klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}`.
            axlist: The axes for the bandstructure plot and the DOS plot. If axlist is None, a new figure
                is created and the two axes are automatically generated.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues and the DOS to have zero energy at the Fermi energy.
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   This should be ok in semiconductors. In metals, however, a better value of the Fermi energy
                   can be obtained from the DOS provided that the k-sampling for the DOS is much denser than
                   the one used to compute the density. See `edos_fermie`.
                - `edos_fermie`: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0

        Returns:
            `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        if axlist is None:
            # Build axes and align bands and DOS.
            gspec = GridSpec(1, 2, width_ratios=[2, 1])
            gspec.update(wspace=0.05)
            ax1 = plt.subplot(gspec[0])
            ax2 = plt.subplot(gspec[1], sharey=ax1)
        else:
            # Take them from axlist.
            ax1, ax2 = axlist

        # Define the zero of energy.
        e0 = self.get_e0(e0) if e0 != "edos_fermie" else dos.fermie
        if not kwargs: kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the band structure
        for spin in self.spins:
            for band in range(self.mband):
                self.plot_ax(ax1, e0, spin=spin, band=band, **kwargs)

        self.decorate_ax(ax1, klabels=klabels)
        set_axlims(ax1, ylims, "y")

        # Plot the DOS
        dos.plot_ax(ax2, e0, exchange_xy=True, **kwargs)
        ax2.grid(True)
        ax2.yaxis.set_ticks_position("right")
        ax2.yaxis.set_label_position("right")
        set_axlims(ax2, ylims, "y")

        fig = plt.gcf()
        return fig

    def to_xmgrace(self, filepath):
        """
        Write xmgrace file with band structure energies and labels for high-symmetry k-points.
        """
        f = open(filepath, "wt")
        def w(s):
            f.write(s)
            f.write("\n")

        emef = np.array(self.eigens - self.fermie)

        import datetime
        w("# Grace project file")
        w("# Generated by abipy on: %s" % str(datetime.datetime.today()))
        w("# Crystalline structure:")
        for s in str(self.structure).splitlines():
            w("# %s" % s)
        w("# mband: %d, nkpt: %d, nsppol: %d, nspinor: %d" % (self.mband, self.nkpt, self.nsppol, self.nspinor))
        w("# nelect: %8.2f, %s" % (self.nelect, str(self.smearing)))
        w("# Energies are in eV. Zero set to efermi, previously it was at: %s [eV]" % self.fermie)
        w("# List of k-points and their index (C notation i.e. count from 0)")
        for ik, kpt in enumerate(self.kpoints):
            w("# %d %s" % (ik, str(kpt.frac_coords)))
        w("@page size 792, 612")
        w("@page scroll 5%")
        w("@page inout 5%")
        w("@link page off")
        w("@with g0")
        w("@world xmin 0.00")
        w('@world xmax %d' % (self.nkpt-1))
        w('@world ymin %s' % emef.min())
        w('@world ymax %s' % emef.max())
        w('@default linewidth 1.5')
        w('@xaxis  tick on')
        w('@xaxis  tick major 1')
        w('@xaxis  tick major color 1')
        w('@xaxis  tick major linestyle 3')
        w('@xaxis  tick major grid on')
        w('@xaxis  tick spec type both')
        w('@xaxis  tick major 0, 0')

        kticks, klabels = self._make_ticks_and_labels(klabels=None)
        w('@xaxis  tick spec %d' % len(kticks))
        for ik, (ktick, klabel) in enumerate(zip(kticks, klabels)):
            w('@xaxis  tick major %d, %d' % (ik, ktick))
            w('@xaxis  ticklabel %d, "%s"' % (ik, klabel))

        w('@xaxis  ticklabel char size 1.500000')
        w('@yaxis  tick major 10')
        w('@yaxis  label "Band Energy [eV]"')
        w('@yaxis  label char size 1.500000')
        w('@yaxis  ticklabel char size 1.500000')
        ii = -1
        for spin in range(self.nsppol):
            for band in range(self.mband):
                ii += 1
                w('@    s%d line color %d' % (ii, spin + 1))

        ii = -1
        for spin in range(self.nsppol):
            for band in range(self.mband):
                ii += 1
                w('@target G0.S%d' % ii)
                w('@type xy')
                for ik in range(self.nkpt):
                    w('%d %.8E' % (ik, emef[spin, ik, band]))
                w('&')

        f.close()

    def to_bxsf(self, filepath):
        """
        Export the full band structure to `filepath` in BXSF format
        suitable for the visualization of the Fermi surface with Xcrysden (xcrysden --bxsf FILE).
        Require k-points in IBZ and gamma-centered k-mesh.
        """
        # Sanity check.
        errors = []; eapp = errors.append
        if np.any(self.nband_sk != self.nband_sk[0,0]):
            eapp("The number of bands in nband must be constant")
        if not self.kpoints.is_ibz:
            eapp("Expecting an IBZ sampling for the Fermi surface but got %s" % type(ebands.kpoints))
        if not self.kpoints.is_mpmesh:
            eapp("Monkhorst-Pack meshes are required for bxsf output.")

        mpdivs, shifts = self.kpoints.mpdivs_shifts
        if shifts is not None and not np.all(shifts == 0.0):
            eapp("Gamma-centered k-meshes are required by Xcrysden.")
        if errors:
            raise ValueError("\n".join(errors))

        # Xcrysden requires points in the unit cell (C-order)
        # and the mesh must include the periodic images hence pbc=True.
        bz2ibz = map_bz2ibz(self.structure, self.kpoints.frac_coords, mpdivs, self.has_timrev, pbc=True)

        # Construct bands in BZ: e_{TSk} = e_{k}
        len_bz = len(bz2ibz)
        emesh_sbk = np.empty((self.nsppol, self.nband, len_bz))
        for ik_bz in range(len_bz):
            ik_ibz = bz2ibz[ik_bz]
            emesh_sbk[:, :, ik_bz] = self.eigens[:, ik_ibz, :]

        # Write BXSF file.
        with open(filepath, "wt") as fh:
            bxsf_write(fh, self.structure, self.nsppol, self.nband, mpdivs+1, emesh_sbk, self.fermie, unit="eV")

    def derivatives(self, spin, band, order=1, acc=4, asmarker=None):
        """
        Compute the derivative of the eigenvalues wrt to k.

        Args:
            spin: Spin index
            band: Band index
            order:
            acc:
            asmarker:

        Returns:
        """
        if self.kpoints.is_path:
            # Extract the branch.
            branch = self.eigens[spin, :, band]
            # Simulate free-electron bands. This will produce all(effective masses == 1)
            #branch = 0.5 * units.Ha_to_eV * np.array([(k.norm * units.bohr_to_ang)**2 for k in self.kpoints])

            # Compute derivatives by finite differences.
            ders_onlines = self.kpoints.finite_diff(branch, order=order, acc=acc)

            if asmarker is not None:
                x, y, s = [], [], []
                for i, line in enumerate(self.kpoints.lines):
                    #print(line)
                    x.extend(line)
                    y.extend(branch[line])
                    s.extend(ders_onlines[i])
                    assert len(x) == len(y) == len(s)
                self.set_marker(asmarker, (x, y, s))

            return ders_onlines

        else:
            raise ValueError("Derivatives on homogeneous k-meshes are not supported yet")

    def effective_masses(self, spin, band, acc=4):
        """
        Compute the effective masses for the given `spin` and `band` index.
        Use finite difference with accuracy `acc`.

        Returns:
            numpy array of size nkpt.
        """
        ders2 = self.derivatives(spin, band, order=2, acc=acc) * (units.eV_to_Ha / units.bohr_to_ang**2)
        return 1. / ders2

    def effmass_line(self, spin, kpoint, band, acc=4):
        """
        Compute the effective masses along a line. Requires band energies on a k-path.

        Args:
            spin: Spin index.
            kpoint: integer or class:`Kpoint` object. Note that if kpoint is not an integer,
                and the path contains duplicated k-points, the first k-point is selected.
            band: Band index.
            acc: accuracy

        Returns:
        """
        if not self.kpoints.is_path:
            raise ValueError("effmass_line requires a k-path.")

        # Find index associate to the k-point
        ik = self.kindex(kpoint)

        # We have to understand if the k-point is a vertex or not.
        # If it's a vertex, indeed, we have to compute the left and right derivative
        # If kpt is inside the line, left and right derivatives are supposed to be equal
        for iline, line in enumerate(self.kpoints.lines):
            if line[-1] >= ik >= line[0]: break
        else:
            raise ValueError("Cannot find k-index %s in lines: %s" % (ik, self.kpoints.lines))

        kpos = line.index(ik)
        is_inside = kpos not in (0, len(line)-1)
        do_right = (not is_inside) and kpos != 0 and iline != len(self.kpoints.lines) - 1

        from abipy.tools.derivatives import finite_diff
        vals_on_line, h_left, vers_left = self._eigens_hvers_iline(spin, band, iline)
        d2line = finite_diff(vals_on_line, h_left, order=2, acc=acc) * (units.eV_to_Ha / units.bohr_to_ang**2)
        em_left = 1. / d2line[kpos]
        em_right = em_left
        h_right, vers_right = h_left, vers_left

        if do_right:
            kpos_right = self.kpoints.lines[iline+1].index(ik)
            assert kpos_right == 0
            vals_on_line, h_right, vers_right = self._eigens_hvers_iline(spin, band, iline+1)
            d2line = finite_diff(vals_on_line, h_right, order=2, acc=acc) * (units.eV_to_Ha / units.bohr_to_ang**2)
            em_right = 1. / d2line[kpos_right]

        return EffectiveMassAlongLine(spin, self.kpoints[ik], band, self.eigens[spin, ik, band],
                                      acc, self.structure.reciprocal_lattice,
                                      is_inside, h_left, vers_left, em_left, h_right, vers_right, em_right)

    def _eigens_hvers_iline(self, spin, band, iline):
        line = self.kpoints.lines[iline]
        vals_on_line = self.eigens[spin, line, band]
        h = self.kpoints.ds[line[0]]
        if not np.allclose(h, self.kpoints.ds[line[:-1]]):
            raise ValueError("For finite difference derivatives, the path must be homogeneous!\n" +
                             str(self.kpoints.ds[line[:-1]]))
        return vals_on_line, h, self.kpoints.versors[line[0]]

    def interpolate(self, lpratio=5, vertices_names=None, line_density=20,
                    kmesh=None, is_shift=None, filter_params=None, verbose=0):
        """
        Interpolate energies in k-space along a k-path and, optionally, in the IBZ for DOS calculations.
        Note that the interpolation will likely fail if there are symmetrical k-points in the input sampling
        so it's recommended to call this method with band structure obtained in the IBZ.

        Args:
            lpratio: Ratio between the number of star functions and the number of ab-initio k-points.
                The default should be OK in many systems, larger values may be required for accurate derivatives.
            vertices_names: Used to specify the k-path for the interpolated band structure
                It's a list of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path. `line_density` defines
                the density of the sampling. If None, the k-path is automatically generated according
                to the point group of the system.
            line_density: Number of points in the smallest segment of the k-path. Used with `vertices_names`.
            kmesh: Used to activate the interpolation on the homogeneous mesh for DOS (uses spglib API).
                kmesh is given by three integers and specifies mesh numbers along reciprocal primitive axis.
            is_shift: three integers (spglib API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.
            filter_params: TO BE described.
            verbose: Verbosity level

        Returns:
                namedtuple with the following attributes:

            ebands_kpath: :class:`ElectronBands` with the interpolated band structure on the k-path.
            ebands_kmesh: :class:`ElectronBands` with the interpolated band structure on the k-mesh.
                None if kmesh is not given.
            interpolator: :class:`SkwInterpolator` object.
        """
        # Get symmetries from abinit spacegroup (read from file).
        abispg = self.structure.abi_spacegroup
        if abispg is not None:
            fm_symrel = [s for (s, afm) in zip(abispg.symrel, abispg.symafm) if afm == 1]
        else:
            msg = ("Ebands object does not have symmetry operations `spacegroup.symrel`\n"
                   "This usually happens when ebands has not been initalized from a netcdf file\n."
                   "Will call spglib to get symmetry operations.")
            print(msg)
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            spglib_data = SpacegroupAnalyzer(self.structure).get_symmetry_dataset()
            fm_symrel = spglib_data["rotations"]

        # Build interpolator.
        from abipy.core.skw import SkwInterpolator
        my_kcoords = [k.frac_coords for k in self.kpoints]
        cell = (self.structure.lattice.matrix, self.structure.frac_coords,
                self.structure.atomic_numbers)

        skw = SkwInterpolator(lpratio, my_kcoords, self.eigens, self.fermie, self.nelect,
                              cell, fm_symrel, self.has_timrev,
                              filter_params=filter_params, verbose=verbose)

        # Generate k-points for interpolation.
        if vertices_names is None:
            vertices_names = [(k.frac_coords, k.name) for k in self.structure.hsym_kpoints]
        kpath = Kpath.from_vertices_and_names(self.structure, vertices_names, line_density=line_density)
        kfrac_coords, knames = kpath.frac_coords, kpath.names

        # Interpolate energies.
        eigens_kpath = skw.interp_kpts(kfrac_coords).eigens

        # Build new ebands object.
        kpts_kpath = Kpath(self.reciprocal_lattice, kfrac_coords, weights=None, names=knames)
        occfacts_kpath = np.zeros(eigens_kpath.shape)
        ebands_kpath = self.__class__(self.structure, kpts_kpath, eigens_kpath, self.fermie, occfacts_kpath,
                                      self.nelect, self.nspinor, self.nspden)

        ebands_kmesh = None
        if kmesh is not None:
            # Get kpts and weights in IBZ.
            kdos = Ktables(self.structure, kmesh, is_shift, self.has_timrev)
            eigens_kmesh = skw.interp_kpts(kdos.ibz).eigens

            # Build new ebands object with k-mesh
            #kptopt = kptopt_from_timrev()
            ksampling = KSamplingInfo.from_mpdivs(mpdivs=kmesh, shifts=[0,0,0], kptopt=1)
            kpts_kmesh = IrredZone(self.structure.reciprocal_lattice, kdos.ibz, weights=kdos.weights,
                                   names=None, ksampling=ksampling)
            occfacts_kmesh = np.zeros(eigens_kmesh.shape)
            ebands_kmesh = self.__class__(self.structure, kpts_kmesh, eigens_kmesh, self.fermie, occfacts_kmesh,
                                          self.nelect, self.nspinor, self.nspden)

        return dict2namedtuple(ebands_kpath=ebands_kpath, ebands_kmesh=ebands_kmesh, interpolator=skw)


class EffectiveMassAlongLine(object):
    """
    Store the value of the effective mass computed along a line.
    """
    def __init__(self, spin, kpoint, band, eig, acc, lattice,
                 is_inside, h_left, vers_left, em_left, h_right, vers_right, em_right):
        self.spin, self.kpoint, self.eig, self.band, self.acc, self.lattice = spin, kpoint, band, eig, acc, lattice,
        self.is_inside, self.h_left, self.vers_left, self.em_left, self.h_right, self.vers_right, self.em_right = \
            is_inside, h_left, vers_left, em_left, h_right, vers_right, em_right

    def __repr__(self):
        return "em_left: %s, em_right: %s" % (self.em_left, self.em_right)

    def __str__(self):
        lines = []; app = lines.append
        app("Effective masses for spin: %s, band: %s, accuracy: %s" % (self.spin, self.band, self.acc))
        app("K-point: %s, eigenvalue: %s [eV]" % (self.kpoint, self.eig))
        app("h_left: %s, h_right %s" % (self.h_left, self.h_right))
        app("is_inside: %s, vers_left: %s, vers_right: %s" % (self.is_inside, self.vers_left, self.vers_right))
        app("em_left: %s, em_right: %s" % (self.em_left, self.em_right))
        return "\n".join(lines)


def frame_from_ebands(ebands_objects, index=None, with_spglib=True):
    """
    Build a pandas dataframe with the most important results available in a list of band structures.

    Args:
        ebands_objects: List of objects that can be converted to structure.
            Support netcdf filenames or :class:`ElectronBands` objects
            See `ElectronBands.as_ebands` for the complete list.
        index: Index of the dataframe.
        with_spglib: If True, spglib is invoked to get the spacegroup symbol and number.

    Return:
        pandas :class:`DataFrame`
    """
    ebands_list = [ElectronBands.as_ebands(obj) for obj in ebands_objects]
    # Use OrderedDict to have columns ordered nicely.
    odict_list = [(ebands.get_dict4frame(with_spglib=with_spglib)) for ebands in ebands_list]

    import pandas as pd
    return pd.DataFrame(odict_list, index=index,
                        columns=list(odict_list[0].keys()) if odict_list else None)


class ElectronBandsPlotter(NotebookWriter):
    """
    Class for plotting electronic band structure and DOSes.
    Supports plots on the same graph or separated plots.

    Usage example:

    .. code-block:: python

        plotter = ElectronBandsPlotter()
        plotter.add_ebands("foo-label", "foo_GSR.nc")
        plotter.add_ebands("bar-label", "bar_WFK.nc")
        fig = plotter.gridplot()

    Dictionary with the mapping label --> edos.
    """
    _LINE_COLORS = ["b", "r",]
    _LINE_STYLES = ["-",":","--","-.",]
    _LINE_WIDTHS = [2,]

    def __init__(self, key_ebands=None, key_edos=None, edos_kwargs=None):
        """
        Args:
            key_ebands: List of (label, ebands) tuples.
                ebands is any object that can be converted into :class:`ElectronBands` e.g. ncfile, path.
            key_edos: List of (label, edos) tuples.
                edos is any object that can be converted into :class:`ElectronDos`
        """
        if key_ebands is None: key_ebands = []
        key_ebands = [(k, ElectronBands.as_ebands(v)) for k, v in key_ebands]
        self.ebands_dict = OrderedDict(key_ebands)

        if key_edos is None: key_edos = []
        key_edos = [(k, ElectronDos.as_edos(v, edos_kwargs)) for k, v in key_edos]
        self.edoses_dict = OrderedDict(key_edos)
        if key_edos:
            if not key_ebands:
                raise ValueError("key_ebands must be specifed when key_dos is not None")
            if len(key_ebands) != len(key_edos):
                raise ValueError("key_ebands and key_edos must have the same number of elements.")

    def __repr__(self):
        """Invoked by repr"""
        lines = []
        app = lines.append
        for i, (label, ebands) in enumerate(self.ebands_dict.items()):
            app("[%d] %s --> %s" % (i, label, repr(ebands)))

        if self.edoses_dict:
            for i, (label, edos) in enumerate(self.edoses_dict.items()):
                app("[%d] %s --> %s" % (i, label, repr(edos)))

        return "\n".join(lines)

    def get_ebands_frame(self, with_spglib=True):
        """
        Build a pandas dataframe with the most important results available in the band structures."""
        return frame_from_ebands(list(self.ebands_dict.values()),
                                 index=list(self.ebands_dict.keys()), with_spglib=with_spglib)

    @property
    def ebands_list(self):
        """"List of `:class:ElectronBands`."""
        return list(self.ebands_dict.values())

    @property
    def edoses_list(self):
        """"List of :class:`ElectronDos`."""
        return list(self.edoses_dict.values())

    def iter_lineopt(self):
        """Generates style options for lines."""
        for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    @deprecated(message="add_ebands_from_file method of ElectronBandsPlotter has been replaced by add_ebands. It will be removed in 0.4")
    def add_ebands_from_file(self, filepath, label=None):
        """
        Adds a band structure for plotting. Reads data from a Netcdfile
        """
        if label is None: label = os.path.abspath(filepath)
        ebands = ElectronBands.as_ebands(filepath)
        self.add_ebands(label, ebands)

    def add_ebands(self, label, bands, dos=None, edos_kwargs=None):
        """
        Adds a band structure and optionally a dos to the plotter.

        Args:
            label: label for the bands. Must be unique.
            bands: :class:`ElectronBands` object.
            dos: :class:`ElectronDos` object.
            edos_kwargs: optional dictionary with the options passed to `get_edos` to compute the electron DOS.
                Used only if `dos` is not None and it not an ElectronDos instance.
        """
        if label in self.ebands_dict:
            raise ValueError("label %s is already in %s" % (label, list(self.ebands_dict.keys())))

        self.ebands_dict[label] = ElectronBands.as_ebands(bands)
        if dos is not None:
            self.edoses_dict[label] = ElectronDos.as_edos(dos, edos_kwargs)

    def bands_statdiff(self, ref=0):
        """
        Compare the reference bands with index ref with the other bands stored in the plotter.
        """
        for i, label in enumerate(self.ebands_dict.keys()):
            if i == ref:
                ref_label = label
                break
        else:
            raise ValueError("ref index %s is > number of bands" % ref)

        ref_bands = self.ebands_dict[ref_label]

        text = []
        for label, bands in self.ebands_dict.items():
            if label == ref_label: continue
            stat = ref_bands.statdiff(bands)
            text.append(str(stat))

        return "\n\n".join(text)

    @add_fig_kwargs
    def combiplot(self, e0="fermie", ylims=None, **kwargs):
        """
        Plot the band structure and the DOS on the same figure.
        Use `gridplot` to plot band structures on different figures.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (ebands.fermie)
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   This should be ok in semiconductors. In metals, however, a better value of the Fermi energy
                   can be obtained from the DOS provided that the k-sampling for the DOS is much denser than
                   the one used to compute the density. See `edos_fermie`.
                - `edos_fermie`: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                   Available only if plotter contains dos objects.
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
        Returns:
            matplotlib figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        if self.edoses_dict:
            # Build grid with two axes.
            gspec = GridSpec(1, 2, width_ratios=[2, 1])
            gspec.update(wspace=0.05)
            # bands and DOS will share the y-axis
            ax1 = plt.subplot(gspec[0])
            ax2 = plt.subplot(gspec[1], sharey=ax1)
            ax_list = [ax1, ax2]
            fig = plt.gcf()
        else:
            # One axis for bands only
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax_list = [ax1]

        for ax in ax_list:
            ax.grid(True)
            set_axlims(ax, ylims, "y")

        # Plot ebands.
        lines, legends = [], []
        my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        for (label, ebands), lineopt in zip(self.ebands_dict.items(), self.iter_lineopt()):
            i += 1
            my_kwargs.update(lineopt)
            opts_label[label] = my_kwargs.copy()

            # Get energy zero.
            if e0 == "edos_fermie":
                mye0 = self.edoses_dict[label].fermie
            else:
                mye0 = ebands.get_e0(e0)

            l = ebands.plot_ax(ax1, mye0, spin=None, band=None, **my_kwargs)
            lines.append(l[0])

            # Use relative paths if label is a file.
            if os.path.isfile(label):
                legends.append("%s" % os.path.relpath(label))
            else:
                legends.append("%s" % label)

            # Set ticks and labels, legends.
            if i == 0:
                ebands.decorate_ax(ax1)

        ax1.legend(lines, legends, loc='upper right', shadow=True)

        # Add DOSes
        if self.edoses_dict:
            ax = ax_list[1]
            for label, edos in self.edoses_dict.items():
                ebands = self.edoses_dict[label]
                mye0 = ebands.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
                edos.plot_ax(ax, mye0, exchange_xy=True, **opts_label[label])

        return fig

    @deprecated(message="plot method of ElectronBandsPlotter has been replaced by combiplot. It will be removed in 0.4")
    def plot(self, *args, **kwargs):
        if "align" in kwargs or "xlim" in kwargs or "ylim" in kwargs:
            raise ValueError("align|xlim|ylim options are not supported anymore.")
        return self.combiplot(*args, **kwargs)

    @add_fig_kwargs
    def gridplot(self, e0="fermie", with_dos=True, ylims=None, **kwargs):
        """
        Plot multiple electron bandstructures and optionally DOSes on a grid.

        Args:
            eb_objects: List of objects from which the band structures are extracted.
                Each item in eb_objects is either a string with the path of the netcdf file,
                or one of the abipy object with an `ebands` attribute or a :class:`ElectronBands` object.
            edos_objects:
                List of objects from which the electron DOSes are extracted.
                Accept filepaths or :class:`ElectronDos` objects. If edos_objects is not None,
                each subplot in the grid contains a band structure with DOS else a simple bandstructure plot.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues and the DOS to have zero energy at the Fermi energy.
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   This should be ok in semiconductors. In metals, however, a better value of the Fermi energy
                   can be obtained from the DOS provided that the k-sampling for the DOS is much denser than
                   the one used to compute the density. See `edos_fermie`.
                - `edos_fermie`: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                   Available only if edos_objects is not None
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            with_dos: True if DOS should be printed.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used

        Returns:
            matplotlib figure.
        """
        titles = list(self.ebands_dict.keys())
        ebands_list, edos_list = self.ebands_list, self.edoses_list

        import matplotlib.pyplot as plt
        nrows, ncols = 1, 1
        numeb = len(ebands_list)
        if numeb > 1:
            ncols = 2
            nrows = numeb // ncols + numeb % ncols

        if not edos_list or not with_dos:
            # Plot grid with bands only.
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, squeeze=False)
            axes = axes.ravel()
            # don't show the last ax if numeb is odd.
            if numeb % ncols != 0: axes[-1].axis("off")

            for i, (ebands, ax) in enumerate(zip(ebands_list, axes)):
                ebands.plot(ax=ax, e0=e0, show=False)
                set_axlims(ax, ylims, "y")
                if titles is not None: ax.set_title(titles[i])
                if i % ncols != 0:
                    ax.set_ylabel("")

        else:
            # Plot grid with bands + DOS. see http://matplotlib.org/users/gridspec.html
            from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
            fig = plt.figure()
            gspec = GridSpec(nrows, ncols)

            for i, (ebands, edos) in enumerate(zip(ebands_list, edos_list)):
                subgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=gspec[i], width_ratios=[2, 1], wspace=0.05)
                # Get axes and align bands and DOS.
                ax1 = plt.subplot(subgrid[0])
                ax2 = plt.subplot(subgrid[1], sharey=ax1)
                set_axlims(ax1, ylims, "y")
                set_axlims(ax2, ylims, "y")

                # Define the zero of energy and plot
                mye0 = ebands.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
                ebands.plot_with_edos(edos, e0=mye0, axlist=(ax1, ax2), show=False)

                if titles is not None: ax1.set_title(titles[i])
                if i % ncols != 0:
                    for ax in (ax1, ax2):
                        ax.set_ylabel("")

        return fig

    @add_fig_kwargs
    def boxplot(self, e0="fermie", brange=None, swarm=False, **kwargs):
        """
        Use seaborn to draw a box plot to show distributions of eigenvalues with respect to the band index.
        Band structures are drawn on different subplots.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            brange: Only bands such as `brange[0] <= band_index < brange[1]` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            kwargs: Keyword arguments passed to seaborn boxplot.
        """
        # Build grid of plots.
        num_plots, ncols, nrows = len(self.ebands_dict), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()
        # don't show the last ax if numeb is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        for (label, ebands), ax in zip(self.ebands_dict.items(), ax_list):
            ebands.boxplot(ax=ax, show=False)
            ax.set_title(label)
        return fig

    @add_fig_kwargs
    def combiboxplot(self, e0="fermie", brange=None, swarm=False, ax=None, **kwargs):
        """
        Use seaborn to draw a box plot comparing the distributions of the eigenvalues
        Band structures are drawn on the same plot.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0

            brange: Only bands such as `brange[0] <= band_index < brange[1]` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            kwargs: Keyword arguments passed to seaborn boxplot.
        """
        spin_polarized = False
        frames = []
        for label, ebands in self.ebands_dict.items():
            # Get the dataframe, select bands and add column with label
            frame = ebands.to_dataframe(e0=e0)
            if brange is not None: frame = frame[brange[0] <= frame["band"] < brange[1]]
            frame["label"] = label
            frames.append(frame)
            if ebands.nsppol == 2: spin_polarized = True

        # Merge frames ignoring index (not meaningful)
        import pandas as pd
        data = pd.concat(frames, ignore_index=True)

        import matplotlib.pyplot as plt
        import seaborn.apionly as sns
        if not spin_polarized:
            ax, fig, plt = get_ax_fig_plt(ax=ax)
            ax.grid(True)
            sns.boxplot(x="band", y="eig", data=data, hue="label", ax=ax, **kwargs)
            if swarm:
                sns.swarmplot(x="band", y="eig", data=data, hue="label", color=".25", ax=ax)
        else:
            if ax is not None:
                raise NotImplementedError("ax == None not implemented when nsppol==2")
            fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, squeeze=False)
            for spin, ax in zip(range(2), axes.ravel()):
                ax.grid(True)
                data_spin = data[data["spin"] == spin]
                sns.boxplot(x="band", y="eig", data=data_spin, hue="label", ax=ax, **kwargs)
                if swarm:
                    sns.swarmplot(x="band", y="eig", data=data_spin, hue="label", color=".25", ax=ax)

        return fig

    def animate(self, e0="fermie", interval=250, savefile=None, show=True):
        """
        Use matplotlib to animate a list of band structure plots (with or without DOS).

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues and the DOS to have zero energy at the Fermi energy.
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   See `edos_fermie`.
                - `edos_fermie`: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            interval: draws a new frame every interval milliseconds.
            savefile: Use e.g. 'myanimation.mp4' to save the animation in mp4 format.
            show: True if the animation should be shown immediately

        Returns:
            Animation object.

        See also http://matplotlib.org/api/animation_api.html
                 http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

        Note:
            It would be nice to have the possibility of animating the title of the plot, unfortunately
            this feature is not available in the present version of matplotlib. See
            http://stackoverflow.com/questions/17558096/animated-title-in-matplotlib
        """
        ebands_list, edos_list = self.ebands_list, self.edoses_list
        if edos_list and len(edos_list) != len(ebands_list):
            raise ValueError("The number of objects for DOS must be equal to the number of bands")

        import matplotlib.pyplot as plt
        fig = plt.figure()
        plotax_kwargs = {"color": "black", "linewidth": 2.0}

        artists = []
        if not edos_list:
            # Animation with band structures
            ax = fig.add_subplot(1, 1, 1)
            ebands_list[0].decorate_ax(ax)
            for i, ebands in enumerate(ebands_list):
                lines = ebands.plot_ax(ax, e0, **plotax_kwargs)
                #if titles is not None:
                #    lines += [ax.set_title(titles[i])]
                artists.append(lines)
        else:
            # Animation with band structures + DOS.
            from matplotlib.gridspec import GridSpec
            gspec = GridSpec(1, 2, width_ratios=[2, 1])
            gspec.update(wspace=0.05)
            ax1 = plt.subplot(gspec[0])
            ax2 = plt.subplot(gspec[1], sharey=ax1)
            ebands_list[0].decorate_ax(ax1)
            ax2.grid(True)
            ax2.yaxis.set_ticks_position("right")
            ax2.yaxis.set_label_position("right")

            for i, (ebands, edos) in enumerate(zip(ebands_list, edos_list)):
                # Define the zero of energy to align bands and dos
                mye0 = ebands.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
                ebands_lines = ebands.plot_ax(ax1, mye0, **plotax_kwargs)
                edos_lines = edos.plot_ax(ax2, mye0, exchange_xy=True, **plotax_kwargs)
                lines = ebands_lines + edos_lines
                #if titles is not None:
                #    lines += [ax.set_title(titles[i])]
                artists.append(lines)

        import matplotlib.animation as animation
        anim = animation.ArtistAnimation(fig, artists, interval=interval,
                                         blit=False, # True is faster but then the movie starts with an empty frame!
                                         #repeat_delay=1000
                                         )

        if savefile is not None: anim.save(savefile)
        if show: plt.show()
        return anim

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        # Use pickle files for data persistence.
        tmpfile = self.pickle_dump()

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("plotter = abilab.ElectronBandsPlotter.pickle_load('%s')" % tmpfile),
            nbv.new_code_cell("print(plotter)"),
            nbv.new_code_cell("frame = plotter.get_ebands_frame()\ndisplay(frame)"),
            nbv.new_code_cell("ylims = (None, None)"),
            nbv.new_code_cell("fig = plotter.gridplot(ylims=ylims)"),
            nbv.new_code_cell("fig = plotter.combiplot(ylims=ylims)"),
            nbv.new_code_cell("fig = plotter.boxplot()"),
            nbv.new_code_cell("fig = plotter.combiboxplot()"),
            nbv.new_code_cell("if False: anim = plotter.animate()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

    #def _can_use_basenames_as_labels(self):
    #    """
    #    Return True if all labels represent valid files and the basenames are unique
    #    In this case one can use the file basename instead of the full path in the plots.
    #    """
    #    if not all(os.path.exists(l) for l in self.ebands_dict): return False
    #    labels = [os.path.basename(l) for l in self.ebands_dict]
    #    return len(set(labels)) == len(labels)


class ElectronsReader(ETSF_Reader, KpointsReaderMixin):
    """
    This object reads band structure data from a netcdf file.
    """
    def read_ebands(self):
        """
        Returns an instance of :class:`ElectronBands`. Main entry point for client code
        """
        return ElectronBands(
            structure=self.read_structure(),
            kpoints=self.read_kpoints(),
            eigens=self.read_eigenvalues(),
            fermie=self.read_fermie(),
            occfacts=self.read_occupations(),
            nelect=self.read_nelect(),
            nspinor=self.read_nspinor(),
            nspden=self.read_nspden(),
            nband_sk=self.read_nband_sk(),
            smearing=self.read_smearing(),
            )

    def read_nband_sk(self):
        """Array with the number of bands indexed by [s,k]."""
        return self.read_value("number_of_states")

    def read_nspinor(self):
        """Number of spinors."""
        return self.read_dimvalue("number_of_spinor_components")

    def read_nspden(self):
        """Number of spin-density components"""
        # FIXME: default 1 is needed for SIGRES files (abinit8)
        return self.read_dimvalue("number_of_components", default=1)

    def read_tsmear(self):
        return self.read_value("smearing_width")

    def read_eigenvalues(self):
        """Eigenvalues in eV."""
        return units.ArrayWithUnit(self.read_value("eigenvalues"), "Ha").to("eV")

    def read_occupations(self):
        """Occupancies."""
        return self.read_value("occupations")

    def read_fermie(self):
        """Fermi level in eV."""
        return units.Energy(self.read_value("fermi_energy"), "Ha").to("eV")

    def read_nelect(self):
        """Number of valence electrons."""
        return self.read_value("number_of_electrons")

    def read_smearing(self):
        """Returns a :class:`Smearing` instance with info on the smearing technique."""
        occopt = int(self.read_value("occopt"))

        try:
            scheme = "".join(c for c in self.read_value("smearing_scheme"))
            scheme = scheme.strip()
        except TypeError as exc:
            scheme = None

        return Smearing(
            scheme=scheme,
            occopt=occopt,
            tsmear_ev=units.Energy(self.read_value("smearing_width"), "Ha").to("eV")
        )


class ElectronDos(object):
    """
    This object stores the electronic density of states.
    It is usually created by calling the get_edos method of :class:`ElectronBands`.
    """

    def __init__(self, mesh, spin_dos, nelect, fermie=None):
        """
        Args:
            mesh: array-like object with the mesh points.
            spin_dos: array-like object with the DOS value for the different spins.
                      spin_dos[1, nw] if spin-unpolarized.
                      spin_dos[2, nw] if spin-polarized case.
            nelect: Number of electrons in the unit cell.
            fermie: Fermi level in eV. If None, fermie is obtained from the idos.

        .. note::

            mesh is given in eV, spin_dos is in states/eV.
        """
        spin_dos = np.atleast_2d(spin_dos)
        self.nsppol = len(spin_dos)
        self.nelect = nelect

        # Save DOS and IDOS for each spin.
        sumv = np.zeros(len(mesh))
        self.spin_dos, self.spin_idos = [], []
        for values in spin_dos:
            sumv += values
            f = Function1D(mesh, values)
            self.spin_dos.append(f)
            self.spin_idos.append(f.integral())

        # Total DOS and IDOS.
        if self.nsppol == 1: sumv = 2 * sumv
        self.tot_dos = Function1D(mesh, sumv)
        self.tot_idos = self.tot_dos.integral()

        if fermie is not None:
            self.fermie = float(fermie)
        else:
            # *Compute* fermie from nelect. Note that this value could differ
            # from the one stored in ElectronBands (coming from the SCF run)
            # The accuracy of self.fermie depends on the number of k-points used for the DOS
            # and the parameters used to call ebands.get_edos.
            try:
                self.fermie = self.find_mu(self.nelect)
            except ValueError:
                print("tot_idos values:")
                print(self.tot_idos)
                raise

    def __str__(self):
        lines = []; app = lines.append
        app("nsppol=%d, nelect=%s" % (self.nsppol, self.nelect))
        app("Fermi energy: %s (recomputed from nelect):" % self.fermie)
        return "\n".join(lines)

    @classmethod
    def as_edos(cls, obj, edos_kwargs):
        """
        Return an instance of :class:`ElectronDos` from a generic obj.
        Supports:

            - instances of cls
            - files (string) that can be open with abiopen and that provide an `ebands` attribute.
            - objects providing an `ebands` or `get_edos` attribute

        Args:
            edos_kwargs: optional dictionary with the options passed to `get_edos` to compute the electron DOS.
            Used when obj is not already an instance of `cls`.
        """
        if edos_kwargs is None: edos_kwargs = {}
        if isinstance(obj, cls):
            return obj
        elif is_string(obj):
            # path?
            if obj.endswith(".pickle"):
                with open(obj, "rb") as fh:
                    return cls.as_edos(pickle.load(fh), edos_kwargs)

            from abipy.abilab import abiopen
            with abiopen(obj) as abifile:
                return abifile.ebands.get_edos(**edos_kwargs)
        elif hasattr(obj, "ebands"):
            return obj.ebands.get_edos(**edos_kwargs)
        elif hasattr(obj, "get_edos"):
            return obj.get_edos(**edos_kwargs)

        raise TypeError("Don't know how to create `ElectronDos` from %s" % type(obj))

    def __eq__(self, other):
        if other is None: return False
        if self.nsppol != other.nsppol: return False
        for f1, f2 in zip(self.spin_dos, other.spin_dos):
            if f1 != f2: return False
        return True

    def __ne__(self, other):
        return not (self == other)

    def dos_idos(self, spin=None):
        """
        Returns DOS and IDOS for given spin. Total DOS and IDOS if spin is None.
        """
        if spin is None:
            return self.tot_dos, self.tot_idos
        else:
            return self.spin_dos[spin], self.spin_idos[spin]

    def find_mu(self, nelect, spin=None):
        """
        Finds the chemical potential given the number of electrons.
        """
        idos = self.tot_idos if spin is None else self.spin_idos[spin]

        # Cannot use bisection because DOS might be negative due to smearing.
        # This one is safer albeit slower.
        for i, (ene, intg) in enumerate(idos):
            if intg > nelect: break
        else:
            raise ValueError("Cannot find I(e) such that I(e) > nelect")

        # Use linear interpolation to find mu (useful if mesh is coarse)
        e0, y0 = idos[i-1]
        e1, y1 = idos[i]

        alpha = (y1 - y0) / (e1 - e0)
        beta = y0 - alpha * e0
        mu = (nelect - beta) / alpha
        #print("idos[i-1]:", idos[i-1], "idos[i]:", idos[i], "intg", intg, "nelect", nelect)
        #print("mu linear", mu)
        return mu

    def get_e0(self, e0):
        """
        e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
        """
        if e0 is None:
            return 0.0
        elif is_string(e0):
            if e0 == "fermie":
                return self.fermie
            elif e0 == "None":
                return 0.0
            else:
                raise ValueError("Wrong value for e0: %s" % e0)
        else:
            # Assume number
            return float(e0)

    def plot_ax(self, ax, e0, spin=None, what="d", fact=1.0, exchange_xy=False, **kwargs):
        """
        Helper function to plot the DOS data on the axis ax.

        Args:
            ax: matplotlib axis.
            e0: Option used to define the zero of energy in the band structure plot.
            spin: selects the spin component, None for total DOS, IDOS.
            what: string selecting what will be plotted:
                  "d" for DOS, "i" for IDOS. chars can be concatenated
                  hence what="id" plots both IDOS and DOS. (default "d").
            fact: Multiplication factor for DOS/IDOS. Usually +-1 for spin DOS
            exchange_xy: True to exchange x-y axes.
            kwargs:
                Options passes to matplotlib.

        Return:
            list of lines added to the axis.
        """
        dosf, idosf = self.dos_idos(spin=spin)
        opts = [c.lower() for c in what]

        e0 = self.get_e0(e0)
        lines = []
        for c in opts:
            if c == "d": f = dosf
            if c == "i": f = idosf
            xx, yy = f.mesh - e0, f.values * fact
            if exchange_xy: xx, yy = yy, xx
            lines.extend(ax.plot(xx, yy, **kwargs))
        return lines

    @add_fig_kwargs
    def plot(self, e0="fermie", spin=None, ax=None, xlims=None, **kwargs):
        """
        Plot DOS

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            spin: Selects the spin component, None if total DOS is wanted.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            kwargs: options passed to ax.plot.

        Returns:
            matplotlib figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)
        ax.set_xlabel('Energy [eV]')
        set_axlims(ax, xlims, "x")
        e0 = self.get_e0(e0)
        if not kwargs:
            kwargs = {"color": "black", "linewidth": 1.0}

        for spin in range(self.nsppol):
            spin_sign = +1 if spin == 0 else -1
            x, y = self.spin_dos[spin].mesh - e0, spin_sign * self.spin_dos[spin].values
            ax.plot(x, y, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_dos_idos(self, e0="fermie", spin=None, xlims=None, **kwargs):
        """
        Plot DOS and IDOS.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            spin: Selects the spin component, None if total DOS is wanted.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            kwargs: options passed to plot_ax

        Returns:
            matplotlib figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        gspec = GridSpec(2, 1, height_ratios=[1, 2])
        gspec.update(wspace=0.05)
        ax1 = plt.subplot(gspec[0])
        ax2 = plt.subplot(gspec[1])

        for ax in (ax1, ax2):
            ax.grid(True)
            set_axlims(ax, xlims, "x")

        ax1.set_ylabel("TOT IDOS" if spin is None else "IDOS (spin %s)" % spin)
        ax2.set_ylabel("TOT DOS" if spin is None else "DOS (spin %s)" % spin)
        ax2.set_xlabel('Energy [eV]')

        self.plot_ax(ax1, e0, spin=spin, what="i", **kwargs)
        self.plot_ax(ax2, e0, spin=spin, what="d", **kwargs)

        fig = plt.gcf()
        return fig

    @add_fig_kwargs
    def plot_up_minus_down(self, e0="fermie", ax=None, xlims=None, **kwargs):
        """
        Plot DOS

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            kwargs: options passed to ax.plot.

        Returns:
            matplotlib figure.
        """
        if self.nsppol == 1: # DOH!
            dos_diff = Function1D.from_constant(self.spin_dos[0].mesh, 0.0)
        else:
            dos_diff = self.spin_dos[0] - self.spin_dos[1]
        idos_diff= dos_diff.integral()

        e0 = self.get_e0(e0)
        if not kwargs:
            kwargs = {"color": "black", "linewidth": 1.0}

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.plot(dos_diff.mesh - e0, dos_diff.values, **kwargs)
        ax.plot(idos_diff.mesh - e0, idos_diff.values, **kwargs)

        ax.grid(True)
        set_axlims(ax, xlims, "x")
        ax.set_xlabel('Energy [eV]')

        return fig


class ElectronDosPlotter(NotebookWriter):
    """
    Class for plotting electronic DOSes.

    Usage example:

    .. code-block:: python

        plotter = ElectronDosPlotter()
        plotter.add_edos("foo dos", "foo.nc")
        plotter.add_edos("bar dos", "bar.nc")
        fig = plotter.gridplot()
    """
    #_LINE_COLORS = ["b", "r",]
    #_LINE_STYLES = ["-",":","--","-.",]
    #_LINE_WIDTHS = [2,]

    def __init__(self, key_edos=None, edos_kwargs=None):
        if key_edos is None: key_edos = []
        key_edos = [(k, ElectronDos.as_edos(v, edos_kwargs)) for k, v in key_edos]
        self.edoses_dict = OrderedDict(key_edos)

    #def iter_lineopt(self):
    #    """Generates style options for lines."""
    #    for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
    #        yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    @property
    def edos_list(self):
        """List of DOSes"""
        return list(self.edoses_dict.values())

    @deprecated(message="add_edos_from_file method of ElectronDosPlotter has been replaced by add_edos. It will be removed in 0.4")
    def add_edos_from_file(self, filepath, label=None, method="gaussian", step=0.1, width=0.2):
        """
        Adds a dos for plotting. Reads data from a Netcdf file
        """
        ebands = ElectronBands.as_ebands(filepath)
        edos = ebands.get_edos(method=method, step=step, width=width)
        if label is None: label = filepath
        self.add_edos(label, edos)

    def add_edos(self, label, edos, edos_kwargs=None):
        """
        Adds a DOS for plotting.

        Args:
            label: label for the DOS. Must be unique.
            dos: :class:`ElectronDos` object.
            edos_kwargs: optional dictionary with the options passed to `get_edos` to compute the electron DOS.
                Used only if `edos` is not an ElectronDos instance.
        """
        if label in self.edoses_dict:
            raise ValueError("label %s is already in %s" % (label, list(self.edoses_dict.keys())))
        self.edoses_dict[label] = ElectronDos.as_edos(edos, edos_kwargs)

    @add_fig_kwargs
    def combiplot(self, ax=None, e0="fermie", xlims=None, **kwargs):
        """
        Plot the the DOSes on the same figure.
        Use `gridplot` to plot DOSes on different figures.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used

        Returns:
            `matplotlib` figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        can_use_basename = self._can_use_basenames_as_labels()

        for label, dos in self.edoses_dict.items():
            if can_use_basename:
                label = os.path.basename(label)
            else:
                # Use relative paths if label is a file.
                if os.path.isfile(label): label = os.path.relpath(label)
            dos.plot_ax(ax, e0, label=label)

        ax.grid(True)
        ax.set_xlabel("Energy [eV]")
        ax.set_ylabel('DOS [states/eV]')
        set_axlims(ax, xlims, "x")
        ax.legend(loc="best")

        return fig

    @deprecated(message="plot method of ElectronDos has been replaced by combiplot. It will be removed in 0.4")
    def plot(self, *args, **kwargs):
        return self.combiplot(*args, **kwargs)

    @add_fig_kwargs
    def gridplot(self, e0="fermie", xlims=None, **kwargs):
        """
        Plot multiple DOSes on a grid.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues and the DOS to have zero energy at the Fermi energy.
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   This should be ok in semiconductors. In metals, however, a better value of the Fermi energy
                   can be obtained from the DOS provided that the k-sampling for the DOS is much denser than
                   the one used to compute the density. See `edos_fermie`.
                - `edos_fermie`: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                   Available only if edos_objects is not None
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used

        Returns:
            matplotlib figure.
        """
        titles = list(self.edoses_dict.keys())
        edos_list = self.edos_list

        import matplotlib.pyplot as plt
        nrows, ncols = 1, 1
        numeb = len(edos_list)
        if numeb > 1:
            ncols = 2
            nrows = numeb // ncols + numeb % ncols

        # Build Grid
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, squeeze=False)
        axes = axes.ravel()
        # don't show the last ax if numeb is odd.
        if numeb % ncols != 0: axes[-1].axis("off")

        for i, (label, edos) in enumerate(self.edoses_dict.items()):
            ax = axes[i]
            edos.plot(ax=ax, e0=e0, show=False)
            ax.set_title(label)
            set_axlims(ax, xlims, "x")
            if i % ncols != 0:
                ax.set_ylabel("")

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        # Use pickle files for data persistence.
        tmpfile = self.pickle_dump()

        nb.cells.extend([
            nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("plotter = abilab.ElectronDosPlotter.pickle_load('%s')" % tmpfile),
            nbv.new_code_cell("print(plotter)"),
            nbv.new_code_cell("xlims = (None, None)"),
            nbv.new_code_cell("fig = plotter.combiplot(xlims=xlims)"),
            nbv.new_code_cell("fig = plotter.gridplot(xlims=xlims)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

    #def animate(self, **kwargs):
    #    animator = Animator()
    #    import tempfile
    #    tmpdir = tempfile.mkdtemp()
    #    for (label, dos) in self.edoses_dict.items():
    #        savefig = os.path.join(tmpdir, label + ".png")
    #        dos.plot(show=False, savefig=savefig)
    #        animator.add_figure(label, savefig)
    #    return animator.animate(**kwargs)

    def _can_use_basenames_as_labels(self):
        """
        Return True if all labels represent valid files and the basenames are unique
        In this case one can use the file basename instead of the full path in the plots.
        """
        if not all(os.path.exists(l) for l in self.edoses_dict): return False
        labels = [os.path.basename(l) for l in self.edoses_dict]
        return len(set(labels)) == len(labels)
