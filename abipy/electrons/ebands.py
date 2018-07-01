# coding: utf-8
"""Classes to analyse electronic structures."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import copy
import itertools
import json
import warnings
import tempfile
import pickle
import numpy as np
import pymatgen.core.units as units

from collections import OrderedDict, namedtuple, Iterable
from monty.string import is_string, list_strings, marquee
from monty.termcolor import cprint
from monty.json import MSONable, MontyEncoder
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
from monty.bisect import find_le, find_gt
try:
    from pymatgen.util.serialization import pmg_serialize
except ImportError:
    from pymatgen.serializers.json_coders import pmg_serialize
from pymatgen.electronic_structure.core import Spin as PmgSpin
from abipy.core.func1d import Function1D
from abipy.core.mixins import Has_Structure, NotebookWriter
from abipy.core.kpoints import (Kpoint, KpointList, Kpath, IrredZone, KSamplingInfo, KpointsReaderMixin,
    Ktables, has_timrev_from_kptopt, map_grid2ibz, kmesh_from_mpdivs)
from abipy.core.structure import Structure
from abipy.iotools import ETSF_Reader
from abipy.tools import gaussian, duck
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
    get_ax3d_fig_plt, rotate_ticklabels, set_visible, plot_unit_cell)

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "ElectronBands",
    "ElectronDos",
    "dataframe_from_ebands",
    "ElectronBandsPlotter",
    "ElectronDosPlotter",
]


class Electron(namedtuple("Electron", "spin kpoint band eig occ kidx")):
    """
    Single-particle state.

    .. Attributes:

        spin: spin index (C convention, i.e >= 0)
        kpoint: |Kpoint| object.
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
    #            self.band == other.band and  # ??
    #            self.eig == other.eig
                 # and self.occ == other.occ
    #            )

    #def __ne__(self, other):
    #    return not (self == other)

    def __str__(self):
        return "spin=%d, kpt=%s, band=%d, eig=%.3f, occ=%.3f" % (
            self.spin, self.kpoint, self.band, self.eig, self.occ)

    @property
    def skb(self):
        """Tuple with (spin, kpoint, band)."""
        return self.spin, self.kpoint, self.band

    def copy(self):
        """Shallow copy."""
        return self.__class__(**{f: copy.copy(getattr(self, f)) for f in self._fields})

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
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app("Energy: %.3f (eV)" % self.energy)
        app("Initial state: %s" % str(self.in_state))
        app("Final state:   %s" % str(self.out_state))

        return "\n".join(lines)

    #def __eq__(self, other):
    #    if other is None: return False
    #    if not isinstance(other, self.__class__): return False
    #    return self.in_state == other.in_state and
    #           self.out_state == other.out_state

    #def __ne__(self, other):
    #    return not (self == other)

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
        Returns a JSON_ string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    @classmethod
    def as_smearing(cls, obj):
        """"
        Convert obj into a Smearing instance.
        Accepts: Smearing instance, None (if info are not available), Dict-like object.
        """
        if isinstance(obj, cls): return obj
        if obj is None:
            return cls(scheme=None, occopt=1, tsmear_ev=0.0)

        # Assume dict-like object.
        try:
            return cls(**obj)
        except Exception as exc:
            raise TypeError("Don't know how to convert %s into Smearing object:\n%s" % (type(obj), str(exc)))

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
        return "mean = %.3f, stdev = %.3f, min = %.3f, max = %.3f (eV)" % (
            self.mean, self.stdev, self.min, self.max)


class ElectronBandsError(Exception):
    """Exceptions raised by ElectronBands."""


class ElectronBands(Has_Structure):
    """
    This object stores the electronic band structure.

    .. attribute:: fermie

            Fermi level in eV. Note that, if the band structure has been computed
            with a NSCF run, fermie corresponds to the fermi level obtained
            in the SCF run that produced the density used for the band structure calculation.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ElectronBands
    """
    Error = ElectronBandsError

    # FIXME
    # Increase a bit the value of fermie used in bisection routines to solve the problem mentioned below
    pad_fermie = 1e-3
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
        Initialize an instance of |ElectronBands| from the netCDF file ``filepath``.
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
        """Reconstruct object from dictionary ``d``."""
        d = d.copy()
        kd = d["kpoints"].copy()
        kd.pop("@module")

        kpoints_cls = KpointList.subclass_from_name(kd.pop("@class"))
        kpoints = kpoints_cls.from_dict(kd)

        # Needed to support old dictionaries
        if "nspden" not in d: d["nspden"] = 1
        if "nspinor" not in d: d["nspinor"] = 1
        return cls(Structure.from_dict(d["structure"]), kpoints,
                   d["eigens"], d["fermie"], d["occfacts"], d["nelect"], d["nspinor"], d["nspden"],
                   nband_sk=d["nband_sk"], smearing=d["smearing"],
                   linewidths=d.get("linewidths", None)
                   )

    @pmg_serialize
    def as_dict(self):
        """Return dictionary with JSON_ serialization."""
        linewidths = None if not self.has_linewidths else self.linewidths.tolist()
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
            linewidths=linewidths,
        )

    @classmethod
    def as_ebands(cls, obj):
        """
        Return an instance of |ElectronBands| from a generic object `obj`.
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

    @classmethod
    def from_mpid(cls, material_id, api_key=None, endpoint=None,
                  nelect=None, has_timerev=True, nspinor=1, nspden=None):
        """
        Read bandstructure data corresponding to a materials project ``material_id``.
        and return Abipy ElectronBands object.

        Args:
            material_id (str): Materials Project material_id (a string, e.g., mp-1234).
            api_key (str): A String API key for accessing the MaterialsProject
                REST interface. Please apply on the Materials Project website for one.
                If this is None, the code will check if there is a `PMG_MAPI_KEY` in
                your .pmgrc.yaml. If so, it will use that environment
                This makes easier for heavy users to simply add
                this environment variable to their setups and MPRester can
                then be called without any arguments.
            endpoint (str): Url of endpoint to access the MaterialsProject REST interface.
                Defaults to the standard Materials Project REST address, but
                can be changed to other urls implementing a similar interface.
            nelect: Number of electrons in the unit cell.
            nspinor: Number of spinor components.
        """
        # Get pytmatgen structure and convert it to abipy structure
        from abipy.core import restapi
        with restapi.get_mprester(api_key=api_key, endpoint=endpoint) as rest:
            pmgb = rest.get_bandstructure_by_material_id(material_id=material_id)

            # Structure is set to None so we have to perform another request and patch the object.
            structure = rest.get_structure_by_material_id(material_id, final=True)
            if pmgb.structure is None: pmgb.structure = structure
            #pmgb = pmgb.__class__.from_dict(pmgb.as_dict())

        if nelect is None:
            # Get nelect from valence band maximum index.
            if pmgb.is_metal():
                raise RuntimeError("Nelect must be specified if metallic bands.")
            else:
                d = pmgb.get_vbm()
                iv_up = max(d["band_index"][PmgSpin.up])
                nelect = (iv_up + 1) * 2
                #print("iv_up", iv_up, "nelect: ", nelect)
                if pmgb.is_spin_polarized:
                    iv_down = max(d["band_index"][PmgSpin.down])
                    assert iv_down == iv_up

        #ksampling = KSamplingInfo.from_kbounds(kbounds)

        return cls.from_pymatgen(pmgb, nelect, weights=None, has_timerev=has_timerev,
                                 ksampling=None, smearing=None, nspinor=nspinor, nspden=nspden)

    def to_json(self):
        """
        Returns a JSON_ string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    def __init__(self, structure, kpoints, eigens, fermie, occfacts, nelect, nspinor, nspden,
                 nband_sk=None, smearing=None, linewidths=None):
        """
        Args:
            structure: |Structure| object.
            kpoints: |KpointList| instance.
            eigens: Array-like object with the eigenvalues (eV) stored as [s, k, b]
                where s: spin , k: kpoint, b: band index
            fermie: Fermi level in eV.
            occfacts: Occupation factors (same shape as eigens)
            nelect: Number of valence electrons in the unit cell.
            nspinor: Number of spinor components
            nspden: Number of independent density components.
            nband_sk: Array-like object with the number of bands treated at each [spin,kpoint]
                      If not given, nband_sk is initialized from eigens.
            smearing: :class:`Smearing` object storing information on the smearing technique.
            linewidths: Array-like object with the linewidths (eV) stored as [s, k, b]
        """
        self._structure = structure

        # Eigenvalues and occupancies are stored in ndarrays ordered by [spin,kpt,band]
        self._eigens = np.atleast_3d(eigens)
        self._occfacts = np.atleast_3d(occfacts)
        assert self._eigens.shape == self._occfacts.shape
        self._linewidths = None
        if linewidths is not None:
            self._linewidths = np.reshape(linewidths, self._eigens.shape)

        self.nsppol, self.nkpt, self.mband = self.eigens.shape
        self.nspinor, self.nspden = nspinor, nspden

        if nband_sk is not None:
            self.nband_sk = np.array(nband_sk)
        else:
            self.nband_sk = np.array(self.nsppol * self.nkpt * [self.mband])
            self.nband_sk.shape = (self.nsppol, self.nkpt)

        self.kpoints = kpoints
        assert self.nkpt == len(self.kpoints)
        assert isinstance(self.kpoints, KpointList)

        self.smearing = {} if smearing is None else smearing
        self.nelect = float(nelect)
        self.fermie = float(fermie)

        # Recompute the Fermi level (in principle should do this only if
        # bands are computed on a BZ mesh with a NSCF run.
        #if self.kpoints.is_ibz: # and iscf < 0
        #    self.recalc_fermie()

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    @lazy_property
    def _auto_klabels(self):
        # Find the k-point names in the pymatgen database.
        # We'll use _auto_klabels to label the point in the matplotlib plot
        # if klabels are not specified by the user.

        _auto_klabels = OrderedDict()
        # If the first or the last k-point are not recognized in findname_in_hsym_stars
        # matplotlib won't show the full band structure along the k-path
        # because the labels are not defined. So we have to make sure that
        # the labels for the extrema of the path are always defined.
        _auto_klabels[0] = " "

        for idx, kpoint in enumerate(self.kpoints):
            name = kpoint.name if kpoint.name is not None else self.structure.findname_in_hsym_stars(kpoint)
            if name is not None:
                _auto_klabels[idx] = name
                if kpoint.name is None: kpoint.set_name(name)

        last = len(self.kpoints) - 1
        if last not in _auto_klabels: _auto_klabels[last] = " "

        return _auto_klabels

    def __repr__(self):
        """String representation (short version)"""
        return "<%s, nk=%d, %s, id=%s>" % (self.__class__.__name__, self.nkpt, self.structure.formula, id(self))

    def __str__(self):
        """String representation"""
        return self.to_string()

    def __add__(self, other):
        """self + other returns a |ElectronBandsPlotter|."""
        if not isinstance(other, (ElectronBands, ElectronBandsPlotter)):
            raise TypeError("Cannot add %s to %s" % (type(self), type(other)))

        if isinstance(other, ElectronBandsPlotter):
            self_key = repr(self)
            other.add_ebands(self_key, self)
            return other
        else:
            plotter = ElectronBandsPlotter()
            self_key = repr(self)
            plotter.add_ebands(self_key, self)
            self_key = repr(self)
            other_key = repr(other)
            plotter.add_ebands(other_key, other)
            return plotter

    __radd__ = __add__

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
        """Eigenvalues in eV. |numpy-array| with shape [nspin, nkpt, mband]."""
        return self._eigens

    @property
    def linewidths(self):
        """linewidths in eV. |numpy-array| with shape [nspin, nkpt, mband]."""
        return self._linewidths

    @linewidths.setter
    def linewidths(self, linewidths):
        """Set the linewidths. Accept real array of shape [nspin, nkpt, mband] or None."""
        if linewidths is not None:
            linewidths = np.reshape(linewidths, self.shape)
        self._linewidths = linewidths

    @property
    def has_linewidths(self):
        """True if bands with linewidths."""
        return getattr(self, "_linewidths", None) is not None

    @property
    def occfacts(self):
        """Occupation factors. |numpy-array| with shape [nspin, nkpt, mband]."""
        return self._occfacts

    @property
    def reciprocal_lattice(self):
        """|Lattice| with the reciprocal lattice vectors in Angstrom."""
        return self.structure.reciprocal_lattice

    @property
    def shape(self):
        """Shape of the array with the eigenvalues."""
        return self.nsppol, self.nkpt, self.mband

    @property
    def has_metallic_scheme(self):
        """True if we are using a metallic scheme for occupancies."""
        return self.smearing.has_metallic_scheme

    #def recalc_fermie(self, nelect=None, method="gaussian", step=0.001, width=0.002):
    #    """
    #    Recompute the Fermi level.
    #    """
    #    if nelect is None: nelect = self.nelect
    #    edos = self.get_edos(method=method, step=step, width=width)
    #    ef = edos.find_mu(nelect)
    #    self.set_fermie(ef)
    #    return ef

    #def set_fermie(self, fermie):
    #    self.fermie = fermie
    #    # TODO: Recalculate occupations.

    def get_dict4pandas(self, with_spglib=True):
        """
        Return a :class:`OrderedDict` with the most important parameters:

            - Chemical formula and number of atoms.
            - Lattice lengths, angles and volume.
            - The spacegroup number computed by Abinit (set to None if not available).
            - The spacegroup number and symbol computed by spglib (set to None not `with_spglib`).

        Useful to construct pandas DataFrames

        Args:
            with_spglib: If True, spglib_ is invoked to get the spacegroup symbol and number.
        """
        odict = OrderedDict([
            ("nsppol", self.nsppol), ("nspinor", self.nspinor), ("nspden", self.nspden),
            ("nkpt", self.nkpt), ("nband", self.nband_sk.min()),
            ("nelect", self.nelect), ("fermie", self.fermie),

        ])
        odict.update(self.structure.get_dict4pandas(with_spglib=with_spglib))
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
        try:
            return self.kpoints.ksampling.kptopt
        except AttributeError:
            cprint("ebands.kpoints.ksampling.kptopt is not defined, assuming kptopt = 1", "red")
            return 1

    @lazy_property
    def has_timrev(self):
        """True if time-reversal symmetry is used in the BZ sampling."""
        return has_timrev_from_kptopt(self.kptopt)

    def kindex(self, kpoint):
        """
        The index of the k-point in the internal list of k-points.
        Accepts: |Kpoint| instance or integer.
        """
        if duck.is_intlike(kpoint):
            return int(kpoint)
        else:
            return self.kpoints.index(kpoint)

    #def sb_iter(self, ik):
    #    """Iterator over (spin, band) indices."""
    #    for spin in self.spins:
    #        for band in range(self.nband_sk[spin, ik]):
    #            yield spin, band

    def skb_iter(self):
        """Iterator over (spin, k, band) indices."""
        for spin in self.spins:
            for ik in self.kidxs:
                for band in range(self.nband_sk[spin, ik]):
                    yield spin, ik, band

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
            kpoint: K-point index or |Kpoint| object
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
        in the energy intervale [e0 - deltae, e0 + deltae]

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

    def get_dataframe(self, e0="fermie"):
        """
        Return a |pandas-DataFrame| with the following columns:

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

    @add_fig_kwargs
    def boxplot(self, ax=None, e0="fermie", brange=None, swarm=False, **kwargs):
        """
        Use seaborn_ to draw a box plot to show distributions of eigenvalues with respect to the band index.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV.
                -  None: Don't shift energies, equivalent to ``e0 = 0``.
            brange: Only bands such as ``brange[0] <= band_index < brange[1]`` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            kwargs: Keyword arguments passed to seaborn boxplot.

        Return: |matplotlib-Figure|
        """
        # Get the dataframe and select bands
        frame = self.get_dataframe(e0=e0)
        if brange is not None:
            frame = frame[(frame["band"] >= brange[0]) & (frame["band"] < brange[1])]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        import seaborn as sns
        hue = None if self.nsppol == 1 else "spin"
        ax = sns.boxplot(x="band", y="eig", data=frame, hue=hue, ax=ax, **kwargs)
        if swarm:
            sns.swarmplot(x="band", y="eig", data=frame, hue=hue, color=".25", ax=ax)

        return fig

    @classmethod
    def from_pymatgen(cls, pmg_bands, nelect, weights=None, has_timerev=True,
                      ksampling=None, smearing=None, nspinor=1, nspden=None):
        """
        Convert a pymatgen bandstructure object to an Abipy |ElectronBands| object.

        Args:
            pmg_bands: pymatgen bandstructure object.
            nelect: Number of electrons in unit cell.
            weights: List of K-points weights (normalized to one, same order as pmg_bands.kpoints).
                This argument is optional but recommended when ``pmg_bands`` represents an IBZ sampling.
                If weights are not provided, Abipy methods requiring integrations in the BZ won't work.
            has_timerev: True if time-reversal symmetry can be used.
            ksampling: dictionary with parameters passed to :class:`KSamplingInfo` defining the k-points sampling.
                If None, hard-coded values are used. This argument is recommended if IBZ sampling.
            smearing: dictionary with parameters passed to :class:`Smearing`
                If None, default hard-coded values are used.
            nspinor: Number of spinor components.
            nspden: Number of independent spin-density components.
                If None, nspden is automatically computed from nsppol

        .. warning::

            The Abipy bandstructure contains more information than the pymatgen object so
            the conversion is not complete, especially if you rely on the default values.
            Please read the docstring and the code carefully and use the optional arguments to pass
            additional data required by Abipy if you need a complete conversion.
        """
        from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine

        # Cast to abipy structure and call spglib to init AbinitSpaceGroup.
        abipy_structure= Structure.as_structure(pmg_bands.structure.copy())
        if not abipy_structure.has_abi_spacegroup:
            abipy_structure.spgset_abi_spacegroup(has_timerev)

        # Get dimensions.
        nsppol = 2 if pmg_bands.is_spin_polarized else 1
        if nspden is None:
            if nspinor == 1: nspden = nsppol
            if nspinor == 2: nspden = 4
        nkpt = len(pmg_bands.kpoints)

        smearing = Smearing.as_smearing(smearing)
        ksampling = KSamplingInfo.as_ksampling(ksampling)

        # Build numpy array with eigenvalues.
        abipy_eigens = np.empty((nsppol, nkpt, pmg_bands.nb_bands))
        abipy_eigens[0] = np.array(pmg_bands.bands[PmgSpin.up]).T.copy()
        if nsppol == 2:
            abipy_eigens[1] = np.array(pmg_bands.bands[PmgSpin.down]).T.copy()

        # Compute occupation factors. Note that pmg bands don't have occfact so
        # I have to compute them from the eigens assuming T=0)
        atol = 1e-4
        abipy_occfacts = np.where(abipy_eigens <= pmg_bands.efermi + atol, 1, 0)
        if nsppol == 1: abipy_occfacts *= 2

        reciprocal_lattice = pmg_bands.structure.lattice.reciprocal_lattice
        frac_coords = np.array([k.frac_coords for k in pmg_bands.kpoints])

        if isinstance(pmg_bands, BandStructureSymmLine):
            abipy_kpoints = Kpath(reciprocal_lattice, frac_coords,
                                  weights=weights, names=None, ksampling=ksampling)

        elif isinstance(pmg_bands, BandStructure):
            abipy_kpoints = IrredZone(reciprocal_lattice, frac_coords,
                                      weights=weights, names=None, ksampling=ksampling)

        else:
            raise TypeError("Don't know how to handle type: %s" % type(pmg_bands))

        # Find names of the kpoints.
        for kpoint in abipy_kpoints:
            name = abipy_structure.findname_in_hsym_stars(kpoint)

        return cls(abipy_structure, abipy_kpoints, abipy_eigens, pmg_bands.efermi, abipy_occfacts,
                   nelect, nspinor, nspden, smearing=smearing)

    def to_pymatgen(self):
        """
        Return a pymatgen bandstructure object from an Abipt |ElectronBands| object.
        """
        from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
        assert np.all(self.nband_sk == self.nband_sk[0, 0])

        # eigenvals is a dict of energies for spin up and spin down
        # {Spin.up:[][], Spin.down:[][]}, the first index of the array
        # [][] refers to the band and the second to the index of the
        # kpoint. The kpoints are ordered according to the order of the
        # kpoints array. If the band structure is not spin polarized, we
        # only store one data set under Spin.up
        eigenvals = {PmgSpin.up: self.eigens[0,:,:].T.copy().tolist()}
        if self.nsppol == 2:
            eigenvals[PmgSpin.down] = self.eigens[1,:,:].T.copy().tolist()

        if self.kpoints.is_path:
            labels_dict = {k.name: k.frac_coords for k in self.kpoints if k.name is not None}
            return BandStructureSymmLine(self.kpoints.frac_coords, eigenvals, self.reciprocal_lattice, self.fermie,
                                         labels_dict, coords_are_cartesian=False, structure=self.structure, projections=None)
        else:
            return BandStructure(self.kpoints.frac_coords, eigenvals, self.reciprocal_lattice, self.fermie,
                                 labels_dict=None, coords_are_cartesian=False, structure=self.structure, projections=None)

    def _electron_state(self, spin, kpoint, band):
        """
        Build an instance of :class:`Electron` from the spin, kpoint and band index
        """
        kidx = self.kindex(kpoint)
        #print("kidx", kidx)
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
            kpoint: Index of the kpoint or |Kpoint| object.
        """
        return self._electron_state(spin, kpoint, 0)

    def homo_sk(self, spin, kpoint):
        """
        Returns the HOMO state for the given spin, kpoint.

        Args:
            spin: Spin index
            kpoint: Index of the kpoint or |Kpoint| object.
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
            kpoint: Index of the kpoint or |Kpoint| object.
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
                b = find_gt(self.eigens[spin, k, :], self.fermie + self.pad_fermie)
                blist.append(b)
                enes.append(self.eigens[spin, k, b])

            #print("enes", enes)
            lumo_kidx = np.array(enes).argmin()
            #print("blist:", blist)
            lumo_band = blist[lumo_kidx]
            #print("lumo_band:", lumo_band)
            #print("type: lumo_kidx", type(lumo_kidx))

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

    def to_string(self, title=None, with_structure=True, with_kpoints=False, verbose=0):
        """
        Human-readable string with useful info such as band gaps, position of HOMO, LOMO...

        Args:
            with_structure: False if structural info shoud not be displayed.
            with_kpoints: False if k-point info shoud not be displayed.
            verbose: Verbosity level.
        """
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))

        if with_structure:
            app(self.structure.to_string(verbose=verbose, title="Structure"))
            app("")

        app("Number of electrons: %s, Fermi level: %.3f (eV)" % (self.nelect, self.fermie))
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
                    # This can fail so we have to catch the exception.
                    try:
                        app("Direct gap:\n%s" % indent(str(self.direct_gaps[spin])))
                        app("Fundamental gap:\n%s" % indent(str(self.fundamental_gaps[spin])))
                    except Exception as exc:
                        app("WARNING: Cannot compute direct and fundamental gap.")
                        if verbose: app("Exception:\n%s" % str(exc))

                app("Bandwidth: %.3f (eV)" % self.bandwidths[spin])
                if verbose:
                    app("Valence minimum located at:\n%s" % indent(str(self.lomos[spin])))
                app("Valence maximum located at:\n%s" % indent(str(self.homos[spin])))
                try:
                    # Cannot assume enough states for this!
                    app("Conduction minimum located at:\n%s" % indent(str(self.lumos[spin])))
                    app("")
                except Exception:
                    pass

        if with_kpoints:
            app(self.kpoints.to_string(verbose=verbose, title="K-points"))
            app("")

        return "\n".join(lines)

    def new_with_irred_kpoints(self, prune_step=None):
        """
        Return a new |ElectronBands| object in which only the irreducible k-points are kept.
        This method is mainly used to prepare the band structure interpolation as the interpolator
        will likely fail if the input k-path contains symmetrical k-points.

        Args:
            prune_step: Optional argument used to select a subset of the irreducible points found.
            If ``prune_step`` is None, all irreducible k-points are used.
        """
        # Get the index of the irreducible kpoints.
        from abipy.core.kpoints import find_irred_kpoints_generic
        nmt = find_irred_kpoints_generic(self.structure, self.kpoints.frac_coords)
        irred_map = nmt.irred_map

        if prune_step is not None:
            irred_map = irred_map[::prune_step].copy()

        # Build new set of k-points
        new_kcoords = self.kpoints.frac_coords[irred_map].copy()
        new_kpoints = KpointList(self.structure.reciprocal_lattice, new_kcoords,
                                 weights=None, names=None, ksampling=self.kpoints.ksampling)

        # Extract eigevanlues and occupation factors associated to irred k-points.
        new_eigens = self.eigens[:, irred_map, :].copy()
        new_occfacts = self.occfacts[:, irred_map, :].copy()

        return self.__class__(self.structure, new_kpoints, new_eigens, self.fermie, new_occfacts,
                              self.nelect, self.nspinor, self.nspden)

    def spacing(self, axis=None):
        """
        Compute the statistical parameters of the energy spacing, i.e. e[b+1] - e[b]

        Returns:
            ``namedtuple`` with the statistical parameters in eV
        """
        ediff = self.eigens[:, :, 1:] - self.eigens[:, :, :self.mband-1]

        return StatParams(mean=ediff.mean(axis=axis), stdev=ediff.std(axis=axis),
                          min=ediff.min(axis=axis), max=ediff.max(axis=axis))

    def statdiff(self, other, axis=None, numpy_op=np.abs):
        """
        Compare the eigenenergies of two bands and compute the
        statistical parameters: mean, standard deviation, min and max
        The bands are aligned wrt to their fermi level.

        Args:
            other: |ElectronBands| object.
            axis:  Axis along which the statistical parameters are computed.
                The default is to compute the parameters of the flattened array.
            numpy_op: Numpy function to apply to the difference of the eigenvalues. The
                      default computes ``|self.eigens - other.eigens|``.

        Returns:
            ``namedtuple`` with the statistical parameters in eV
        """
        ediff = numpy_op(self.eigens - self.fermie - other.eigens + other.fermie)
        return StatParams(mean=ediff.mean(axis=axis), stdev=ediff.std(axis=axis),
                          min=ediff.min(axis=axis), max=ediff.max(axis=axis))

    def ipw_edos_widget(self): # pragma: no cover
        """
        Return an ipython widget with controllers to compute the electron DOS.
        """
        def plot_dos(method, step, width):
            edos = self.get_edos(method=method, step=step, width=width)
            edos.plot()

        import ipywidgets as ipw
        return ipw.interact_manual(
                plot_dos,
                method=["gaussian", "tetra"],
                step=ipw.FloatSlider(value=0.1, min=1e-6, max=1, step=0.05, description="Step of linear mesh (eV)"),
                width=ipw.FloatSlider(value=0.2, min=1e-6, max=1, step=0.05, description="Gaussian broadening (eV)"),
            )

    def get_edos(self, method="gaussian", step=0.1, width=0.2):
        """
        Compute the electronic DOS on a linear mesh.

        Args:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Returns: |ElectronDos| object.
        """
        self.kpoints.check_weights()

        # Compute linear mesh.
        epad = 3.0 * width
        e_min = self.enemin() - epad
        e_max = self.enemax() + epad
        nw = int(1 + (e_max - e_min) / step)
        mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)

        # TODO: Write cython version.
        dos = np.zeros((self.nsppol, nw))
        if method == "gaussian":
            for spin in self.spins:
                for k, kpoint in enumerate(self.kpoints):
                    weight = kpoint.weight
                    for band in range(self.nband_sk[spin,k]):
                        e = self.eigens[spin,k,band]
                        dos[spin] += weight * gaussian(mesh, width, center=e)

        else:
            raise NotImplementedError("Method %s is not supported" % method)

        # Use fermie from Abinit if we are not using metallic scheme for occopt.
        fermie = None
        #if self.smearing["occopt"] == 1:
        #    print("using fermie from GSR")
        #    fermie = self.fermie
        edos = ElectronDos(mesh, dos, self.nelect, fermie=fermie)
        #print("ebands.fermie", self.fermie, "edos.fermie", edos.fermie)
        return edos

    def compare_gauss_edos(self, widths, step=0.1):
        """
        Compute the electronic DOS with the Gaussian method for different values
        of the broadening. Return plotter object.

        Args:
            widths: List with the tandard deviation (eV) of the gaussian.
            step: Energy step (eV) of the linear mesh.

        Return: |ElectronDosPlotter|
        """
        edos_plotter = ElectronDosPlotter()
        for width in widths:
           edos = self.get_edos(method="gaussian", step=0.1, width=width)
           label=r"$\sigma = %s$ (eV)" % width
           edos_plotter.add_edos(label, edos)

        return edos_plotter

    @add_fig_kwargs
    def plot_transitions(self, omega_ev, qpt=(0, 0, 0), atol_ev=0.1, atol_kdiff=1e-4,
                         ylims=None, ax=None, alpha=0.4, **kwargs):
        """
        Plot energy bands with arrows signaling possible k --> k + q indipendent-particle transitions
        of energy ``omega_ev`` connecting occupied to empty states.

        Args:
            omega_ev: Transition energy in eV.
            qpt: Q-point in reduced coordinates.
            atol_ev: Absolute tolerance for energy difference in eV
            atol_kdiff: Tolerance used to compare k-points.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                or scalar e.g. `left`. If left (right) is None, default values are used
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        e0 = self.get_e0("fermie")
        self.plot(ax=ax, e0=e0, ylims=ylims, show=False)

        # Pre-compute mapping k_index --> (k + q)_index, g0
        k2kqg = self.kpoints.get_k2kqg_map(qpt, atol_kdiff=atol_kdiff)

        # Add arrows to the plot (different colors for spin up/down)
        from matplotlib.patches import FancyArrowPatch
        for spin in self.spins:
            cachek = {}
            arrow_opts = {"color": "k"} if spin == 0 else {"color": "red"}
            for ik, (ikq, g0) in k2kqg.items():
                dx = ikq - ik
                ek = self.eigens[spin, ik]
                ekq = self.eigens[spin, ikq]
                # Find rightmost value less than or equal to fermie.
                nv_k = cachek.get(ik)
                if nv_k is None:
                    nv_k = find_le(ek, self.fermie + self.pad_fermie)
                    cachek[ik] = nv_k

                if ik == ikq:
                    nv_kq = nv_k
                else:
                    nv_kq = cachek.get(ikq)
                    if nv_kq is None:
                        nv_kq = find_le(ekq, self.fermie + self.pad_fermie)
                        cachek[ikq] = nv_kq

                #print("nv_k:", nv_k, "nc_kq", nv_kq)
                for v_k in range(nv_k):
                    for c_kq in range(nv_kq + 1, self.nband):
                        dy = self.eigens[spin, ikq, c_kq] - self.eigens[spin, ik, v_k]
                        if abs(dy - omega_ev) > atol_ev: continue
                        y = self.eigens[spin, ik, v_k] - e0
                        # http://matthiaseisen.com/matplotlib/shapes/arrow/
                        p = FancyArrowPatch((ik, y), (ik + dx, y + dy),
                                connectionstyle='arc3', mutation_scale=20,
                                alpha=alpha, **arrow_opts)
                        ax.add_patch(p)
        return fig

    def get_ejdos(self, spin, valence, conduction, method="gaussian", step=0.1, width=0.2, mesh=None):
        r"""
        Compute the join density of states at q == 0.

            :math:`\sum_{kbv} f_{vk} (1 - f_{ck}) \delta(\omega - E_{ck} + E_{vk})`

        .. warning::

            The present implementation assumes an energy gap

        Args:
            spin: Spin index.
            valence: Int or iterable with the valence indices.
            conduction: Int or iterable with the conduction indices.
            method (str): String defining the integraion method.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            mesh: Frequency mesh to use. If None, the mesh is computed automatically from the eigenvalues.

        Returns: |Function1D| object.
        """
        # TODO: Generalize to k+q with
        # k2kqg = self.kpoints.get_k2kqg_map(qpt, atol_kdiff=atol_kdiff)
        self.kpoints.check_weights()
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
                    ec = self.eigens[spin, k, c]
                    fc = 1.0 - self.occfacts[spin,k,c] / full
                    for v in valence:
                        ev = self.eigens[spin, k, v]
                        fv = self.occfacts[spin, k, v] / full
                        fact = weight * fv * fc
                        jdos += fact * gaussian(mesh, width, center=ec-ev)

        else:
            raise NotImplementedError("Method %s is not supported" % str(method))

        return Function1D(mesh, jdos)

    @add_fig_kwargs
    def plot_ejdosvc(self, vrange, crange, method="gaussian", step=0.1, width=0.2, colormap="jet",
                     cumulative=True, ax=None, alpha=0.7, fontsize=12, **kwargs):
        """
        Plot the decomposition of the joint-density of States (JDOS).

        .. warning::

            The present implementation assumes an energy gap

        Args:
            vrange: Int or `Iterable` with the indices of the valence bands to consider.
            crange: Int or `Iterable` with the indices of the conduction bands to consider.
            method: String defining the method.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            colormap: Have a look at the colormaps here and decide which one you like:
                http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html
            cumulative: True for cumulative plots (default).
            ax: |matplotlib-Axes| or None if a new figure should be created.
            alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        if not isinstance(crange, Iterable): crange = [crange]
        if not isinstance(vrange, Iterable): vrange = [vrange]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)
        ax.set_xlabel('Energy (eV)')
        cmap = plt.get_cmap(colormap)
        lw = kwargs.pop("lw", 1.0)

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
                    color = cmap(float(i) / num_plots)
                    x, y = jdos.mesh, jdos.values
                    ax.plot(x, cumulative + y, lw=lw, label=label, color=color)
                    ax.fill_between(x, cumulative, cumulative + y, facecolor=color, alpha=alpha)
                    cumulative += jdos.values
                    i += 1
            else:
                num_plots, i = len(jdos_vc), 0
                for (v, c), jdos in jdos_vc.items():
                    color = cmap(float(i) / num_plots)
                    jdos.plot_ax(ax, color=color, lw=lw,
                        label=r"$v=%s \rightarrow c=%s, \sigma=%s$" % (v, c, s))
                    i += 1

            tot_jdos.plot_ax(ax, color="k", lw=lw, label=r"Total JDOS, $\sigma=%s$" % s)

        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def apply_scissors(self, scissors):
        """
        Modify the band structure with the scissors operator.

        Args:
            scissors: An instance of :class:`Scissors`.

        Returns:
            New instance of |ElectronBands| with modified energies.
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
            nband_sk=self.nband_sk, smearing=self.smearing)

    @add_fig_kwargs
    def plot(self, spin=None, band_range=None, klabels=None, e0="fermie", ax=None, ylims=None,
	     points=None, **kwargs):
        r"""
        Plot the electronic band structure.

        Args:
            spin: Spin index. None to plot both spins.
            band_range: Tuple specifying the minimum and maximum band to plot (default: all bands are plotted)
            klabels: dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels. e.g.
                ``klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0):"L"}``.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            points:

        Returns: |matplotlib-Figure|
        """
        # Select spins
        spin_list = self.spins if spin is None else [spin]

        # Select the band range.
        if band_range is None:
            band_list = list(range(self.mband))
        else:
            # This does not work in py2.7 because range is not a class
            #if not isinstance(band_range, range):
            #    band_list = list(band_range)
            band_list = list(range(band_range[0], band_range[1], 1))

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, klabels=klabels)
        set_axlims(ax, ylims, "y")

        # Plot the band energies.
        for spin in spin_list:
            if spin == 0:
                opts = {"color": "black", "linewidth": 2.0}
            else:
                opts = {"color": "red", "linewidth": 2.0}

            for band in band_list:
                self.plot_ax(ax, e0, spin=spin, band=band, **opts)

        if points is not None:
            e0 = self.get_e0(e0)
            ax.scatter(points.x, np.array(points.y) - e0, s=np.abs(points.s),
                       marker="o", c="b")

        return fig

    @add_fig_kwargs
    def plot_scatter3d(self, band, spin=0, e0="fermie", colormap="jet", ax=None, **kwargs):
        r"""
        Use matplotlib ``scatter3D`` to produce a scatter plot of the eigenvalues in 3D.
        The color of the points gives the energy of the state wrt to the Fermi level.

        Args:
            band: Band index
            spin: Spin index.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            colormap: Have a look at the colormaps here and decide which one you like:
                <http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html>
            ax: matplotlib :class:`Axes3D` or None if a new figure should be created.
        """
        kcart_coords = self.kpoints.get_cart_coords()
        c = self.eigens[spin, :, band] - self.get_e0(e0)

        ax, fig, plt = get_ax3d_fig_plt(ax)
        cmap = plt.get_cmap(colormap)
        #ax.scatter3D(xs, ys, zs, s=6, alpha=0.8, marker=',', facecolors=cmap(N), lw=0)
        p = ax.scatter3D(kcart_coords[:, 0], kcart_coords[:, 1], zs=kcart_coords[:, 2], zdir='z',
                         s=20, c=c, depthshade=True, cmap=cmap)

        #self.structure.plot_bz(ax=ax, pmg_path=False, with_labels=False, show=False, linewidth=0)
        from pymatgen.electronic_structure.plotter import plot_wigner_seitz
        plot_wigner_seitz(self.structure.reciprocal_lattice, ax=ax, linewidth=1)
        ax.set_xlabel("$K_x$")
        ax.set_ylabel("$K_y$")
        ax.set_zlabel("$K_z$")
        fig.colorbar(p)

        #ax.set_title(structure.composition.formula)
        ax.set_axis_off()

        return fig

    def decorate_ax(self, ax, **kwargs):
        """
        Add k-labels, title and unit name to axis ax.

        Args:
            title:
            fontsize
            klabels:
            klabel_size:
        """
        title = kwargs.pop("title", None)
        fontsize = kwargs.pop("fontsize", 12)
        if title is not None: ax.set_title(title, fontsize=fontsize)

        ax.grid(True)
        ax.set_ylabel("Energy (eV)")
        ax.set_xlabel("Wave vector")

        # Set ticks and labels.
        klabels = kwargs.pop("klabels", None)
        ticks, labels = self._make_ticks_and_labels(klabels)
        if ticks:
            # Don't show label if previous k-point is the same.
            for il in range(1, len(labels)):
                if labels[il] == labels[il-1]: labels[il] = ""
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))

            #print("ticks", len(ticks), ticks)
            ax.set_xlim(ticks[0], ticks[-1])

    def get_e0(self, e0):
        """
        e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``.
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
            ax: |matplotlib-Axes|.
            e0: Option used to define the zero of energy in the band structure plot.
            spin: Spin index. If None, all spins are plotted.
            band: Band index, If None, all bands are plotted.
            kwargs: Passed to ax.plot

        Return: matplotlib lines
        """
        spin_range = range(self.nsppol) if spin is None else [spin]
        band_range = range(self.mband) if band is None else [band]

        label = kwargs.pop("label", None)
        # Handle linewidths
        with_linewidths = kwargs.pop("with_linewidths", False) and self.has_linewidths
        if with_linewidths:
            lw_opts = kwargs.pop("lw_opts", dict(alpha=0.6))
            lw_fact = lw_opts.pop("fact", 2.0)

        xx, lines = np.arange(self.nkpt), []
        e0 = self.get_e0(e0)
        for spin in spin_range:
            for band in band_range:
                yy = self.eigens[spin, :, band] - e0

                # Set label only at the first iteration
                lines.extend(ax.plot(xx, yy, label=label, **kwargs))
                label = None

                if with_linewidths:
                    w = self.linewidths[spin, :, band] * lw_fact / 2
                    lw_color = lines[-1].get_color()
                    ax.fill_between(xx, yy - w, yy + w, facecolor=lw_color, **lw_opts)
                    #, alpha=self.alpha, facecolor=self.l2color[l])

        return lines

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
    def plot_with_edos(self, edos, klabels=None, ax_list=None, e0="fermie", ylims=None, width_ratios=(2, 1), **kwargs):
        r"""
        Plot the band structure and the DOS.

        Args:
            edos: An instance of |ElectronDos|.
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. ``klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            ax_list: The axes for the bandstructure plot and the DOS plot. If ax_list is None, a new figure
                is created and the two axes are automatically generated.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            e0: Option used to define the zero of energy in the band structure plot. Possible values::

                * ``fermie``: shift all eigenvalues and the DOS to have zero energy at the Fermi energy.
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   This should be ok in semiconductors. In metals, however, a better value of the Fermi energy
                   can be obtained from the DOS provided that the k-sampling for the DOS is much denser than
                   the one used to compute the density. See ``edos_fermie``.
                * ``edos_fermie``: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                *  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                *  None: Don't shift energies, equivalent to ``e0 = 0``

            width_ratios: Defines the ratio between the band structure plot and the dos plot.

        Return: |matplotlib-Figure|
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        if ax_list is None:
            # Build axes and align bands and DOS.
            fig = plt.figure()
            gspec = GridSpec(nrows=1, ncols=2, width_ratios=width_ratios, wspace=0.05)
            ax0 = plt.subplot(gspec[0])
            ax1 = plt.subplot(gspec[1], sharey=ax0)
        else:
            # Take them from ax_list.
            ax0, ax1 = ax_list
            fig = plt.gcf()

        # Define the zero of energy.
        e0 = self.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
        #if not kwargs: kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the band structure
        for spin in self.spins:
            if spin == 0:
                opts = {"color": "black", "linewidth": 2.0}
            else:
                opts = {"color": "red", "linewidth": 2.0}

            for band in range(self.mband):
                self.plot_ax(ax0, e0, spin=spin, band=band, **opts)

        self.decorate_ax(ax0, klabels=klabels)
        set_axlims(ax0, ylims, "y")

        # Plot the DOS
        if self.nsppol == 1:
            opts = {"color": "black", "linewidth": 2.0}
            edos.plot_ax(ax1, e0, exchange_xy=True, **opts)
        else:
            for spin in self.spins:
                if spin == 0:
                    opts = {"color": "black", "linewidth": 2.0}
                else:
                    opts = {"color": "red", "linewidth": 2.0}
                edos.plot_ax(ax1, e0, spin=spin, exchange_xy=True, **opts)

        ax1.grid(True)
        ax1.yaxis.set_ticks_position("right")
        ax1.yaxis.set_label_position("right")
        set_axlims(ax1, ylims, "y")

        return fig

    @add_fig_kwargs
    def plot_lws_vs_e0(self, ax=None, e0="fermie", exchange_xy=False,
                       xlims=None, ylims=None, fontsize=12, **kwargs):
        r"""
        Plot the electronic linewidths vs KS energy.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            exchange_xy: True to exchange x-y axis.
            xlims, ylims: Set the data limits for the x-axis or the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: fontsize for titles and legend.

        Returns: |matplotlib-Figure|
        """
        if not self.has_linewidths: return None
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        xlabel = r"$\epsilon_{KS}\;(eV)$"
        if e0 is not None:
            xlabel = r"$\epsilon_{KS}-\epsilon_F\;(eV)$"

        # DSU sort to get lw(e) with sorted energies.
        e0mesh, lws = zip(*sorted(zip(self.eigens.flat, self.linewidths.flat), key=lambda t: t[0]))
        e0 = self.get_e0(e0)
        e0mesh = np.array(e0mesh) - e0

        kw_linestyle = kwargs.pop("linestyle", "o")
        #kw_lw = kwargs.pop("lw", 1)
        #kw_lw = kwargs.pop("markersize", 5)
        kw_color = kwargs.pop("color", "red")
        kw_label = kwargs.pop("label", None)

        xx, yy = e0mesh, lws
        if exchange_xy: xx, yy = yy, xx
        ax.plot(xx, yy, kw_linestyle, color=kw_color, label=kw_label, **kwargs)
        #ax.scatter(xx, yy)

        ax.grid(True)
        ax.set_ylabel("Linewidth")
        ax.set_xlabel(xlabel)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        if kw_linestyle:
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def to_xmgrace(self, filepath):
        """
        Write xmgrace_ file with band structure energies and labels for high-symmetry k-points.

        Args:
            filepath: String with filename or stream.
        """
        is_stream = hasattr(filepath, "write")
        if is_stream:
            f = filepath
        else:
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
        w("# Energies are in eV. Zero set to efermi, previously it was at: %s (eV)" % self.fermie)
        w("# List of k-points and their index (C notation i.e. count from 0)")
        for ik, kpt in enumerate(self.kpoints):
            w("# %d %s" % (ik, str(kpt.frac_coords)))
        w("@page size 792, 612")
        w("@page scroll 5%")
        w("@page inout 5%")
        w("@link page off")
        w("@with g0")
        w("@world xmin 0.00")
        w('@world xmax %d' % (self.nkpt - 1))
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
        w('@yaxis  label "Band Energy (eV)"')
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

        if not is_stream:
            f.close()

    def to_bxsf(self, filepath):
        """
        Export the full band structure to ``filepath`` in BXSF format
        suitable for the visualization of isosurfaces with xcrysden_ (xcrysden --bxsf FILE).
        Require k-points in IBZ and gamma-centered k-mesh.
        """
        self.get_ebands3d().to_bxsf(filepath)

    def get_ebands3d(self):
        return ElectronBands3D(self.structure, self.kpoints, self.has_timrev, self.eigens, self.fermie)

    def derivatives(self, spin, band, order=1, acc=4):
        """
        Compute the derivative of the eigenvalues wrt to k.

        Args:
            spin: Spin index
            band: Band index
            order:
            acc:

        Returns:
        """
        if self.kpoints.is_path:
            # Extract the energy branch.
            ebranch = self.eigens[spin, :, band]
            # Simulate free-electron bands. This will produce all(effective masses == 1)
            #ebranch = 0.5 * units.Ha_to_eV * np.array([(k.norm * units.bohr_to_ang)**2 for k in self.kpoints])

            # Compute derivatives by finite differences.
            ders_onlines = self.kpoints.finite_diff(ebranch, order=order, acc=acc)
            return ders_onlines

        else:
            raise NotImplementedError("Derivatives on homogeneous k-meshes are not supported yet")

    def effective_masses(self, spin, band, acc=4):
        """
        Compute the effective masses for the given ``spin`` and ``band`` index.
        Use finite difference with accuracy ``acc``.

        Returns:
            |numpy-array| of size self.nkpt with effective masses.
        """
        ders2 = self.derivatives(spin, band, order=2, acc=acc) * (units.eV_to_Ha / units.bohr_to_ang**2)
        return 1. / ders2

    def effmass_line(self, spin, kpoint, band, acc=4):
        """
        Compute the effective masses along a line. Requires band energies on a k-path.

        Args:
            spin: Spin index.
            kpoint: integer or |Kpoint| object. Note that if kpoint is not an integer,
                and the path contains duplicated k-points, the first k-point is selected.
            band: Band index.
            acc: accuracy
        """
        if not self.kpoints.is_path:
            raise ValueError("effmass_line requires points along a path.")

        warnings.warn("This code is still under development. API may change!")

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
        evals_on_line, h_left, vers_left = self._eigens_hvers_iline(spin, band, iline)
        d2line = finite_diff(evals_on_line, h_left, order=2, acc=acc) * (units.eV_to_Ha / units.bohr_to_ang**2)
        em_left = 1. / d2line[kpos]
        em_right = em_left
        h_right, vers_right = h_left, vers_left

        if do_right:
            kpos_right = self.kpoints.lines[iline+1].index(ik)
            assert kpos_right == 0
            evals_on_line, h_right, vers_right = self._eigens_hvers_iline(spin, band, iline+1)
            d2line = finite_diff(evals_on_line, h_right, order=2, acc=acc) * (units.eV_to_Ha / units.bohr_to_ang**2)
            em_right = 1. / d2line[kpos_right]

        return EffectiveMassAlongLine(spin, self.kpoints[ik], band, self.eigens[spin, ik, band],
                                      acc, self.structure.reciprocal_lattice,
                                      is_inside, h_left, vers_left, em_left, h_right, vers_right, em_right)

    def _eigens_hvers_iline(self, spin, band, iline):
        line = self.kpoints.lines[iline]
        evals_on_line = self.eigens[spin, line, band]
        h = self.kpoints.ds[line[0]]

        if not np.allclose(h, self.kpoints.ds[line[:-1]]):
            raise ValueError("For finite difference derivatives, the path must be homogeneous!\n" +
                             str(self.kpoints.ds[line[:-1]]))

        return evals_on_line, h, self.kpoints.versors[line[0]]

    def interpolate(self, lpratio=5, vertices_names=None, line_density=20,
                    kmesh=None, is_shift=None, filter_params=None, verbose=0):
        """
        Interpolate energies in k-space along a k-path and, optionally, in the IBZ for DOS calculations.
        Note that the interpolation will likely fail if there are symmetrical k-points in the input set of k-points
        so it's recommended to call this method with band structure obtained in the IBZ.

        Args:
            lpratio: Ratio between the number of star functions and the number of ab-initio k-points.
                The default should be OK in many systems, larger values may be required for accurate derivatives.
            vertices_names: Used to specify the k-path for the interpolated band structure
                It's a list of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path. ``line_density`` defines
                the density of the sampling. If None, the k-path is automatically generated according
                to the point group of the system.
            line_density: Number of points in the smallest segment of the k-path. Used with ``vertices_names``.
            kmesh: Used to activate the interpolation on the homogeneous mesh for DOS (uses spglib_ API).
                kmesh is given by three integers and specifies mesh numbers along reciprocal primitive axis.
            is_shift: three integers (spglib_ API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshited mesh.
            filter_params: TO BE described.
            verbose: Verbosity level

        Returns:
                namedtuple with the following attributes::

                    ebands_kpath: |ElectronBands| with the interpolated band structure on the k-path.
                    ebands_kmesh: |ElectronBands| with the interpolated band structure on the k-mesh.
                        None if ``kmesh`` is not given.
                    interpolator: |SkwInterpolator| object.
        """
        # Get symmetries from abinit spacegroup (read from file).
        abispg = self.structure.abi_spacegroup
        if abispg is None:
            abispg = self.structure.spgset_abi_spacegroup(has_timerev=self.has_timrev)

        fm_symrel = [s for (s, afm) in zip(abispg.symrel, abispg.symafm) if afm == 1]

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
            ksampling = KSamplingInfo.from_mpdivs(mpdivs=kmesh, shifts=[0, 0, 0], kptopt=1)
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
        app("K-point: %s, eigenvalue: %s (eV)" % (self.kpoint, self.eig))
        app("h_left: %s, h_right %s" % (self.h_left, self.h_right))
        app("is_inside: %s, vers_left: %s, vers_right: %s" % (self.is_inside, self.vers_left, self.vers_right))
        app("em_left: %s, em_right: %s" % (self.em_left, self.em_right))
        return "\n".join(lines)


def dataframe_from_ebands(ebands_objects, index=None, with_spglib=True):
    """
    Build a pandas dataframe with the most important results available in a list of band structures.

    Args:
        ebands_objects: List of objects that can be converted to structure.
            Support netcdf filenames or |ElectronBands| objects
            See ``ElectronBands.as_ebands`` for the complete list.
        index: Index of the dataframe.
        with_spglib: If True, spglib is invoked to get the spacegroup symbol and number.

    Return: |pandas-DataFrame|
    """
    ebands_list = [ElectronBands.as_ebands(obj) for obj in ebands_objects]
    # Use OrderedDict to have columns ordered nicely.
    odict_list = [(ebands.get_dict4pandas(with_spglib=with_spglib)) for ebands in ebands_list]

    import pandas as pd
    return pd.DataFrame(odict_list, index=index)
                        #columns=list(odict_list[0].keys()) if odict_list else None)


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

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ElectronBandsPlotter
    """
    # Used in iter_lineopt to generate matplotlib linestyles.
    _LINE_COLORS = ["b", "r", "g", "m", "y", "k"]
    _LINE_STYLES = ["-",":","--","-.",]
    _LINE_WIDTHS = [2,]

    def __init__(self, key_ebands=None, key_edos=None, edos_kwargs=None):
        """
        Args:
            key_ebands: List of (label, ebands) tuples.
                ebands is any object that can be converted into |ElectronBands| e.g. ncfile, path.
            key_edos: List of (label, edos) tuples.
                edos is any object that can be converted into |ElectronDos|.
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
        return self.to_string(func=repr)

    def __str__(self):
        """Invoked by str"""
        return self.to_string(func=str)

    def add_plotter(self, other):
        """Merge two plotters, return new plotter."""
        if not isinstance(other, self.__class__):
            raise TypeError("Don't know to to add %s to %s" % (other.__class__, self.__class__))

        key_ebands = list(self.ebands_dict.items()) + list(other.ebands_dict.items())
        key_edos = list(self.edoses_dict.items()) + list(other.edoses_dict.items())

        return self.__class__(key_ebands=key_ebands, key_edos=key_edos)

    def to_string(self, func=str, verbose=0):
        """String representation."""
        lines = []
        app = lines.append
        for i, (label, ebands) in enumerate(self.ebands_dict.items()):
            app("[%d] %s --> %s" % (i, label, func(ebands)))

        if self.edoses_dict:
            for i, (label, edos) in enumerate(self.edoses_dict.items()):
                app("[%d] %s --> %s" % (i, label, func(edos)))

        return "\n".join(lines)

    def get_ebands_frame(self, with_spglib=True):
        """
        Build a |pandas-DataFrame| with the most important results available in the band structures.
        Useful to analyze band-gaps.
        """
        return dataframe_from_ebands(list(self.ebands_dict.values()),
                                     index=list(self.ebands_dict.keys()), with_spglib=with_spglib)

    @property
    def ebands_list(self):
        """"List of |ElectronBands| objects."""
        return list(self.ebands_dict.values())

    @property
    def edoses_list(self):
        """"List of |ElectronDos| objects."""
        return list(self.edoses_dict.values())

    def iter_lineopt(self):
        """Generates matplotlib linestyles."""
        for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def add_ebands(self, label, bands, edos=None, dos=None, edos_kwargs=None):
        """
        Adds a band structure and optionally a edos to the plotter.

        Args:
            label: label for the bands. Must be unique.
            bands: |ElectronBands| object.
            edos: |ElectronDos| object.
            edos_kwargs: optional dictionary with the options passed to ``get_edos`` to compute the electron DOS.
                Used only if ``edos`` is not None and it's not an |ElectronDos| instance.
        """
        if dos is not None:
            warnings.warn("dos has been replaced by edos! This argument will be removed in v0.4")
            assert edos is None
            edos = dos

        if label in self.ebands_dict:
            raise ValueError("label %s is already in %s" % (label, list(self.ebands_dict.keys())))

        self.ebands_dict[label] = ElectronBands.as_ebands(bands)
        if edos is not None:
            self.edoses_dict[label] = ElectronDos.as_edos(edos, edos_kwargs)

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

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        for mname in ("gridplot", "boxplot"):
            yield getattr(self, mname)(show=False)

    @add_fig_kwargs
    def combiplot(self, e0="fermie", ylims=None, width_ratios=(2, 1), fontsize=8, **kwargs):
        """
        Plot the band structure and the DOS on the same figure.
        Use ``gridplot`` to plot band structures on different figures.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values::

                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (ebands.fermie)
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   This should be ok in semiconductors. In metals, however, a better value of the Fermi energy
                   can be obtained from the DOS provided that the k-sampling for the DOS is much denser than
                   the one used to compute the density. See `edos_fermie`.
                - ``edos_fermie``: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                   Available only if plotter contains dos objects.
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0

            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            width_ratios: Defines the ratio between the band structure plot and the dos plot.
                Used when there are DOS stored in the plotter.
            fontsize: fontsize for titles and legend.

        Returns: |matplotlib-Figure|.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        fig = plt.figure()

        if self.edoses_dict:
            # Build grid with two axes.
            gspec = GridSpec(nrows=1, ncols=2, width_ratios=width_ratios, wspace=0.05)
            # bands and DOS will share the y-axis
            ax0 = plt.subplot(gspec[0])
            ax1 = plt.subplot(gspec[1], sharey=ax0)
            ax_list = [ax0, ax1]
        else:
            # One axis for bands only
            ax0 = fig.add_subplot(111)
            ax_list = [ax0]

        for ax in ax_list:
            ax.grid(True)
            set_axlims(ax, ylims, "y")

        # Plot ebands.
        lines, legends = [], []
        my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        nkpt_list = [ebands.nkpt for ebands in self.ebands_dict.values()]
        if any(nk != nkpt_list[0] for nk in nkpt_list):
            cprint("WARNING: Bands have different number of k-points:\n%s" % str(nkpt_list), "yellow")

        for (label, ebands), lineopt in zip(self.ebands_dict.items(), self.iter_lineopt()):
            i += 1
            my_kwargs.update(lineopt)
            opts_label[label] = my_kwargs.copy()

            # Get energy zero.
            if e0 == "edos_fermie":
                mye0 = self.edoses_dict[label].fermie
            else:
                mye0 = ebands.get_e0(e0)

            l = ebands.plot_ax(ax0, mye0, spin=None, band=None, **my_kwargs)
            lines.append(l[0])

            # Use relative paths if label is a file.
            if os.path.isfile(label):
                legends.append("%s" % os.path.relpath(label))
            else:
                legends.append("%s" % label)

            # Set ticks and labels, legends.
            if i == 0:
                ebands.decorate_ax(ax0)

        ax0.legend(lines, legends, loc='upper right', fontsize=fontsize, shadow=True)

        # Add DOSes
        if self.edoses_dict:
            ax = ax_list[1]
            for label, edos in self.edoses_dict.items():
                ebands = self.edoses_dict[label]
                mye0 = ebands.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
                edos.plot_ax(ax, mye0, exchange_xy=True, **opts_label[label])

        return fig

    def plot(self, *args, **kwargs):
        """An alias for combiplot."""
        if "align" in kwargs or "xlim" in kwargs or "ylim" in kwargs:
            raise ValueError("align|xlim|ylim options are not supported anymore.")
        return self.combiplot(*args, **kwargs)

    @add_fig_kwargs
    def gridplot(self, e0="fermie", with_dos=True, ylims=None, fontsize=8, **kwargs):
        """
        Plot multiple electron bandstructures and optionally DOSes on a grid.

        Args:
            eb_objects: List of objects from which the band structures are extracted.
                Each item in eb_objects is either a string with the path of the netcdf file,
                or one of the abipy object with an ``ebands`` attribute or a |ElectronBands| object.
            edos_objects: List of objects from which the electron DOSes are extracted.
                Accept filepaths or |ElectronDos| objects. If edos_objects is not None,
                each subplot in the grid contains a band structure with DOS else a simple bandstructure plot.
            e0: Option used to define the zero of energy in the band structure plot. Possible values::

                - ``fermie``: shift all eigenvalues and the DOS to have zero energy at the Fermi energy.
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   This should be ok in semiconductors. In metals, however, a better value of the Fermi energy
                   can be obtained from the DOS provided that the k-sampling for the DOS is much denser than
                   the one used to compute the density. See `edos_fermie`.
                - ``edos_fermie``: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                   Available only if edos_objects is not None
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0

            with_dos: True if DOS should be printed.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ```(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: fontsize for titles and legend.

        Returns: |matplotlib-Figure|
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
            fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharey=True, squeeze=False)
            ax_list = ax_list.ravel()
            # don't show the last ax if numeb is odd.
            if numeb % ncols != 0: ax_list[-1].axis("off")

            for i, (ebands, ax) in enumerate(zip(ebands_list, ax_list)):
                irow, icol = divmod(i, ncols)
                ebands.plot(ax=ax, e0=e0, show=False)
                set_axlims(ax, ylims, "y")
                if titles is not None: ax.set_title(titles[i], fontsize=fontsize)
                if (irow, icol) != (0, 0):
                    set_visible(ax, False, "ylabel")

        else:
            # Plot grid with bands + DOS. see http://matplotlib.org/users/gridspec.html
            from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
            fig = plt.figure()
            gspec = GridSpec(nrows, ncols)

            for i, (ebands, edos) in enumerate(zip(ebands_list, edos_list)):
                subgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=gspec[i], width_ratios=[2, 1], wspace=0.05)
                # Get axes and align bands and DOS.
                ax0 = plt.subplot(subgrid[0])
                ax1 = plt.subplot(subgrid[1], sharey=ax0)
                set_axlims(ax0, ylims, "y")
                set_axlims(ax1, ylims, "y")

                # Define the zero of energy and plot
                mye0 = ebands.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
                ebands.plot_with_edos(edos, e0=mye0, ax_list=(ax0, ax1), show=False)

                if titles is not None: ax0.set_title(titles[i], fontsize=fontsize)
                if i % ncols != 0:
                    for ax in (ax0, ax1):
                        ax.set_ylabel("")

        return fig

    @add_fig_kwargs
    def boxplot(self, e0="fermie", brange=None, swarm=False, fontsize=8, **kwargs):
        """
        Use seaborn_ to draw a box plot to show distributions of eigenvalues with respect to the band index.
        Band structures are drawn on different subplots.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            brange: Only bands such as ``brange[0] <= band_index < brange[1]`` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            fontsize: Fontsize for title.
            kwargs: Keyword arguments passed to seaborn boxplot.
        """
        # Build grid of plots.
        num_plots, ncols, nrows = len(self.ebands_dict), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=True, squeeze=False)
        ax_list = ax_list.ravel()

        # don't show the last ax if numeb is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        for (label, ebands), ax in zip(self.ebands_dict.items(), ax_list):
            ebands.boxplot(ax=ax, brange=brange, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def combiboxplot(self, e0="fermie", brange=None, swarm=False, ax=None, **kwargs):
        """
        Use seaborn_ to draw a box plot comparing the distributions of the eigenvalues
        Band structures are drawn on the same plot.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0

            brange: Only bands such as ``brange[0] <= band_index < brange[1]`` are included in the plot.
            swarm: True to show the datapoints on top of the boxes
            ax: |matplotlib-Axes| or None if a new figure should be created.
            kwargs: Keyword arguments passed to seaborn boxplot.
        """
        spin_polarized = False
        frames = []
        for label, ebands in self.ebands_dict.items():
            # Get the dataframe, select bands and add column with label
            frame = ebands.get_dataframe(e0=e0)
            if brange is not None:
                frame = frame[(frame["band"] >= brange[0]) & (frame["band"] < brange[1])]
            frame["label"] = label
            frames.append(frame)
            if ebands.nsppol == 2: spin_polarized = True

        # Merge frames ignoring index (not meaningful)
        import pandas as pd
        data = pd.concat(frames, ignore_index=True)

        import matplotlib.pyplot as plt
        import seaborn as sns
        if not spin_polarized:
            ax, fig, plt = get_ax_fig_plt(ax=ax)
            ax.grid(True)
            sns.boxplot(x="band", y="eig", data=data, hue="label", ax=ax, **kwargs)
            if swarm:
                sns.swarmplot(x="band", y="eig", data=data, hue="label", color=".25", ax=ax)
        else:
            # Generate two subplots for spin-up / spin-down channels.
            if ax is not None:
                raise NotImplementedError("ax == None not implemented when nsppol==2")
            fig, ax_list = plt.subplots(nrows=2, ncols=1, sharex=True, squeeze=False)
            for spin, ax in zip(range(2), ax_list.ravel()):
                ax.grid(True)
                data_spin = data[data["spin"] == spin]
                sns.boxplot(x="band", y="eig", data=data_spin, hue="label", ax=ax, **kwargs)
                if swarm:
                    sns.swarmplot(x="band", y="eig", data=data_spin, hue="label", color=".25", ax=ax)

        return fig

    def animate(self, e0="fermie", interval=500, savefile=None, width_ratios=(2, 1), show=True):
        """
        Use matplotlib_ to animate a list of band structure plots (with or without DOS).

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values::

                * ``fermie``: shift all eigenvalues and the DOS to have zero energy at the Fermi energy.
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   See `edos_fermie`.
                * ``edos_fermie``: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                *  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                *  None: Don't shift energies, equivalent to e0=0

            interval: draws a new frame every interval milliseconds.
            savefile: Use e.g. 'myanimation.mp4' to save the animation in mp4 format.
            width_ratios: Defines the ratio between the band structure plot and the dos plot.
                Used when there are DOS stored in the plotter.
            show: True if the animation should be shown immediately

        Returns: Animation object.

        .. See also::

            http://matplotlib.org/api/animation_api.html
            http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

        .. Note::

            It would be nice to animate the title of the plot, unfortunately
            this feature is not available in the present version of matplotlib.
            See: http://stackoverflow.com/questions/17558096/animated-title-in-matplotlib
        """
        ebands_list, edos_list = self.ebands_list, self.edoses_list
        if edos_list and len(edos_list) != len(ebands_list):
            raise ValueError("The number of objects for DOS must be equal to the number of bands")
        #titles = list(self.ebands_dict.keys())

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
                #if titles is not None: lines += [ax.set_title(titles[i])]
                artists.append(lines)
        else:
            # Animation with band structures + DOS.
            from matplotlib.gridspec import GridSpec
            gspec = GridSpec(nrows=1, ncols=2, width_ratios=width_ratios, wspace=0.05)
            ax0 = plt.subplot(gspec[0])
            ax1 = plt.subplot(gspec[1], sharey=ax0)
            ebands_list[0].decorate_ax(ax0)
            ax1.grid(True)
            ax1.yaxis.set_ticks_position("right")
            ax1.yaxis.set_label_position("right")

            for i, (ebands, edos) in enumerate(zip(ebands_list, edos_list)):
                # Define the zero of energy to align bands and dos
                mye0 = ebands.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
                ebands_lines = ebands.plot_ax(ax0, mye0, **plotax_kwargs)
                edos_lines = edos.plot_ax(ax1, mye0, exchange_xy=True, **plotax_kwargs)
                lines = ebands_lines + edos_lines
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

    def _repr_html_(self):
        """Integration with jupyter_ notebooks."""
        return self.ipw_select_plot()

    def ipw_select_plot(self): # pragma: no cover
        """
        Return an ipython widget with controllers to select the plot.
        """
        def plot_callback(plot_type, e0):
            r = getattr(self, plot_type)(e0=e0, show=True)
            if plot_type == "animate": return r

        import ipywidgets as ipw
        return ipw.interact_manual(
                plot_callback,
                plot_type=["combiplot", "gridplot", "boxplot", "combiboxplot", "animate"],
                e0=["fermie", "0.0"],
            )

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        # Use pickle file for data persistence.
        tmpfile = self.pickle_dump()

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("plotter = abilab.ElectronBandsPlotter.pickle_load('%s')" % tmpfile),
            nbv.new_code_cell("print(plotter)"),
            nbv.new_code_cell("frame = plotter.get_ebands_frame()\ndisplay(frame)"),
            nbv.new_code_cell("ylims = (None, None)"),
            nbv.new_code_cell("plotter.gridplot(ylims=ylims);"),
            nbv.new_code_cell("plotter.combiplot(ylims=ylims);"),
            nbv.new_code_cell("plotter.boxplot();"),
            nbv.new_code_cell("plotter.combiboxplot();"),
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
    This object reads band structure data from a netcdf_ file.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ElectronReader
    """
    def read_ebands(self):
        """
        Returns an instance of |ElectronBands|. Main entry point for client code
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
        """|numpy-array| with the number of bands indexed by [s, k]."""
        return self.read_value("number_of_states")

    def read_nspinor(self):
        """Number of spinors."""
        return self.read_dimvalue("number_of_spinor_components")

    def read_nsppol(self):
        """Number of independent spins (collinear case)."""
        return self.read_dimvalue("number_of_spins")

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
        scheme = self.read_string("smearing_scheme")

        return Smearing(
            scheme=scheme,
            occopt=occopt,
            tsmear_ev=units.Energy(self.read_value("smearing_width"), "Ha").to("eV")
        )


class ElectronDos(object):
    """
    This object stores the electronic density of states.
    It is usually created by calling the get_edos method of |ElectronBands|.
    """

    def __init__(self, mesh, spin_dos, nelect, fermie=None):
        """
        Args:
            mesh: array-like object with the mesh points in eV.
            spin_dos: array-like object with the DOS value for the different spins.
                      spin_dos[1, nw] if spin-unpolarized.
                      spin_dos[2, nw] if spin-polarized case.
            nelect: Number of electrons in the unit cell.
            fermie: Fermi level in eV. If None, fermie is obtained from the idos integral.

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
                print("tot_idos values:\n", self.tot_idos)
                raise

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app("nsppol: %d, nelect: %s" % (self.nsppol, self.nelect))
        app("Fermi energy: %s (eV) (recomputed from nelect):" % self.fermie)
        return "\n".join(lines)

    @classmethod
    def as_edos(cls, obj, edos_kwargs):
        """
        Return an instance of |ElectronDos| from a generic object ``obj``.
        Supports:

            - instances of cls
            - files (string) that can be open with abiopen and that provide an `ebands` attribute.
            - objects providing an `ebands` or `get_edos` attribute

        Args:
            edos_kwargs: optional dictionary with the options passed to `get_edos` to compute the electron DOS.
            Used when obj is not already an instance of ``cls``.
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
            # If the mesh is not large enough, we never cross nelect
            # If the last point in IDOS is sufficiently close to nelect
            # use it as Fermi level.
            if abs(idos.values[-1] - nelect) < 1e-3:
                i = len(idos) -1
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

    @lazy_property
    def up_minus_down(self):
        """
        Function1D with dos_up - dos_down
        """
        if self.nsppol == 1: # DOH!
            return Function1D.from_constant(self.spin_dos[0].mesh, 0.0)
        else:
            return self.spin_dos[0] - self.spin_dos[1]

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
                try:
                    return float(e0)
                except:
                    raise TypeError("Wrong value for e0: %s" % str(e0))
        else:
            # Assume number
            return float(e0)

    def plot_ax(self, ax, e0, spin=None, what="dos", fact=1.0, exchange_xy=False, **kwargs):
        """
        Helper function to plot the DOS data on the axis ``ax``.

        Args:
            ax: |matplotlib-Axes|.
            e0: Option used to define the zero of energy in the band structure plot.
            spin: selects the spin component, None for total DOS, IDOS.
            what: string selecting what will be plotted. "dos" for DOS, "idos" for IDOS
            fact: Multiplication factor for DOS/IDOS. Usually +-1 for spin DOS
            exchange_xy: True to exchange x-y axis.
            kwargs: Options passed to matplotlib ``ax.plot``

        Return: list of lines added to the axis ax.
        """
        dosf, idosf = self.dos_idos(spin=spin)
        e0 = self.get_e0(e0)

        w2f = {"dos": dosf, "idos": idosf}
        if what not in w2f:
            raise ValueError("Unknown value for what: `%s`" % str(what))
        f = w2f[what]

        xx, yy = f.mesh - e0, f.values * fact
        if exchange_xy: xx, yy = yy, xx
        lines = []
        lines.extend(ax.plot(xx, yy, **kwargs))

        return lines

    @add_fig_kwargs
    def plot(self, e0="fermie", spin=None, ax=None, xlims=None, **kwargs):
        """
        Plot electronic DOS

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                - Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                - None: Don't shift energies, equivalent to ``e0 = 0``.
            spin: Selects the spin component, None if total DOS is wanted.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            kwargs: options passed to ``ax.plot``.

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        e0 = self.get_e0(e0)

        for spin in range(self.nsppol):
            if spin == 0:
                opts = {"color": "black", "linewidth": 1.0}
            else:
                opts = {"color": "red", "linewidth": 1.0}
            opts.update(kwargs)
            spin_sign = +1 if spin == 0 else -1
            x, y = self.spin_dos[spin].mesh - e0, spin_sign * self.spin_dos[spin].values
            ax.plot(x, y, **opts)

        ax.grid(True)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('DOS (states/eV)')
        set_axlims(ax, xlims, "x")

        return fig

    @add_fig_kwargs
    def plot_dos_idos(self, e0="fermie", ax_list=None, xlims=None, height_ratios=(1, 2), **kwargs):
        """
        Plot electronic DOS and Integrated DOS on two different subplots.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``.
            ax_list: The axes for the DOS and IDOS plot. If ax_list is None, a new figure
                is created and the two axes are automatically generated.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            height_ratios:
            kwargs: options passed to ``plot_ax``

        Return: |matplotlib-Figure|
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        if ax_list is None:
            fig = plt.figure()
            gspec = GridSpec(nrows=2, ncols=1, height_ratios=height_ratios, wspace=0.05)
            ax0 = plt.subplot(gspec[0])
            ax1 = plt.subplot(gspec[1], sharex=ax0)
            ax_list = [ax0, ax1]

            for ax in ax_list:
                ax.grid(True)
                set_axlims(ax, xlims, "x")

            ax0.set_ylabel("TOT IDOS")
            ax1.set_ylabel("TOT DOS")
            ax1.set_xlabel('Energy (eV)')
        else:
            fig = ax_list[0].get_figure()

        for spin in range(self.nsppol):
            if spin == 0:
                opts = {"color": "black", "linewidth": 1.0}
            else:
                opts = {"color": "red", "linewidth": 1.0}
            # Plot Total dos if unpolarized.
            if self.nsppol == 1: spin = None
            self.plot_ax(ax_list[0], e0, spin=spin, what="idos", **opts)
            self.plot_ax(ax_list[1], e0, spin=spin, what="dos", **opts)

        return fig

    @add_fig_kwargs
    def plot_up_minus_down(self, e0="fermie", ax=None, xlims=None, **kwargs):
        """
        Plot Dos_up -Dow_down

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            kwargs: options passed to ``ax.plot``.

        Return: |matplotlib-Figure|
        """
        dos_diff = self.up_minus_down
        idos_diff = dos_diff.integral()

        e0 = self.get_e0(e0)
        if not kwargs:
            kwargs = {"color": "black", "linewidth": 1.0}

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.plot(dos_diff.mesh - e0, dos_diff.values, **kwargs)
        ax.plot(idos_diff.mesh - e0, idos_diff.values, **kwargs)

        ax.grid(True)
        set_axlims(ax, xlims, "x")
        ax.set_ylabel('Dos_up - Dos_down (states/eV)')
        ax.set_xlabel('Energy (eV)')

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
    # TODO: down-up option animate?

    def __init__(self, key_edos=None, edos_kwargs=None):
        if key_edos is None: key_edos = []
        key_edos = [(k, ElectronDos.as_edos(v, edos_kwargs)) for k, v in key_edos]
        self.edoses_dict = OrderedDict(key_edos)

    def __len__(self):
        return len(self.edoses_dict)

    @property
    def edos_list(self):
        """List of DOSes"""
        return list(self.edoses_dict.values())

    def add_edos(self, label, edos, edos_kwargs=None):
        """
        Adds a DOS for plotting.

        Args:
            label: label for the DOS. Must be unique.
            edos: |ElectronDos| object.
            edos_kwargs: optional dictionary with the options passed to ``get_edos`` to compute the electron DOS.
                Used only if ``edos`` is not an ElectronDos instance.
        """
        if label in self.edoses_dict:
            raise ValueError("label %s is already in %s" % (label, list(self.edoses_dict.keys())))
        self.edoses_dict[label] = ElectronDos.as_edos(edos, edos_kwargs)

    @add_fig_kwargs
    def combiplot(self, what_list="dos", spin_mode="total", e0="fermie",
                  ax_list=None,  xlims=None, fontsize=8, **kwargs):
        """
        Plot the the DOSes on the same figure. Use ``gridplot`` to plot DOSes on different figures.

        Args:
            what_list: Selects quantities to plot e.g. ["dos", "idos"] to plot DOS and integrated DOS.
                "dos" for DOS only and "idos" for IDOS only
            spin_mode: "total" for total (I)DOS, "resolved" for plotting individual contributions.
                Meaningful only if nsppol == 2.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize (int): fontsize for titles and legend

        Return: |matplotlib-Figure|
        """
        what_list = list_strings(what_list)
        nrows, ncols = len(what_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        can_use_basename = self._can_use_basenames_as_labels()
        for i, (what, ax) in enumerate(zip(what_list, ax_list)):
            for label, edos in self.edoses_dict.items():
                if can_use_basename:
                    label = os.path.basename(label)
                else:
                    # Use relative paths if label is a file.
                    if os.path.isfile(label): label = os.path.relpath(label)

                # Here I handle spin and spin_mode.
                if edos.nsppol == 1 or spin_mode == "total":
                    # Plot total values
                    edos.plot_ax(ax, e0, what=what, spin=None, label=label)

                elif spin_mode == "resolved":
                    # Plot spin resolved quantiies with sign.
                    # Note get_color to have same color for both spins.
                    for spin in range(edos.nsppol):
                        fact = 1 if spin == 0 else -1
                        lines = edos.plot_ax(ax, e0, what=what, spin=spin, fact=fact,
                            color=None if spin == 0 else lines[0].get_color(),
                            label=label if spin == 0 else None)
                else:
                    raise ValueError("Wrong value for spin_mode: `%s`:" % str(spin_mode))

            ax.grid(True)
            if i == len(what_list) - 1:
                ax.set_xlabel("Energy (eV)")
            ax.set_ylabel('DOS (states/eV)' if what == "dos" else "IDOS")
            set_axlims(ax, xlims, "x")
            ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    # An alias for combiplot.
    plot = combiplot

    @add_fig_kwargs
    def gridplot(self, what="dos", spin_mode="total", e0="fermie",
                 sharex=True, sharey=True, xlims=None, fontsize=8, **kwargs):
        """
        Plot multiple DOSes on a grid.

        Args:
            what: "dos" to plot DOS, "idos" for integrated DOS.
            spin_mode: "total" for total (I)DOS, "resolved" for plotting individual contributions.
                Meaningful only if nsppol == 2.
            e0: Option used to define the zero of energy in the band structure plot. Possible values::

                - ``fermie``: shift all eigenvalues and the DOS to have zero energy at the Fermi energy.
                   Note that, by default, the Fermi energy is taken from the band structure object
                   i.e. the Fermi energy computed at the end of the SCF file that produced the density.
                   This should be ok in semiconductors. In metals, however, a better value of the Fermi energy
                   can be obtained from the DOS provided that the k-sampling for the DOS is much denser than
                   the one used to compute the density. See ``edos_fermie``.
                - ``edos_fermie``: Use the Fermi energy computed from the DOS to define the zero of energy in both subplots.
                   Available only if edos_objects is not None
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``.

            sharex, sharey: True if x (y) axis should be shared.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Label and title fontsize.

        Return: |matplotlib-Figure|
        """
        titles = list(self.edoses_dict.keys())
        edos_list = self.edos_list

        nrows, ncols = 1, 1
        numeb = len(edos_list)
        if numeb > 1:
            ncols = 2
            nrows = numeb // ncols + numeb % ncols

        # Build Grid
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=sharex, sharey=sharey, squeeze=False)
        ax_list = ax_list.ravel()

        # don't show the last ax if numeb is odd.
        if numeb % ncols != 0: ax_list[-1].axis("off")

        for i, ((label, edos), ax) in enumerate(zip(self.edoses_dict.items(), ax_list)):
            irow, icol = divmod(i, ncols)

            # Here I handle spin and spin_mode.
            if edos.nsppol == 1 or spin_mode == "total":
                opts = {"color": "black", "linewidth": 1.0}
                edos.plot_ax(ax, e0=e0, what=what, spin=None, **opts)

            elif spin_mode == "resolved":
                # Plot spin resolved quantiies with sign.
                # Note get_color to have same color for both spins.
                for spin in range(edos.nsppol):
                    fact = 1 if spin == 0 else -1
                    lines = edos.plot_ax(ax, e0, what=what, spin=spin, fact=fact,
                        color=None if spin == 0 else lines[0].get_color(),
                        label=label if spin == 0 else None)
            else:
                raise ValueError("Wrong value for spin_mode: `%s`:" % str(spin_mode))

            ax.grid(True)
            ax.set_title(label, fontsize=fontsize)
            set_axlims(ax, xlims, "x")
            if (irow, icol) == (0, 0):
                ax.set_ylabel('DOS (states/eV)' if what == "dos" else "IDOS")
            if irow == nrows - 1:
                ax.set_xlabel("Energy (eV)")

            #ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    def ipw_select_plot(self): # pragma: no cover
        """
        Return an ipython widget with controllers to select the plot.
        """
        def plot_callback(plot_type, e0):
            getattr(self, plot_type)(e0=e0, show=True)

        import ipywidgets as ipw
        return ipw.interact_manual(
                plot_callback,
                plot_type=["combiplot", "gridplot"],
                e0=["fermie", "0.0"],
            )

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.combiplot(show=False)
        yield self.gridplot(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If nbpath is None, a temporay file in the current
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
            nbv.new_code_cell("plotter.combiplot(xlims=xlims);"),
            nbv.new_code_cell("plotter.gridplot(xlims=xlims);"),
        ])

        return self._write_nb_nbpath(nb, nbpath)

    def _can_use_basenames_as_labels(self):
        """
        Return True if all labels represent valid files and the basenames are unique
        In this case one can use the file basename instead of the full path in the plots.
        """
        if not all(os.path.exists(l) for l in self.edoses_dict): return False
        labels = [os.path.basename(l) for l in self.edoses_dict]
        return len(set(labels)) == len(labels)


class Bands3D(Has_Structure):

    def __init__(self, structure, ibz, has_timrev, eigens, fermie):
        """
        This object reconstructs by symmetry the eigenvalues in the full BZ starting from the IBZ.
        Provides methods to extract and visualize isosurfaces.

        Args:
            structure:
            ibz:
            has_timrev:
            eigens:
            fermie
        """
        self.ibz = ibz
        self._structure = structure
        self.reciprocal_lattice = structure.lattice.reciprocal_lattice
        self.has_timrev = has_timrev
        self.fermie = fermie
        self.eigens = np.atleast_3d(eigens)
        self.nsppol, _, self.nband = self.eigens.shape

        # Sanity check.
        errors = []; eapp = errors.append
        if not self.ibz.is_ibz:
            eapp("Expecting an IBZ sampling but got %s" % type(self.ibz))
        if not self.ibz.is_mpmesh:
            eapp("Monkhorst-Pack meshes are required.\nksampling: %s" % str(self.ibz.ksampling))

        mpdivs, shifts = self.ibz.mpdivs_shifts
        if shifts is not None and not np.all(shifts == 0.0):
            eapp("Gamma-centered k-meshes are required by Xcrysden.")
        if errors:
            raise ValueError("\n".join(errors))

        # Xcrysden requires points in the unit cell (C-order)
        # and the mesh must include the periodic images hence pbc=True.
        self.uc2ibz = map_grid2ibz(self.structure, self.ibz.frac_coords, mpdivs, self.has_timrev, pbc=True)
        self.mpdivs = mpdivs
        self.kdivs = mpdivs + 1
        self.spacing = 1.0 / mpdivs
        self.ucdata_shape = (self.nsppol, self.nband) + tuple(self.kdivs)
        #self.ibzdata_shape = (self.nsppol, self.nband, len(self.ibz))

        # Construct energy bands on unit cell grid: e_{TSk} = e_{k}
        self.ucdata_sbk = self.symmetrize_ibz_scalars(self.eigens)

        self.ucell_scalars = OrderedDict()
        self.ucell_vectors = OrderedDict()
        #if reference_sb is None:
        #self.reference_sb = [[] for _ in self.spins]
        #for spin in self.spins:
        #    self.reference_sb[spin] = {band: band for band in self.bands}
        #else:

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    # Handy variables used to loop
    @property
    def spins(self):
        return range(self.nsppol)

    @property
    def bands(self):
        return range(self.nband)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []
        app = lines.append
        # TODO: Finalize implementation
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")

        return "\n".join(lines)

    def add_ucell_scalars(self, name, scalars):
        """
        Add scalar quantities given in the unit cell.

        Args:
            name: keyword used to store scalars.
            scalars:
        """
        self.ucell_scalars[name] = np.reshape(scalars, self.ucdata_shape)

    def add_ibz_scalars(self, name, scalars, inshape="skb"):
        """
        Add scalar quantities given in the IBZ i.e. symmetrize values to get array in unit cell.

        Args:
            name: keyword used to store symmetrized values.
            scalars: scalars in IBZ. See ``inshape`` for shape
            inshape: shape of input scalars. "skb" if (nsppol, nkibz, nband)
            "sbk" for (nsppol, nband, nkibz).
        """
        self.add_ucell_scalars(name, self.symmetrize_ibz_scalars(scalars, inshape=inshape))

    def symmetrize_ibz_scalars(self, scalars, inshape="skb"):
        """
        Symmetrize scalar quantities given in the IBZ.

        Args:
            scalars: scalars in IBZ. See `inshape` for shape
            inshape: shape of input scalars. "skb" if (nsppol, nkibz, nband)
            "sbk" for (nsppol, nband, nkibz).

        Return:
            |numpy-array| with scalars in unit cell. shape is **always**: (nsppol, nband, nkbz)
        """
        # Symmetrize scalars unit cell grid: e_{TSk} = e_{k}
        ucdata_sbk = np.empty((self.nsppol, self.nband, len(self.uc2ibz)))

        if inshape == "skb":
            scalars = np.reshape(scalars, (self.nsppol, len(self.ibz), self.nband))
            for ikuc, ik_ibz in enumerate(self.uc2ibz):
                ucdata_sbk[:, :, ikuc] = scalars[:, ik_ibz, :]
        elif inshape == "sbk":
            scalars = np.reshape(scalars, (self.nsppol, self.nband, len(self.ibz)))
            for ikuc, ik_ibz in enumerate(self.uc2ibz):
                ucdata_sbk[:, :, ikuc] = scalars[:, :, ik_ibz]
        else:
            raise ValueError("Wrong inshape: %s" % str(insp))

        return ucdata_sbk

    #def add_ucell_vectors(self, name, vectors, inshape="skb"):
    #    self.ucell_vectors[name] = np.reshape(vectors, self.ucdata + (3,))

    #def add_ibz_vectors(self, name, scalars, inshape="skb")
    #    self.add_ucell_vectors(name, self.symmetrize_ibz_vectors(vectors, inshape=inshape))

    #def wsmap(self):
        #ws = -np.ones(ngkpt, dtype=np.int)
        #for i in range(ngkpt[0]):
        #    ki = (i - ngkpt[0] // 2)
        #    if ki < 0: ki += ngkpt[0]
        #    for j in range(ngkpt[1]):
        #        kj = (j - ngkpt[1] // 2)
        #        if kj < 0: kj += ngkpt[1]
        #        for k in range(ngkpt[2]):
        #            kz = (k - ngkpt[2] // 2)
        #            if kz < 0: kz += ngkpt[2]
        #            #bzgrid2ibz[gp_bz[0], gp_bz[1], gp_bz[2]] = ik_ibz
        #            ws[i, j, k] = bzgrid2ibz[ki, kj, kz]
        #bzgrid2ibz = ws

    def get_isobands(self, e0):
        """Return index of the bands crossing ``e0``in eV. None if no band is found."""
        isobands = [[] for _ in self.spins]
        for spin in self.spins:
            for band in self.bands:
                emin, emax = self.eigens[spin, :, band].min(), self.eigens[spin, :, band].max()
                if isobands[spin] and e0 > emax: break
                if emax >= e0 >= emin: isobands[spin].append(band)
        if all(not l for l in isobands): return None
        return isobands

    def xcrysden_view(self):  # pragma: no cover
        """
        Visualize electron energy isosurfaces with xcrysden_.
        """
        _, tmp_filepath = tempfile.mkstemp(suffix=".bxsf", text=True)
        #print("Producing BXSF file in:", tmp_filepath)
        self.to_bxsf(tmp_filepath, unit="eV")
        from abipy.iotools.visualizer import Xcrysden
        return Xcrysden(tmp_filepath)()

    def to_bxsf(self, filepath, unit="eV"):
        """
        Export the full band structure to ``filepath`` in BXSF format
        suitable for the visualization of the Fermi surface with xcrysden_ (use ``xcrysden --bxsf FILE``).
        Require k-points in IBZ and gamma-centered k-mesh.

        Args:
            filepath: BXSF filename or stream.
            unit: Input energies are in unit ``unit``.
        """
        from abipy.iotools import bxsf_write
        if hasattr(filepath, "write"):
            return bxsf_write(filepath, self.structure, self.nsppol, self.nband, self.kdivs,
                              self.ucdata_sbk, self.fermie, unit=unit)
        else:
            with open(filepath, "wt") as fh:
                bxsf_write(fh, self.structure, self.nsppol, self.nband, self.kdivs,
                           self.ucdata_sbk, self.fermie, unit=unit)
                return filepath

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

    @add_fig_kwargs
    def plot_isosurfaces(self, e0="fermie", verbose=0, **kwargs):
        """
        Plot isosurface with matplotlib_

        .. warning::

            Requires scikit-image package, matplotlib rendering is usually slow.

        Args:
            e0: Isolevel in eV. Default: Fermi energy.
            verbose: verbosity level.

        Return: |matplotlib-Figure|
        """
        try:
            import skimage
        except ImportError:
            raise ImportError("scikit-image not installed.\n"
                "Please install with it with `conda install scikit-image` or `pip install scikit-image`")

        try:
            from skimage.measure import marching_cubes_lewiner as marching_cubes
        except ImportError:
            from skimage.measure import marching_cubes

        e0 = self.get_e0(e0)
        isobands = self.get_isobands(e0)
        if isobands is None: return None
        if verbose: print("Bands for isosurface:", isobands)

        from pymatgen.electronic_structure.plotter import plot_lattice_vectors, plot_wigner_seitz
        ax, fig, plt = get_ax3d_fig_plt(ax=None)
        plot_unit_cell(self.reciprocal_lattice, ax=ax, color="k", linewidth=1)
        #plot_wigner_seitz(self.reciprocal_lattice, ax=ax, color="k", linewidth=1)

        for spin in self.spins:
            for band in isobands[spin]:
                # From http://scikit-image.org/docs/stable/api/skimage.measure.html#marching-cubes
                # verts: (V, 3) array
                #   Spatial coordinates for V unique mesh vertices. Coordinate order matches input volume (M, N, P).
                # faces: (F, 3) array
                #   Define triangular faces via referencing vertex indices from verts.
                #   This algorithm specifically outputs triangles, so each face has exactly three indices.
                # normals: (V, 3) array
                #   The normal direction at each vertex, as calculated from the data.
                # values: (V, ) array
                #   Gives a measure for the maximum value of the data in the local region near each vertex.
                #   This can be used by visualization tools to apply a colormap to the mesh
                voldata = np.reshape(self.ucdata_sbk[spin, band], self.kdivs)
                verts, faces, normals, values = marching_cubes(voldata, level=e0, spacing=tuple(self.spacing))
                #verts, faces, normals, values = marching_cubes_lewiner(voldata, level=e0, spacing=tuple(self.spacing))
                verts = self.reciprocal_lattice.get_cartesian_coords(verts)
                ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2]) #, cmap='Spectral', lw=1, antialiased=True)
                # mayavi package:
                #mlab.triangular_mesh([v[0] for v in verts], [v[1] for v in verts], [v[2] for v in verts], faces) #, color=(0, 0, 0))

        ax.set_axis_off()

        return fig

    def mvplot_isosurfaces(self, e0="fermie", verbose=0, figure=None, show=True):  # pragma: no cover
        """
        Plot isosurface with mayavi_

        Args:
            e0:
            verbose:
            show:
        """
        # Find bands crossing e0.
        e0 = self.get_e0(e0)
        isobands = self.get_isobands(e0)
        if isobands is None: return None
        if verbose: print("Bands for isosurface:", isobands)

        #from pymatgen.electronic_structure.plotter import plot_fermi_surface
        #spin, band = 0, 4
        #for i, band in enumerate(isobands[spin]):
        #data = np.reshape(isoenes[0][band], mpdivs + 1 if pbc else mpdivs)
        #plot_fermi_surface(data, self.structure, False, energy_levels=[e0]) interative=not (i == len(isobands[spin]) - 1))

        # Plot isosurface with mayavi.
        from abipy.display import mvtk
        figure, mlab = mvtk.get_fig_mlab(figure=figure)
        mvtk.plot_unit_cell(self.reciprocal_lattice, figure=figure)
        mvtk.plot_wigner_seitz(self.reciprocal_lattice, figure=figure)
        cell = self.reciprocal_lattice.matrix

        for spin in self.spins:
            for band in isobands[spin]:
                data = np.reshape(self.ucdata_sbk[spin, band], self.kdivs)
                cp = mlab.contour3d(data, contours=[e0], transparent=True,
                                    #colormap='hot', color=(0, 0, 1), opacity=1.0, figure=figure)
                                    colormap='Set3', opacity=0.9, figure=figure)

                polydata = cp.actor.actors[0].mapper.input
                pts = np.array(polydata.points) #  - 1  # TODO this + mpdivs should be correct
                if verbose: print("shape:", pts.shape, pts)
                polydata.points = np.dot(pts, cell / np.array(data.shape)[:, np.newaxis])
                #polydata.points = np.dot(pts, cell / np.array(self.mpdivs)[:, np.newaxis])
                mlab.view(distance="auto", figure=figure)

        # Add k-point labels.
        labels = {k.name: k.frac_coords for k in self.structure.hsym_kpoints}
        mvtk.plot_labels(labels, lattice=self.structure.reciprocal_lattice, figure=figure)

        if show: mlab.show()
        return figure

    @add_fig_kwargs
    def plot_contour(self, band, spin=0, plane="xy", elevation=0, ax=None, fontsize=8, **kwargs):
        """
        Contour plot with matplotlib_.

        Args:
            band: Band index
            spin: Spin index.
            plane:
            elevation:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Label and title fontsize.

        Return: |matplotlib-Figure|
        """
        data = np.reshape(self.ucdata_sbk[spin, band], self.kdivs) - self.fermie

        x = np.arange(self.kdivs[0]) / (self.kdivs[0] - 1)
        y = np.arange(self.kdivs[1]) / (self.kdivs[1] - 1)
        fxy = data[:, :, elevation]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        x, y = np.meshgrid(x, y)
        c = ax.contour(x, y, fxy, **kwargs)
        ax.clabel(c, inline=1, fontsize=fontsize)
        kvert = dict(xy="z", xz="y", yz="x")[plane]
        ax.set_title(r"Band %s in %s plane at $K_{%s}=%d$" % (band, plane, kvert, elevation), fontsize=fontsize)
        ax.grid(True)
        ax.set_xlabel("$K_%s$" % plane[0])
        ax.set_ylabel("$K_%s$" % plane[1])

        return fig

    #def interpolate(self, densify_mpdivs):
        #densify_mpdivs = np.array(densify_mpdivs)
        #if np.any(densify_mpdivs > 1):
        #dense_mpdivs = densify_mpdivs * mpdivs
        #dense_kpts = kmesh_from_mpdivs(dense_mpdivs, shifts=(0, 0, 0), pbc=pbc, order="unit_cell")
        #from scipy.interpolate import RegularGridInterpolator
        #x = np.arange(0, mpdivs[0] + 1) / mpdivs[0]
        #y = np.arange(0, mpdivs[1] + 1) / mpdivs[1]
        #z = np.arange(0, mpdivs[2] + 1) / mpdivs[2]
        #for spin in self.spins:
        #    for band in isobands[spin]:
        #        interp = RegularGridInterpolator((x, y, z), isoenes[spin][band], method='linear')
        #        isoenes[spin][band] = interp(dense_mpdivs)

    #def mvplot_surf(self):
        #spin, band = 0, 2
        #data = np.reshape(isoenes[spin][band], mpdivs + 1 if pbc else mpdivs)
        #x, y = np.arange(data.shape[0]), np.arange(data.shape[1])
        #cp = mlab.surf(x, y, data[0,:,:], figure=figure)
        #polydata = cp.actor.actors[0].mapper.input
        #pts = np.array(polydata.points) # - 1
        #xs, ys, zs = pts.T
        #print(pts.shape)
        #print(zs)
        #polydata.points = np.dot(pts, cell / np.array(data.shape)[:, np.newaxis])

        #data = np.reshape(isoenes[spin][band+1], mpdivs + 1 if pbc else mpdivs)
        #cp = mlab.surf(x, y, data[0,:,:], figure=figure)
        #polydata = cp.actor.actors[0].mapper.input
        #pts = np.array(polydata.points) # - 1
        #polydata.points = np.dot(pts, cell / np.array(data.shape)[:, np.newaxis])
        #mlab.view(distance="auto", figure=figure)
        #if show: mlab.show()
        #return

    def mvplot_cutplanes(self, band, spin=0, figure=None, show=True, **kwargs): # pragma: no cover
        """Plot cutplanes with mayavi_."""
        data = np.reshape(self.ucdata_sbk[spin, band], self.kdivs) - self.fermie
        contours = [-1.0, 0.0, 1.0]

        from abipy.display import mvtk
        figure, mlab = mvtk.get_fig_mlab(figure=figure)
        src = mlab.pipeline.scalar_field(data)

        mlab.pipeline.image_plane_widget(src, plane_orientation='x_axes', slice_index=self.kdivs[0]//2)
        mlab.pipeline.image_plane_widget(src, plane_orientation='y_axes', slice_index=self.kdivs[1]//2)
        mlab.pipeline.image_plane_widget(src, plane_orientation='z_axes', slice_index=self.kdivs[2]//2)
        mlab.pipeline.iso_surface(src, contours=contours) #, opacity=0.1)
        #mlab.pipeline.iso_surface(src, contours=[data.min()+ 0.1 * data.ptp()], opacity=0.1)
        mlab.outline()

        if show: mlab.show()
        return figure

    #def write_data(self, workdir, fmt="cube", rmdir=False)


class ElectronBands3D(Bands3D):
    pass

#class PhononBands3D(Bands3D):
#    pass


class RobotWithEbands(object):
    """
    Mixin class for robots associated to files with |ElectronBands|.
    """
    def combiplot_ebands(self, **kwargs):
        """Wraps combiplot method of |ElectronBandsPlotter|. kwargs passed to combiplot."""
        return self.get_ebands_plotter().combiplot(**kwargs)

    def gridplot_ebands(self, **kwargs):
        """Wraps gridplot method of |ElectronBandsPlotter|. kwargs passed to gridplot."""
        return self.get_ebands_plotter().gridplot(**kwargs)

    def boxplot_ebands(self, **kwargs):
        """Wraps boxplot method of |ElectronBandsPlotter|. kwargs passed to boxplot."""
        return self.get_ebands_plotter().boxplot(**kwargs)

    def combiboxplot_ebands(self, **kwargs):
        """Wraps combiboxplot method of |ElectronDosPlotter|. kwargs passed to combiboxplot."""
        return self.get_ebands_plotter().combiboxplot(**kwargs)

    def combiplot_edos(self, **kwargs):
        """Wraps combiplot method of |ElectronDosPlotter|. kwargs passed to combiplot."""
        return self.get_edos_plotter().combiplot(**kwargs)

    def gridplot_edos(self, **kwargs):
        """Wraps gridplot method of |ElectronDosPlotter|. kwargs passed to gridplot."""
        return self.get_edos_plotter().gridplot(**kwargs)

    def get_ebands_plotter(self, filter_abifile=None, cls=None):
        """
        Build and return an instance of |ElectronBandsPlotter| or a subclass is ``cls`` is not None.

        Args:
            filter_abifile: Function that receives an ``abifile`` object and returns
                True if the file should be added to the plotter.
            cls: subclass of |ElectronBandsPlotter|.
        """
        plotter = ElectronBandsPlotter() if cls is None else cls()

        for label, abifile in self.items():
            if filter_abifile is not None and not filter_abifile(abifile): continue
            plotter.add_ebands(label, abifile.ebands)

        return plotter

    def get_edos_plotter(self, cls=None, filter_abifile=None, **kwargs):
        """
        Build and return an instance of |ElectronDosPlotter| or a subclass is cls is not None.

        Args:
            filter_abifile: Function that receives an ``abifile` object and returns
                True if the file should be added to the plotter.
            cls: subclass of |ElectronDosPlotter|.
            kwargs: Arguments passed to ebands.get_edos
        """
        plotter = ElectronDosPlotter() if cls is None else cls()

        for label, abifile in self.items():
            if filter_abifile is not None and not filter_abifile(abifile): continue
            if not abifile.ebands.kpoints.is_ibz:
                cprint("Skipping %s because kpoint sampling not IBZ" % abifile.filepath, "magenta")
                continue
            plotter.add_edos(label, abifile.ebands.get_edos(**kwargs))

        return plotter

    #def get_ebands_dataframe(self, with_spglib=True):
    #    return dataframe_from_ebands(self.ncfiles, index=list(self.keys()), with_spglib=with_spglib)

    @add_fig_kwargs
    def plot_egaps(self, sortby=None, hue=None, fontsize=6, **kwargs):
        """
        Plot the convergence of the direct and fundamental gaps
        wrt to the ``sortby`` parameter. Values can optionally be grouped by ``hue``.

        Args:
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        # Note: Handling nsppol > 1 and the case in which we have abifiles with different nsppol is a bit tricky
        # hence we have to handle the different cases explicitly (see get_xy)
        if not self.abifiles: return None
        max_nsppol = max(f.nsppol for f in self.abifiles)

        items = ["fundamental_gaps", "direct_gaps", "bandwidths"]

        def get_xy(item, spin, all_xvals, all_abifiles):
            """
            Extract (xvals, yvals) from all_abifiles for given (item, spin) and initial all_xvals.
            Here we handle the case in which we have files with different nsppol.
            """
            xvals, yvals = [], []

            for i, af in enumerate(all_abifiles):
                if spin > af.nsppol - 1: continue
                xvals.append(all_xvals[i])
                if callable(item):
                    yy = float(item(af.ebands))
                else:
                    yy = getattr(af.ebands, item)
                    if item in ("fundamental_gaps", "direct_gaps"):
                        yy = yy[spin].energy
                    else:
                        yy = yy[spin]

                yvals.append(yy)

            return xvals, yvals

        # Build grid plot.
        nrows, ncols = len(items), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # Sort and group files if hue.
        if hue is None:
            labels, ncfiles, params = self.sortby(sortby, unpack=True)
        else:
            groups = self.group_and_sortby(hue, sortby)

        marker_spin = {0: "^", 1: "v"}
        for i, (ax, item) in enumerate(zip(ax_list, items)):
            for spin in range(max_nsppol):
                if hue is None:
                    # Extract data.
                    xvals, yvals = get_xy(item, spin, params, self.abifiles)
                    if not is_string(xvals[0]):
                        ax.plot(xvals, yvals, marker=marker_spin[spin], **kwargs)
                    else:
                        # Must handle list of strings in a different way.
                        xn = range(len(xvals))
                        ax.plot(xn, yvals, marker=marker_spin[spin], **kwargs)
                        ax.set_xticks(xn)
                        ax.set_xticklabels(xvals, fontsize=fontsize)
                else:
                    for g in groups:
                        # Extract data.
                        xvals, yvals = get_xy(item, spin, g.xvalues, g.abifiles)
                        label = "%s: %s" % (self._get_label(hue), g.hvalue)
                        ax.plot(xvals, yvals, label=label, marker=marker_spin[spin], **kwargs)

            ax.grid(True)
            ax.set_ylabel(self._get_label(item))
            if i == len(items) - 1:
                ax.set_xlabel("%s" % self._get_label(sortby))
                if sortby is None: rotate_ticklabels(ax, 15)
            if i == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)
                #ax.legend(loc='best', fontsize=fontsize, shadow=True, fancybox=True, framealpha=0.5)

        return fig

    def get_ebands_code_cells(self, title=None):
        """Return list of notebook cells."""
        nbformat, nbv = self.get_nbformat_nbv()
        title = "## Code to compare multiple ElectronBands objects" if title is None else str(title)
        # Try not pollute namespace with lots of variables.
        return [
            nbv.new_markdown_cell(title),
            nbv.new_code_cell("robot.get_ebands_plotter().ipw_select_plot();"),
            nbv.new_code_cell("robot.get_edos_plotter().ipw_select_plot();"),
            nbv.new_code_cell("#robot.plot_egaps(sorby=None, hue=None);"),
        ]

    @add_fig_kwargs
    def gridplot_with_hue(self, hue, ylims=None, fontsize=8, sharex=False, sharey=False, **kwargs):
        """
        Plot multiple electron bandstructures on a grid. Group bands by ``hue``.

        Example:

            robot.gridplot_with_hue("nkpt")

        Args:
            hue: Variable that define subsets of the phonon bands, which will be drawn on separate plots.
                Accepts callable or string
                If string, it's assumed that `abifile has an attribute with the same name and getattr is invoked.
                Dot notation is also supported e.g. hue="structure.formula" --> abifile.structure.formula
                If callable, the output of hue(abifile) is used.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                or scalar e.g. `left`. If left (right) is None, default values are used
            fontsize: legend and title fontsize.
            sharex, sharey: True if X and Y axes should be shared.

        Returns: |matplotlib-Figure|
        """
        # Group abifiles by hue.
        groups = self.group_and_sortby(hue, func_or_string=None)
        nrows, ncols = len(groups), 1

        # Plot grid with phonon bands only.
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=sharex, sharey=sharey, squeeze=False)
        ax_list = ax_list.ravel()
        e0 = "fermie"  # Each ebands is aligned with respect to its Fermi energy.

        for ax, grp in zip(ax_list, groups):
            ax.grid(True)
            ebands_list = [abifile.ebands for abifile in grp.abifiles]
            ax.set_title("%s = %s" % (self._get_label(hue), grp.hvalue), fontsize=fontsize)

            nkpt_list = [ebands.nkpt for ebands in ebands_list]
            if any(nk != nkpt_list[0] for nk in nkpt_list):
                cprint("WARNING: Bands have different number of k-points:\n%s" % str(nkpt_list), "yellow")

            for i, (ebands, lineopts) in enumerate(zip(ebands_list, self.iter_lineopt())):
                # Plot all branches with lineopts and set the label of the last line produced.
                ebands.plot_ax(ax, e0, **lineopts)
                ax.lines[-1].set_label("%s" % grp.labels[i])

                # Set ticks and labels
                # (NB: we do this only for the first ebands, in principle ebands
                # in the group could have different k-points but there's need to be so strict here.
                if i == 0:
                    ebands.decorate_ax(ax, klabels=None)

            # Set legends.
            ax.legend(loc='best', fontsize=fontsize, shadow=True)
            set_axlims(ax, ylims, "y")

        return fig