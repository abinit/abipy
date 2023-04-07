# coding: utf-8
"""Classes to analyse electronic structures."""
from __future__ import annotations

import os
import copy
import itertools
import json
import warnings
import tempfile
import pickle
import numpy as np
import pandas as pd
import pymatgen.core.units as units
import abipy.core.abinit_units as abu

from collections import OrderedDict, namedtuple
from collections.abc import Iterable
from typing import List, Any
from monty.string import is_string, list_strings, marquee
from monty.termcolor import cprint
from monty.json import MontyEncoder
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
from monty.bisect import find_le, find_gt
from pymatgen.electronic_structure.core import Spin as PmgSpin
from abipy.tools.serialization import pmg_serialize
from abipy.core.func1d import Function1D
from abipy.core.mixins import Has_Structure, NotebookWriter
from abipy.core.kpoints import (Kpoint, KpointList, Kpath, IrredZone, KSamplingInfo, KpointsReaderMixin,
    Ktables, has_timrev_from_kptopt, map_grid2ibz) #, kmesh_from_mpdivs)
from abipy.core.structure import Structure
from abipy.iotools import ETSF_Reader
from abipy.tools import duck
from abipy.tools.numtools import gaussian
from abipy.tools.decorators import memoized_method
from abipy.tools.context_managers import Timer
from abipy.tools.plotting import (set_axlims, add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt,
    get_ax3d_fig_plt, rotate_ticklabels, set_visible, plot_unit_cell, set_ax_xylabels, get_figs_plotly,
    get_fig_plotly, add_plotly_fig_kwargs, PlotlyRowColDesc, plotly_klabels, plotly_set_lims)


__all__ = [
    "ElectronBands",
    "ElectronDos",
    "dataframe_from_ebands",
    "ElectronBandsPlotter",
    "ElectronDosPlotter",
]


SUBSCRIPT_UNICODE = {
                "0": "₀",
                "1": "₁",
                "2": "₂",
                "3": "₃",
                "4": "₄",
                "5": "₅",
                "6": "₆",
                "7": "₇",
                "8": "₈",
                "9": "₉",
            }


class Electron(namedtuple("Electron", "spin kpoint band eig occ kidx")):
    """
    Single-particle state.

    .. Attributes:

        spin: spin index (C convention, i.e >= 0)
        kpoint: |Kpoint| object.
        band: band index. (C convention, i.e >= 0)
        eig: KS eigenvalue.
        occ: Occupation factor.
        kidx: Index of the k-point in the initial array.

    .. note::

        Energies are in eV.
    """
    def __eq__(self, other):
        if other is None: return False
        return self.spin == other.spin and self.kpoint == other.kpoint and self.band == other.band

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        return "spin: %d, kpt: %s, band: %d, eig: %.3f, occ: %.3f" % (
            self.spin, self.kpoint.to_string(verbose=verbose), self.band, self.eig, self.occ)

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
        return super()._asdict()

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
    def __init__(self, in_state, out_state, all_kinds=None):
        """
        Args:
            in_state, out_state: Initial and finale state (:class:`Electron` instances).
            all_kinds: List of tuple. Each tuple gives the index of the k-point of the (initial, final) state.
                Used to plot e.g. all the optical gaps when there are equivalent k-points along the path.
        """
        self.in_state = in_state
        self.out_state = out_state
        if all_kinds is None:
            self.all_kinds = [(self.in_state.kidx, self.out_state.kidx)]
        else:
            # Provide default.
            self.all_kinds = all_kinds

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app("Energy: %.3f (eV)" % self.energy)
        app("Initial state: %s" % self.in_state.to_string(verbose=verbose))
        app("Final state:   %s" % self.out_state.to_string(verbose=verbose))

        return "\n".join(lines)

    def __eq__(self, other):
        if other is None: return False
        return self.in_state == other.in_state and self.out_state == other.out_state

    def __ne__(self, other):
        return not (self == other)

    @lazy_property
    def energy(self):
        """Transition energy in eV."""
        return self.out_state.eig - self.in_state.eig

    @lazy_property
    def qpoint(self):
        """k_final - k_initial"""
        return self.out_state.kpoint - self.in_state.kpoint

    @lazy_property
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
        super().__init__(*args, **kwargs)
        for mkey in self._MANDATORY_KEYS:
            if mkey not in self:
                raise ValueError("Mandatory key %s must be provided" % str(mkey))

    def __str__(self):
        return "smearing scheme: %s (occopt %d), tsmear_eV: %.3f, tsmear Kelvin: %.1f" % (
                    self.scheme, self.occopt, self.tsmear_ev, self.tsmear_ev / abu.kb_eVK )

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
    Object storing the electron band structure.

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
    def from_file(cls, filepath: str) -> ElectronBands:
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
    def from_binary_string(cls, bstring) -> ElectronBands:
        """
        Build object from a binary string with the netcdf data.
        Useful for implementing GUIs in which widgets returns binary data.
        """
        workdir = tempfile.mkdtemp()
        fd, tmp_path = tempfile.mkstemp(suffix=".nc")
        with open(tmp_path, "wb") as fh:
            fh.write(bstring)
            return cls.from_file(tmp_path)

    @classmethod
    def from_dict(cls, d: dict) -> ElectronBands:
        """Reconstruct object from the dictionary in MSONable format produced by as_dict."""
        d = d.copy()
        kd = d["kpoints"].copy()
        kd.pop("@module")

        kpoints_cls = KpointList.subclass_from_name(kd.pop("@class"))
        kpoints = kpoints_cls.from_dict(kd)

        # Needed to support old dictionaries
        if "nspden" not in d: d["nspden"] = 1
        if "nspinor" not in d: d["nspinor"] = 1

        smearing = d['smearing']
        if smearing is not None:
            # Handle metallic occupation scheme
            smearing = Smearing.from_dict(smearing)

        return cls(Structure.from_dict(d["structure"]), kpoints,
                   d["eigens"], d["fermie"], d["occfacts"], d["nelect"], d["nspinor"], d["nspden"],
                   nband_sk=d["nband_sk"], smearing=smearing,
                   linewidths=d.get("linewidths", None)
                   )

    @pmg_serialize
    def as_dict(self) -> dict:
        """Return dictionary with JSON serialization."""
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
    def as_ebands(cls, obj: Any) -> ElectronBands:
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

            if obj.endswith("_EBANDS.nc"):
                return cls.from_file(obj)

            from abipy.abilab import abiopen, abifile_subclass_from_filename
            try:
                _ = abifile_subclass_from_filename(obj)
                use_abiopen = True
            except ValueError:
                # This is needed to treat the case in which we are trying to read ElectronBands
                # from a nc file that is not known to AbiPy.
                use_abiopen = False

            if use_abiopen:
                with abiopen(obj) as abifile:
                    return abifile.ebands
            else:
                return cls.from_file(obj)

        elif hasattr(obj, "ebands"):
            # object with ebands
            return obj.ebands

        raise TypeError("Don't know how to extract ebands from object `%s`" % type(obj))

    @classmethod
    def from_mpid(cls, material_id, api_key=None, endpoint=None,
                  nelect=None, has_timerev=True,
                  nspinor=1, nspden=None, line_mode=True) -> ElectronBands:
        """
        Read band structure data corresponding to a materials project ``material_id``.
        and return Abipy ElectronBands object. Return None if bands are not available.

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
                If None, this value is automatically computed using the Fermi level (if metal)
                or the VBM indices reported in the JSON document sent by the MP database.
                It is recommended to pass nelect explictly if this quantity is known in the caller.
            nspinor: Number of spinor components.
            line_mode (bool): If True, fetch a BandStructureSymmLine object
                (default). If False, return the uniform band structure.
        """
        material_id = str(material_id)
        if not material_id.startswith("mp-"):
            raise ValueError("Materials project ID should start with mp-")

        # Get pytmatgen structure and convert it to an AbiPy structure
        from abipy.core import restapi
        with restapi.get_mprester(api_key=api_key, endpoint=endpoint) as rest:
            pmgb = rest.get_bandstructure_by_material_id(material_id=material_id, line_mode=line_mode)
            if pmgb is None: return None

            if pmgb.structure is None:
                # Structure is set to None so we have to perform another request and patch the object.
                structure = rest.get_structure_by_material_id(material_id, final=True)
                pmgb.structure = structure

        #ksampling = KSamplingInfo.from_kbounds(kbounds)
        return cls.from_pymatgen(pmgb, nelect, weights=None, has_timerev=has_timerev,
                                 ksampling=None, smearing=None, nspinor=nspinor, nspden=nspden)

    def to_json(self) -> str:
        """
        Returns a JSON string representation of the MSONable object.
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

        # Now I try to understand if we are dealing with a semiconductor or a metal
        # by using meta variables and the input fermie.
        # The algorithm may fail if we are dealing with a semimetal but in the internal implementation
        # we only discern between systems with/without an energy gap.
        # Note that in the case of NSCF calculations the fermie is taken from the previous SCF run
        # hence we assume the user used the same value of occopt both in the SCF and in the NSCF run
        # although the NSCF run does not need occopt.
        # This is automatically enforced if you use AbiPy flows but it's not guaranteed if the user
        # uses his/her own input files.

        if self.smearing and self.smearing.occopt == 1 and self.nsppol == 1 and self.nspden == 1:
            # This is the simplest case as the calculation has been done assuming a non-magnetic semiconductor
            # so there must be a gap.
            # On the other hand, in a NSCF run with a k-path it may happen that the HOMO energy
            # taken from the GS SCF run underestimates the real HOMO level.
            # For instance, the GS-SCF may have been done with a shifted k-mesh whereas
            # the true HOMO derived from the k-path is at Gamma.
            # Calling set_fermie_to_vbm should fix this.

            self.is_metal = False
            self.set_fermie_to_vbm()

        else:
            #print("This is a wannabe metal")
            # This flag is true if the system must be metallic in a single particle theory
            is_bloch_metal = self.nsppol == 1 and self.nspinor == 1 and self.nelect % 2 != 0

            self.is_metal = True
            if not is_bloch_metal:
                # Use the input fermie to understand if there's an energy gap.
                max_occ = 2.0 / (self.nsppol * self.nspinor)
                trial_occs = np.where(self.eigens > self.fermie, 0, max_occ)
                self.is_metal = False

                for ik in range(self.nkpt):
                    nele_k = trial_occs[:, ik].sum()
                    #print(nele_k)
                    if nele_k != self.nelect:
                        self.is_metal = True
                        break

    @property
    def structure(self) -> Structure:
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
            #if name is not None:
            if name:
                _auto_klabels[idx] = name
                if kpoint.name is None: kpoint.set_name(name)

        last = len(self.kpoints) - 1
        if last not in _auto_klabels: _auto_klabels[last] = " "

        return _auto_klabels

    def __repr__(self):
        """String representation (short version)"""
        #return "<%s, nk=%d, %s, id=%s>" % (self.__class__.__name__, self.nkpt, self.structure.formula, id(self))
        return "<%s, nk=%d, %s>" % (self.__class__.__name__, self.nkpt, self.structure.formula)

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

    def get_plotter_with(self, self_key, other_key, other_ebands):
        """
        Build and |ElectronBandsPlotter| from self and other, use self_key and other_key as keywords
        """
        plotter = ElectronBandsPlotter()
        plotter.add_ebands(self_key, self)
        plotter.add_ebands(other_key, other_ebands)

        return plotter

    #def get_panel(self, **kwargs):
    #    from abipy.panels.core import PanelWithElectronBands
    #    return PanelWithElectronBands(ebands=self).get_panel(**kwargs)

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
            self._nband = self.nband_sk[0, 0]
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
        if self.smearing:
            return self.smearing.has_metallic_scheme
        else:
            cprint("ebands.smearing is not defined, assuming has_metallic_scheme = False", "red")
            return False

    def set_fermie_to_vbm(self) -> float:
        """
        Set the Fermi energy to the valence band maximum (VBM).
        Useful when the initial fermie energy comes from a GS-SCF calculation
        that may underestimate the Fermi energy because e.g. the IBZ sampling
        is shifted whereas the true VMB is at Gamma.

        Return: New Fermi energy in eV.

        .. warning:

            Assume spin-unpolarized band energies.
        """
        if self.nsppol == 2:
            raise ValueError(f"set_fermie_to_vbm assumes nsppol == 1 while it is: {self.nsppol}")
        if self.nspden != 1:
            raise ValueError(f"set_fermie_to_vbm assumes nspden == 1 while it is: {self.nspden}")

        iv = int(self.nelect * self.nspinor) // 2 - 1
        if iv >= self.mband:
            iv = self.mband - 1
            cprint("mband = {self.mband} is not large enough to include the VBM", "red")

        new_fermie = self.eigens[:, :, iv].max()
        return self.set_fermie(new_fermie)

    #def set_fermie_to_midgap(self) -> float:
    #    new_fermie =
    #    return self.set_fermie(new_fermie)

    def set_fermie_from_edos(self, edos, nelect=None) -> float:
        """
        Set the Fermi level using the integrated DOS computed in edos.

         Args:
            edos: |ElectronDos| object.
            nelect: Number of electrons. If None, the number of electrons in self. is used

        Return: New fermi energy in eV.
        """
        if nelect is None:
            new_fermie = edos.find_mu(self.nelect)
        else:
            new_fermie = edos.find_mu(nelect)

        return self.set_fermie(new_fermie)

    def set_fermie(self, new_fermie):
        """Set the new fermi energy. Return new value"""
        self.fermie = new_fermie
        # TODO change occfacts
        return self.fermie

    def with_points_along_path(self, frac_bounds=None, knames=None, dist_tol=1e-12):
        """
        Build new |ElectronBands| object containing the k-points along the
        input k-path specified by `frac_bounds`. Useful to extract energies along a path
        from calculation performed in the IBZ.

        Args:
            frac_bounds: [M, 3] array  with the vertexes of the k-path in reduced coordinates.
                If None, the k-path is automatically selected from the structure.
            knames: List of strings with the k-point labels defining the k-path. It has precedence over frac_bounds.
            dist_tol: A point is considered to be on the path if its distance from the line
                is less than dist_tol.

        Return:
            namedtuple with the following attributes::

                ebands: |ElectronBands| object.
                ik_new2prev: Correspondence between the k-points in the new ebands and the kpoint
                    of the previous band structure (self).
        """
        # Construct the stars of the k-points for all k-points in self.
        # In principle, the input k-path is arbitrary and not necessarily in the IBZ used for self
        # so we have to build the k-stars and find the k-points lying along the path and keep
        # track of the mapping kpt --> star --> kgw
        # TODO: This part becomes a bottleneck for large nk!
        stars = [kpoint.compute_star(self.structure.abi_spacegroup.fm_symmops) for kpoint in self.kpoints]
        cart_coords, back2istar = [], []
        for istar, star in enumerate(stars):
            cart_coords.extend([k.cart_coords for k in star])
            back2istar.extend([istar] * len(star))
        cart_coords = np.reshape(cart_coords, (-1, 3))

        if knames is not None:
            assert frac_bounds is None
            frac_bounds = self.structure.get_kcoords_from_names(knames)
        else:
            if frac_bounds is None:
                frac_bounds = self.structure.calc_kptbounds()

        # Find (star) k-points on the path.
        cart_bounds = self.structure.reciprocal_lattice.get_cartesian_coords(frac_bounds)
        from abipy.core.kpoints import find_points_along_path
        p = find_points_along_path(cart_bounds, cart_coords, dist_tol=dist_tol)
        if len(p.ikfound) == 0:
            raise ValueError("Find zero points lying on the input k-path. Try to increase dist_tol")

        new_eigens = np.zeros((self.nsppol, len(p.ikfound), self.mband))
        new_occfacts = np.zeros_like(new_eigens)
        new_linewidths = None if self.linewidths is None else np.zeros_like(new_eigens)
        new_frac_coords = []

        # Correspondence new.kpoints --> self.ebands.kpoints
        # Useful if client code has to rearrange other arrays ordered according to self.ebands.kpoints.
        ik_new2prev = []
        for ik, ik_new in enumerate(p.ikfound):
            # Stars are ordered as self.kpoints to this is the index we need to access self.eigens.
            # and trasfer the data from self to new
            ik_self = back2istar[ik_new]
            ik_new2prev.append(ik_self)
            fcs = self.structure.reciprocal_lattice.get_fractional_coords(cart_coords[ik_new])
            #print("fcs", fcs, "dist", p.dist_list[ik])
            new_frac_coords.append(fcs)
            for spin in range(self.nsppol):
                new_eigens[spin, ik] = self.eigens[spin, ik_self]
                new_occfacts[spin, ik] = self.occfacts[spin, ik_self]
                if self.linewidths is not None:
                    new_linewidths[spin, ik] = self.linewidths[spin, ik_self]

        new_kpoints = Kpath(self.structure.reciprocal_lattice, new_frac_coords, weights=None, names=None)

        new_ebands = self.__class__(self.structure, new_kpoints, new_eigens, self.fermie, new_occfacts,
                             self.nelect, self.nspinor, self.nspden,
                             smearing=self.smearing, linewidths=new_linewidths)

        return dict2namedtuple(ebands=new_ebands, ik_new2prev=ik_new2prev)

    #def select_bands(self, bands, kinds=None):
    #    """Build new ElectronBands object by selecting bands via band_slice (slice object)."""
    #    bands = np.array(bands)
    #    kinds = np.array(kinds) if kinds is not None else np.array(range(self.nkpt))
    #    # This won't work because I need a KpointList object.
    #    new_kpoints = self.kpoints[kinds]
    #    new_eigens = self.eigens[:, kinds, bands].copy()
    #    new_occfacts = self.occupation[:, kinds, bands].copy()
    #    new_linewidths = None if not self.linewidths else self.linewidths[:, kinds, bands].copy()

    #    return self.__class__(self.structure, new_kpoints, new_eigens, self.fermie, new_occfacts,
    #                          self.nelect, self.nspinor, self.nspden,
    #                          smearing=self.smearing, linewidths=new_linewidths)

    @classmethod
    def empty_with_ibz(cls, ngkpt, structure, fermie, nelect, nsppol, nspinor, nspden, mband,
                       shiftk=(0, 0, 0), kptopt=1,
                       smearing=None, linewidths=None) -> ElectronBands:

        from abipy.abio.factories import gs_input
        from abipy.data.hgh_pseudos import HGH_TABLE
        gsinp = gs_input(structure, HGH_TABLE, spin_mode="unpolarized")
        ibz = gsinp.abiget_ibz(ngkpt=ngkpt, shiftk=shiftk, kptopt=kptopt)
        ksampling = KSamplingInfo.from_mpdivs(ngkpt, shiftk, kptopt)

        kpoints = IrredZone(structure.reciprocal_lattice, ibz.points, weights=ibz.weights,
                            names=None, ksampling=ksampling)

        new_eigens = np.zeros((nsppol, len(kpoints), mband))
        new_occfacts = np.zeros_like(new_eigens)

        return cls(structure, kpoints, new_eigens, fermie, new_occfacts,
                   nelect, nspinor, nspden,
                   smearing=smearing, linewidths=linewidths)

    def get_dict4pandas(self, with_geo=True, with_spglib=True) -> dict:
        """
        Return a :class:`OrderedDict` with the most important parameters:

            - Chemical formula and number of atoms.
            - Lattice lengths, angles and volume.
            - The spacegroup number computed by Abinit (set to None if not available).
            - The spacegroup number and symbol computed by spglib (set to None not `with_spglib`).

        Useful to construct pandas DataFrames

        Args:
            with_geo: True if structure info should be added to the dataframe
            with_spglib: If True, spglib_ is invoked to get the spacegroup symbol and number.
        """
        odict = OrderedDict([
            ("nsppol", self.nsppol), ("nspinor", self.nspinor), ("nspden", self.nspden),
            ("nkpt", self.nkpt), ("nband", self.nband_sk.min()),
            ("nelect", self.nelect), ("fermie", self.fermie),

        ])

        # Add info on structure.
        if with_geo:
            odict.update(self.structure.get_dict4pandas(with_spglib=with_spglib))

        odict.update(self.smearing)

        bws = self.bandwidths
        for spin in self.spins:
            odict["bandwidth_spin%d" % spin] = bws[spin]

        enough_bands = (self.mband > self.nspinor * self.nelect // 2)
        if enough_bands:
            for spin in self.spins:
                odict["fundgap_spin%d" % spin] = self.fundamental_gaps[spin].energy
            for spin in self.spins:
                odict["dirgap_spin%d" % spin] = self.direct_gaps[spin].energy

            # Select min gap over spins.
            min_fgap = self.fundamental_gaps[0]
            min_dgap = self.direct_gaps[0]
            if self.nsppol == 2:
                fgap0, fgap1 = self.fundamental_gaps[0], self.fundamental_gaps[1]
                min_fgap = fgap0 if fgap0.energy < fgap1.energy else fgap1
                dgap0, dgap1 = self.direct_gaps[0], self.direct_gaps[1]
                min_dgap = dgap0 if dgap0.energy < dgap1.energy else dgap1

            # These quantities are not spin-dependent.
            odict["gap_type"] = "direct" if min_fgap.is_direct else "indirect"
            odict["fundgap_kstart"] = repr(min_fgap.in_state.kpoint)
            odict["fundgap_kend"] = repr(min_fgap.out_state.kpoint)
            odict["dirgap_kstart"] = repr(min_dgap.in_state.kpoint)
            odict["dirgap_kend"] = repr(min_dgap.out_state.kpoint)

        return odict

    @lazy_property
    def has_bzmesh(self) -> bool:
        """True if the k-point sampling is homogeneous."""
        return isinstance(self.kpoints, IrredZone)

    @lazy_property
    def has_bzpath(self) -> bool:
        """True if the bands are computed on a k-path."""
        return isinstance(self.kpoints, Kpath)

    @lazy_property
    def kptopt(self) -> int:
        """The value of the kptopt input variable."""
        try:
            return self.kpoints.ksampling.kptopt
        except AttributeError:
            cprint("ebands.kpoints.ksampling.kptopt is not defined, assuming kptopt = 1", "red")
            return 1

    @lazy_property
    def has_timrev(self) -> bool:
        """True if time-reversal symmetry is used in the BZ sampling."""
        return has_timrev_from_kptopt(self.kptopt)

    def isnot_ibz_sampling(self, require_gamma_centered=False):
        """
        Test whether the k-points in the band structure represent an IBZ with an associated k-mesh
        Return string with error message if the condition is not fullfilled.

        Args:
            require_gamma_centered: True if the k-mesh should be Gamma-centered
        """
        errors = []
        eapp = errors.append

        if not self.kpoints.is_ibz:
            eapp("Expecting an IBZ sampling but got type `%s`" % type(self.kpoints))

        if not self.kpoints.is_mpmesh:
            eapp("Note that homogeneous k-meshes are required for the FS.\nksampling:` %s`" % str(self.kpoints.ksampling))

        if self.kpoints.is_mpmesh and require_gamma_centered:
            mpdivs, shifts = self.kpoints.mpdivs_shifts
            if shifts is not None and not np.all(shifts == 0.0):
                eapp(f"k-mesh should be gamma-centered but shifts: {shifts}")

        return "\n".join(errors)

    def kindex(self, kpoint):
        """
        The index of the k-point in the internal list of k-points.
        Accepts: |Kpoint| instance or integer.
        """
        if duck.is_intlike(kpoint):
            return int(kpoint)
        else:
            return self.kpoints.index(kpoint)

    def skb_iter(self):
        """Iterator over (spin, k, band) indices."""
        for spin in self.spins:
            for ik in self.kidxs:
                for band in range(self.nband_sk[spin, ik]):
                    yield spin, ik, band

    def deepcopy(self):
        """Deep copy of the ElectronBands object."""
        new = copy.deepcopy(self)
        # Don't know why but we need to copy the smearing manually.
        new.smearing = Smearing(**self.smearing)
        return new

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

    #@memoized_method(maxsize=5, typed=False)
    def get_dataframe(self, e0="fermie", brange=None, ene_range=None) -> pd.DataFrame:
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
                The Fermi energy is saved in df.fermie
            brange: If not None, only bands such as ``brange[0] <= band_index < brange[1]`` are included.
            ene_range: If not None, only bands whose energy in inside [erange[0],  [erange[1]] are included.
                brange and ene_range are mutually exclusive. Note that e0 is taken into account.
                when computing the energy.
        """
        if brange and ene_range:
            raise ValueError("brange and ene_range are mutually exclusive")

        e0 = self.get_e0(e0)

        if ene_range:
            # Compute brange from ene_range
            min_band, max_band = self.mband, -1
            for spin in self.spins:
                for k, kpoint in enumerate(self.kpoints):
                    e_ks = self.eigens[spin, k] - e0
                    bands = np.where((e_ks >= ene_range[0]) & (e_ks <= ene_range[1]))[0]
                    if bands.size:
                        min_band = min(min_band, int(bands[0]))
                        max_band = max(max_band, int(bands[-1]))
                        #raise ValueError(f"No band found in ene_range: {ene_range} with e0: {e0}")

            if max_band == min_band: max_band += 1
            if min_band > max_band:
                # Wrong interval provided by user --> Show everything
                min_band, max_band = 0, self.mband

            brange = (min_band, max_band)
            #print(min_band, max_band)

        rows = []
        for spin in self.spins:
            for k, kpoint in enumerate(self.kpoints):
                bands = range(self.nband_sk[spin, k]) if brange is None else \
                        range(brange[0], brange[1])

                for band in bands:
                    eig = self.eigens[spin, k, band] - e0
                    rows.append(OrderedDict([
                               ("spin", spin),
                               ("kidx", k),
                               ("band", band),
                               ("eig", eig),
                               ("occ", self.occfacts[spin, k, band]),
                               ("kpoint", self.kpoints[k]),
                            ]))

        df = pd.DataFrame(rows, columns=list(rows[0].keys()))
        df.fermie = e0

        return df

    @add_fig_kwargs
    def boxplot(self, ax=None, e0="fermie", brange=None, ene_range=None, swarm=False, **kwargs):
        """
        Use seaborn to draw a box plot to show the distributions of the eigenvalues
        with respect to the band index.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV.
                -  None: Don't shift energies, equivalent to ``e0 = 0``.
            brange: Only bands such as ``brange[0] <= band_index < brange[1]`` are included in the plot.
            ene_range: If not None, only bands whose energy in inside [erange[0],  [erange[1]] are included.
                brange and ene_range are mutually exclusive. Note that e0 is taken into account.
                when computing the energy window.
            swarm: True to show the datapoints on top of the boxes.
            kwargs: Keyword arguments passed to seaborn boxplot.

        Return: |matplotlib-Figure|
        """
        df = self.get_dataframe(e0=e0, brange=brange, ene_range=ene_range)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        import seaborn as sns
        hue = None if self.nsppol == 1 else "spin"
        ax = sns.boxplot(x="band", y="eig", data=df, hue=hue, ax=ax, **kwargs)
        if swarm:
            sns.swarmplot(x="band", y="eig", data=df, hue=hue, color=".25", ax=ax)

        return fig

    @add_plotly_fig_kwargs
    def boxplotly(self, e0="fermie", brange=None, ene_range=None, swarm=False, fig=None, rcd=None, **kwargs):
        """
        Use ployly to draw a box plot to show the distributions of the eigenvalues
        with respect to the band index.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV.
                -  None: Don't shift energies, equivalent to ``e0 = 0``.
            brange: Only bands such as ``brange[0] <= band_index < brange[1]`` are included in the plot.
            ene_range: If not None, only bands whose energy in inside [erange[0],  [erange[1]] are included.
                brange and ene_range are mutually exclusive. Note that e0 is taken into account.
                when computing the energy window.
            swarm: True to show the datapoints on top of the boxes.
            fig: plotly figure or None if a new figure should be created.
            rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col)
                of the subplot in the grid.
            kwargs: Keyword arguments passed to plotly px.box.

        Returns: |plotly.graph_objects.Figure|
        """
        df = self.get_dataframe(e0=e0, brange=brange, ene_range=ene_range)

        import plotly.express as px
        hue = None if self.nsppol == 1 else "spin"
        points = 'outliers' if not swarm else "all"
        px_fig = px.box(df, x="band", y="eig", color=hue, points=points, **kwargs)

        if rcd is None: return px_fig

        # Add px_fig traces to input fig with subplot.
        rcd = PlotlyRowColDesc.from_object(rcd)
        ply_row, ply_col, iax = rcd.ply_row, rcd.ply_col, rcd.iax
        for trace in px_fig.data:
            fig.add_trace(trace, row=ply_row, col=ply_col)

        return fig

    @classmethod
    def from_pymatgen(cls, pmg_bands, nelect, weights=None, has_timerev=True,
                      ksampling=None, smearing=None, nspinor=1, nspden=None) -> ElectronBands:
        """
        Convert a pymatgen band structure object to an Abipy |ElectronBands| object.

        Args:
            pmg_bands: pymatgen band structure object.
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

            The Abipy band structure contains more information than the pymatgen object so
            the conversion is not complete, especially if you rely on the default values.
            Please read carefylly the docstring and the code and use the optional arguments to pass
            additional data required by AbiPy if you need a complete conversion.
        """
        from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine

        # Cast to abipy structure and call spglib to init AbinitSpaceGroup.
        abipy_structure = Structure.as_structure(pmg_bands.structure.copy())
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

        fermie = pmg_bands.efermi
        #print("fermie energy from pymatgen bands:")

        # Compute occupation factors. Note that pmg bands don't have occfact so
        # I have to compute them from the eigens assuming T=0)
        atol = 1e-4
        abipy_occfacts = np.where(abipy_eigens <= fermie + atol, 1, 0)
        if nsppol == 1: abipy_occfacts *= 2

        if nelect is None:

            if pmg_bands.is_metal():
                nelect = np.rint(abipy_occfacts.sum() / nkpt)
                cprint(f"Using approximated method to get nelect in metals: {nelect}", color="yellow")
                #raise ValueError("Nelect must be specified if metallic bands.")

            else:
                #
                # Get nelect from valence band maximum index.
                #
                # - "band_index": A dict with spin keys pointing to a list of the
                # indices of the band containing the VBM (please note that you
                # can have several bands sharing the VBM) {Spin.up:[],
                # Spin.down:[]}

                d = pmg_bands.get_vbm()

                iv_up = max(d["band_index"][PmgSpin.up])
                homo_up = abipy_eigens[0, :, iv_up].max()
                homo = homo_up

                nelect = (iv_up + 1) * 2
                #print("iv_up", iv_up, "nelect: ", nelect)

                if pmg_bands.is_spin_polarized:
                    vbands_down = d["band_index"][PmgSpin.down]
                    iv_down = None
                    if vbands_down:
                        iv_down = max(vbands_down)
                        homo_down = abipy_eigens[1, :, iv_down].max()
                        homo = max(homo_up, homo_down)
                    nelect = np.count_nonzero(abipy_eigens[:,0,:] <= homo)

                    #cprint("Using approximated method to get nelect in metals: {nelect}", color="yellow")

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

        return cls(abipy_structure, abipy_kpoints, abipy_eigens, fermie, abipy_occfacts,
                   nelect, nspinor, nspden, smearing=smearing)

    def to_pymatgen(self):
        """
        Return a pymatgen band structure object from an Abipy |ElectronBands| object.
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
                                         labels_dict, coords_are_cartesian=False,
                                         structure=self.structure, projections=None)
        else:
            return BandStructure(self.kpoints.frac_coords, eigenvals, self.reciprocal_lattice, self.fermie,
                                 labels_dict=None, coords_are_cartesian=False,
                                 structure=self.structure, projections=None)

    def _electron_state(self, spin, kpoint, band):
        """
        Build an instance of :class:`Electron` from the spin, kpoint and band index
        """
        kidx = self.kindex(kpoint)
        return Electron(spin=spin,
                        kpoint=self.kpoints[kidx],
                        band=band,
                        eig=self.eigens[spin, kidx, band],
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
        """homo states for each spin channel as a list of nsppol :class:`Electron`."""
        homos = self.nsppol * [None]

        for spin in self.spins:
            blist, enes = [], []
            for k in self.kidxs:
                # Find rightmost value less than or equal to fermie.
                # Well, it's possible to have all eigens > fermie for particular k-points e.g. Al.
                try:
                    b = find_le(self.eigens[spin,k,:], self.fermie + self.pad_fermie)
                except ValueError:
                    #print("fermie + pad:", self.fermie + self.pad_fermie)
                    #print("eigens[spin,k,:]", self.eigens[spin,k,:])
                    continue

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

            #enes = np.array(enes)
            #kinds = np.where(enes == enes.min())[0]
            #lumo_kidx = np.asscalar(kinds[len(kinds) // 2])
            lumo_kidx = np.array(enes).argmin()
            lumo_band = blist[lumo_kidx]

            # Build Electron instance.
            lumos[spin] = self._electron_state(spin, lumo_kidx, lumo_band)

        return lumos

    #def get_nkpts_in_ewin(self, e0: float, ewin_ev: float = 1e-4):
    #   """Return the number of k-points inside an energy window centered at e0 of width `ewin_ev`."""
    #    count = 0
    #    for spin in self.spins:
    #        for ik in self.kidxs:
    #            if any(abs(self.eigens[spin, ik] - e0) < tol_ev): count += 1
    #    return count


    def get_edge_state(self, vbm_or_cbm, spin=None):

        if spin is None:
            # Returm max/min over spins (if any)

            if vbm_or_cbm == "vbm":
                ord_states = sorted(self.homos, key=lambda state: state.eig)
                return ord_states[-1]
            if vbm_or_cbm == "cbm":
                ord_states = sorted(self.lumos, key=lambda state: state.eig)
                return ord_states[0]

            raise ValueError(f"Invalid value for vbm_or_cbm: {vbm_or_cbm}")

        else:
            if vbm_or_cbm == "vbm": return self.homos[spin]
            if vbm_or_cbm == "cbm": return self.lumos[spin]
            raise ValueError(f"Invalid value for vbm_or_cbm: {vbm_or_cbm}")

    @property
    def bandwidths(self):
        """The bandwidth for each spin channel i.e. the energy difference (homo - lomo)."""
        return [self.homos[spin].eig - self.lomos[spin].eig for spin in self.spins]

    @property
    def fundamental_gaps(self) -> List[ElectronTransition]:
        """List of :class:`ElectronTransition` with info on the fundamental gaps for each spin."""
        return [ElectronTransition(self.homos[spin], self.lumos[spin]) for spin in self.spins]

    @property
    def direct_gaps(self) -> List[ElectronTransition]:
        """List of `nsppol` :class:`ElectronTransition` with info on the direct gaps for each spin."""
        dirgaps = self.nsppol * [None]
        for spin in self.spins:
            gaps = []
            for k in self.kidxs:
                homo_sk = self.homo_sk(spin, k)
                lumo_sk = self.lumo_sk(spin, k)
                gaps.append(lumo_sk.eig - homo_sk.eig)

            # Find the index of the k-point where the direct gap is located.
            # If there are multiple k-points along the path, prefer the one in the center
            # If not possible e.g. direct at G with G-X-L-G path avoid points on the right border of the graph
            gaps = np.array(gaps)
            kinds = np.where(gaps == gaps.min())[0]
            kdir = kinds[0]
            all_kinds = list(zip(kinds, kinds))
            #kdir = kinds[len(kinds) // 2]
            #kdir = np.array(gaps).argmin()
            dirgaps[spin] = ElectronTransition(self.homo_sk(spin, kdir), self.lumo_sk(spin, kdir), all_kinds=all_kinds)

        return dirgaps

    def get_gaps_string(self, with_latex=True, unicode=False) -> str:
        """
        Return string with info about fundamental and direct gap (if not metallic scheme)

        Args:
            with_latex: True to get latex symbols for the gap names and formula else text.
            unicode: True to get unicode symbols for the formula else text.
        """
        enough_bands = (self.mband > self.nspinor * self.nelect // 2)

        if with_latex:
            dg_name, fg_name = "$E^{dir}_{gap}$", "$E^{fund}_{gap}$"
            formula = self.structure.latex_formula
        else:
            dg_name, fg_name = "direct gap", "fundamental gap"
            formula = self.structure.formula

        if unicode:
            import re
            numl = re.findall(r'\d', formula)
            for s in numl:
                formula = formula.replace(s, SUBSCRIPT_UNICODE[s])

        #if enough_bands and not self.has_metallic_scheme:
        if enough_bands and not self.is_metal:
            if self.nsppol == 1:
                s = "%s: %s = %.2f, %s = %.2f (eV)" % (
                    formula,
                    dg_name, self.direct_gaps[0].energy,
                    fg_name, self.fundamental_gaps[0].energy)
            else:
                dgs = [t.energy for t in self.direct_gaps]
                fgs = [t.energy for t in self.fundamental_gaps]
                s = "%s: %s = %.2f spin ↑ (%.2f spin ↓), %s = %.2f spin ↑ (%.2f spin ↓) (eV)" % (
                    formula,
                    dg_name, dgs[0], dgs[1],
                    fg_name, fgs[0], fgs[1])
        else:
            s = ""

        return s

    def get_kpoints_and_band_range_for_edges(self):
        """
        Find the reduced coordinates and the band indice associate to the band edges.
        Important: Call set_fermie_to_vbm() to set the Fermi level to the VBM before calling this method.

        Return: (k0_list, effmass_bands_f90) (Fortran notation)
        """
        from collections import defaultdict
        k0_list, effmass_bands_f90 = [], []

        for spin in self.spins:
            d = defaultdict(lambda: [np.inf, -np.inf])
            homo, lumo = self.homos[spin], self.lumos[spin]
            k = tuple(homo.kpoint.frac_coords)
            d[k][0] = min(d[k][0], homo.band + 1) # C --> F index
            k = tuple(lumo.kpoint.frac_coords)
            d[k][1] = max(d[k][1], lumo.band + 1)

            for k in d:
                if d[k][0] == np.inf: d[k][0] = d[k][1]
                if d[k][1] == -np.inf: d[k][1] = d[k][0]
                if d[k][0] == np.inf or d[k][1] == -np.inf:
                    raise RuntimeError("Cannot find band extrema, dict:\n%s:" % str(d))

            for k, v in d.items():
                k0_list.append(k)
                effmass_bands_f90.append(v)

        # Set small values to zero.
        k0_list = np.reshape(k0_list, (-1, 3))
        k0_list = np.where(np.abs(k0_list) > 1e-12, k0_list, 0.0)

        effmass_bands_f90 = np.reshape(effmass_bands_f90, (-1, 2))
        #print("k0_list:\n", k0_list, "\neffmass_bands_f90:\n", effmass_bands_f90)

        return k0_list, effmass_bands_f90

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

        #if not self.has_metallic_scheme:
        if not self.is_metal:
            enough_bands = (self.mband > self.nspinor * self.nelect // 2)
            for spin in self.spins:
                if self.nsppol == 2: app(">>> For spin %s" % spin)
                if enough_bands:
                    # This can fail so we have to catch the exception.
                    try:
                        s = indent(self.direct_gaps[spin].to_string(verbose=verbose))
                        app(f"Direct gap:\n{s}")
                        s = indent(self.fundamental_gaps[spin].to_string(verbose=verbose))
                        app(f"Fundamental gap:\n{s}")
                    except Exception as exc:
                        app("WARNING: Cannot compute direct and fundamental gap.")
                        if verbose: app("Exception:\n%s" % str(exc))

                app("Bandwidth: %.3f (eV)" % self.bandwidths[spin])
                if verbose:
                    lomo = self.lomos[spin]
                    s = indent(lomo.to_string(verbose=verbose))
                    app(f"Valence minimum located at kpt index {lomo.kidx}:\n{s}")

                homo = self.homos[spin]
                s = indent(homo.to_string(verbose=verbose))
                app(f"Valence maximum located at kpt index {homo.kidx}:\n{s}")

                try:
                    # Cannot assume enough states for this!
                    lumo = self.lumos[spin]
                    s = indent(lumo.to_string(verbose=verbose))
                    app(f"Conduction minimum located at kpt index {lumo.kidx}:\n{s}\n")
                except Exception:
                    pass

            #app("TIP: Call set_fermie_to_vbm() to set the Fermi level to the VBM if this is a non-magnetic semiconductor\n")

            if not verbose:
                app("TIP: Use `--verbose` to print k-point coordinates with more digits")

        if with_kpoints:
            app(self.kpoints.to_string(verbose=verbose, title="K-points"))

        return "\n".join(lines)

    def new_with_irred_kpoints(self, prune_step=None) -> ElectronBands:
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

        # Extract eigenvalues and occupation factors associated to irred k-points.
        new_eigens = self.eigens[:, irred_map, :].copy()
        new_occfacts = self.occfacts[:, irred_map, :].copy()

        return self.__class__(self.structure, new_kpoints, new_eigens, self.fermie, new_occfacts,
                              self.nelect, self.nspinor, self.nspden, smearing=self.smearing)

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

    def get_edos(self, method="gaussian", step=0.1, width=0.2) -> ElectronDos:
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
            raise NotImplementedError(f"{method} method is not supported")

        # Use fermie from Abinit if we are not using metallic scheme for occopt.
        fermie = None
        #if self.smearing["occopt"] == 1:
        #    print("using fermie from GSR")
        #    fermie = self.fermie
        edos = ElectronDos(mesh, dos, self.nelect, fermie=fermie)
        #print("ebands.fermie", self.fermie, "edos.fermie", edos.fermie)
        return edos

    def compare_gauss_edos(self, widths, step=0.1) -> ElectronDosPlotter:
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
            label = r"$\sigma = %s$ (eV)" % width
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
            arrow_opts.update(dict(lw=2, arrowstyle="-|>",))
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

    def get_ejdos(self, spin, valence, conduction,
                  method="gaussian", step=0.1, width=0.2, mesh=None) -> Function1D:
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
                    jdos.plot_ax(ax, color=color, lw=lw, label=r"$v=%s \rightarrow c=%s, \sigma=%s$" % (v, c, s))
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
             points=None, with_gaps=False, max_phfreq=None, fontsize=8, **kwargs):
        r"""
        Plot the electronic band structure with matplotlib.

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
            points: Marker object with the position and the size of the marker.
                Used for plotting purpose e.g. QP energies, energy derivatives...
            with_gaps: True to add markers and arrows showing the fundamental and the direct gap.
                IMPORTANT: If the gaps are now showed correctly in a non-magnetic semiconductor,
                    call `ebands.set_fermie_to_vbm()` to align the Fermi level at the top of the valence
                    bands before executing `ebands.plot().
                    The Fermi energy stored in the object, indeed, comes from the GS calculation
                    that produced the DEN file. If the k-mesh used for the GS and the CBM is e.g. at Gamma,
                    the Fermi energy will be underestimated and a manual aligment is needed.
            max_phfreq: Max phonon frequency in eV to activate scatterplot showing
                possible phonon absorption/emission processes based on energy-conservation alone.
                All final states whose energy is within +- max_phfreq of the initial state are included.
                By default, the four electronic states defining the fundamental and the direct gaps
                are considered as initial state (not available for metals).
            fontsize: fontsize for legends and titles

        Returns: |matplotlib-Figure|
        """
        # Select spins
        spin_list = self.spins if spin is None else [spin]

        # Select the band range.
        if band_range is None:
            band_list = list(range(self.mband))
        else:
            band_list = list(range(band_range[0], band_range[1], 1))

        e0 = self.get_e0(e0)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, klabels=klabels)
        set_axlims(ax, ylims, "y")

        # Plot the band energies.
        for spin in spin_list:
            opts = {"color": "black", "linewidth": 2.0} if spin == 0 else \
                   {"color": "red", "linewidth": 2.0}
            # This to pass kwargs to plot_ax and avoid both lw and linewidth in opts
            if "lw" in kwargs: opts.pop("linewidth")
            opts.update(kwargs)

            for ib, band in enumerate(band_list):
                if ib != 0: opts.pop("label", None)
                self.plot_ax(ax, e0, spin=spin, band=band, **opts)

        if points is not None:
            ax.scatter(points.x, np.array(points.y) - e0, s=np.abs(points.s), marker="o", c="b")

        if with_gaps and (self.mband > self.nspinor * self.nelect // 2):
            # Show fundamental and direct gaps for each spin.
            from matplotlib.patches import FancyArrowPatch
            for spin in self.spins:
                f_gap = self.fundamental_gaps[spin]
                d_gap = self.direct_gaps[spin]
                # Need arrows only if fundamental and direct gaps for this spin are different.
                need_arrows = f_gap != d_gap

                arrow_opts = {"color": "k"} if spin == 0 else {"color": "red"}
                arrow_opts.update(lw=2, alpha=0.6, arrowstyle="->", connectionstyle='arc3',
                                  mutation_scale=20, zorder=1000)
                scatter_opts = {"color": "blue"} if spin == 0 else {"color": "green"}
                scatter_opts.update(marker="o", alpha=1.0, s=80, zorder=100, edgecolor='black')

                # Fundamental gap.
                mgap = -1
                for ik1, ik2 in f_gap.all_kinds:
                    posA = (ik1, f_gap.in_state.eig - e0)
                    posB = (ik2, f_gap.out_state.eig - e0)
                    mgap = max(mgap, posA[1], posB[1])
                    ax.scatter(posA[0], posA[1], **scatter_opts)
                    ax.scatter(posB[0], posB[1], **scatter_opts)
                    if need_arrows:
                        ax.add_patch(FancyArrowPatch(posA=posA, posB=posB, **arrow_opts))

                if d_gap != f_gap:
                    # Direct gap.
                    for ik1, ik2 in d_gap.all_kinds:
                        posA = (ik1, d_gap.in_state.eig - e0)
                        posB = (ik2, d_gap.out_state.eig - e0)
                        mgap = max(mgap, posA[1], posB[1])
                        ax.scatter(posA[0], posA[1], **scatter_opts)
                        ax.scatter(posB[0], posB[1], **scatter_opts)
                        if need_arrows:
                            ax.add_patch(FancyArrowPatch(posA=posA, posB=posB, **arrow_opts))

            # Try to set nice limits if not given by user.
            if ylims is None:
                set_axlims(ax, (-mgap - 5, +mgap + 5), "y")

            gaps_string = self.get_gaps_string()
            if gaps_string:
                ax.set_title(gaps_string, fontsize=fontsize)

        if max_phfreq is not None and (self.mband > self.nspinor * self.nelect // 2):
            # Add markers showing phonon absorption/emission processes.
            for spin in self.spins:
                #scatter_opts = {"color": "steelblue"} if spin == 0 else {"color": "teal"}
                scatter_opts = dict(alpha=0.4, s=40, zorder=10)
                items = (["fundamental_gaps", "direct_gaps"], ["in_state", "out_state"])
                items = list(enumerate(itertools.product(*items)))
                for i, (gap_name, state_name) in items:
                    # Use getattr to extract gaps, equivalent to:
                    #   gap = self.fundamental_gaps[spin]
                    #   e_start = gap.out_state.eig
                    gap = getattr(self, gap_name)[spin]
                    e_start = getattr(gap, state_name).eig
                    scatter_opts["marker"] = "o"
                    scatter_opts["color"] = plt.get_cmap("cool" if spin == 0 else "summer")(i/len(items))

                    for band in range(self.mband):
                        eks = self.eigens[spin, :, band]
                        where = np.where(np.abs(e_start - eks) <= max_phfreq)[0]
                        if not np.any(where): continue
                        ax.scatter(where, eks[where] - e0, **scatter_opts)

        return fig

    @add_plotly_fig_kwargs
    def plotly(self, spin=None, band_range=None, klabels=None, e0="fermie", fig=None, rcd=None, ylims=None,
               points=None, with_gaps=False, max_phfreq=None, fontsize=12, **kwargs):
        r"""
        Plot the electronic band structure with plotly.

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
            fig: plotly figure or None if a new figure should be created.
            rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col) of the subplot in the grid.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
            points: Marker object with the position and the size of the marker.
                Used for plotting purpose e.g. QP energies, energy derivatives.
            with_gaps: True to add markers and arrows showing the fundamental and the direct gap.
                IMPORTANT: If the gaps are now showed correctly in a non-magnetic semiconductor,
                call `ebands.set_fermie_to_vbm()` to align the Fermi level at the top of the valence
                bands before executing `ebands.plot().
                The Fermi energy stored in the object, indeed, comes from the GS calculation
                that produced the DEN file. If the k-mesh used for the GS and the CBM is e.g. at Gamma,
                the Fermi energy will be underestimated and a manual aligment is needed.
            max_phfreq: Max phonon frequency in eV to activate scatterplot showing
                possible phonon absorption/emission processes based on energy-conservation alone.
                All final states whose energy is within +- max_phfreq of the initial state are included.
                By default, the four electronic states defining the fundamental and the direct gaps
                are considered as initial state (not available for metals).
            fontsize: fontsize for legends and titles
            kwargs: Passed to fig.add_scatter method.

        Returns: |plotly.graph_objects.Figure|
        """
        # Select spins
        spin_list = self.spins if spin is None else [spin]

        # Select the band range.
        if band_range is None:
            band_list = list(range(self.mband))
        else:
            band_list = list(range(band_range[0], band_range[1], 1))

        e0 = self.get_e0(e0)
        fig, _ = get_fig_plotly(fig=fig)
        rcd = PlotlyRowColDesc.from_object(rcd)
        ply_row, ply_col, iax = rcd.ply_row, rcd.ply_col, rcd.iax

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_plotly(fig, klabels=klabels, iax=iax)
        plotly_set_lims(fig, ylims, "y")

        # Plot the band energies.
        for spin in spin_list:
            lw = kwargs.pop("lw", 2.0)
            line_opts = {"color": "black", "width": lw} if spin == 0 else {"color": "red", "width": lw}

            for ib, band in enumerate(band_list):
                if ib != 0: kwargs.pop("label", None)
                self.plotly_traces(fig, e0, rcd=rcd, spin=spin, band=band, line_opts=line_opts, **kwargs)

        if points is not None:
            fig.add_scatter(x=points.x, y=np.array(points.y) - e0, mode='markers', showlegend=False, row=ply_row,
                            col=ply_col, marker=dict(color='blue', size=np.abs(points.s), opacity=0.6, line_width=0))

        if with_gaps and (self.mband > self.nspinor * self.nelect // 2):
            # Show fundamental and direct gaps for each spin.
            from plotly.figure_factory import create_quiver
            for spin in self.spins:
                f_gap = self.fundamental_gaps[spin]
                d_gap = self.direct_gaps[spin]
                # Need arrows only if fundamental and direct gaps for this spin are different.
                need_arrows = f_gap != d_gap

                arrow_opts = {"color": "gray"} if spin == 0 else {"color": "orange"}
                scatter_opts = {"color": "blue"} if spin == 0 else {"color": "green"}
                scatter_opts.update(opacity=0.9, size=12, line_width=2)

                # Fundamental gap.
                mgap = -1
                for ik1, ik2 in f_gap.all_kinds:
                    posA = (ik1, f_gap.in_state.eig - e0)
                    posB = (ik2, f_gap.out_state.eig - e0)
                    mgap = max(mgap, posA[1], posB[1])
                    fig.add_scatter(x=[posA[0], posB[0]], y=[posA[1], posB[1]], mode='markers', name='',
                                    showlegend=False, marker=scatter_opts, row=ply_row, col=ply_col)
                    if need_arrows:
                        figcq = create_quiver(x=[posA[0]], y=[posA[1]], u=[posB[0]-posA[0]], v=[posB[1]-posA[1]],
                                              name='', scale=1, arrow_scale=0.2, showlegend=False, hoverinfo='none',
                                              marker=arrow_opts, line=dict(width=2))
                        fig.add_trace(figcq.data[-1], row=ply_row, col=ply_col)

                if d_gap != f_gap:
                    # Direct gap.
                    for ik1, ik2 in d_gap.all_kinds:
                        posA = (ik1, d_gap.in_state.eig - e0)
                        posB = (ik2, d_gap.out_state.eig - e0)
                        mgap = max(mgap, posA[1], posB[1])
                        fig.add_scatter(x=[posA[0],posB[0]], y=[posA[1],posB[1]], mode='markers', name='',
                                                 showlegend=False, marker=scatter_opts, row=ply_row, col=ply_col)
                        if need_arrows:
                            figcq = create_quiver(x=[posA[0]], y=[posA[1]], u=[posB[0]-posA[0]], v=[posB[1]-posA[1]],
                                                  name='', scale=1, arrow_scale=0.2, showlegend=False, hoverinfo='none',
                                                  marker=arrow_opts, line=dict(width=2))
                            fig.add_trace(figcq.data[-1], row=ply_row, col=ply_col)

            # Try to set nice limits if not given by user.
            if ylims is None:
                plotly_set_lims(fig, (-mgap - 5, +mgap + 5), "y")

            gaps_string = self.get_gaps_string(with_latex=False, unicode=True)
            if gaps_string:
                if fig.layout.annotations == ():
                    fig.layout.annotations = [dict(text=gaps_string, font_size=fontsize, x=0, xref='paper',
                                               xanchor='left', y=1, yref='paper', yanchor='bottom', showarrow=False)]
                else:
                    fig.layout.annotations[iax-1].text = gaps_string
                    fig.layout.annotations[iax-1].font.size = fontsize

        if max_phfreq is not None and (self.mband > self.nspinor * self.nelect // 2):
            # Add markers showing phonon absorption/emission processes.
            for spin in self.spins:
                #scatter_opts = {"color": "steelblue"} if spin == 0 else {"color": "teal"}
                scatter_opts = dict(opacity=0.4, size=8)
                items = (["fundamental_gaps", "direct_gaps"], ["in_state", "out_state"])
                items = list(enumerate(itertools.product(*items)))
                for i, (gap_name, state_name) in items:
                    # Use getattr to extract gaps, equivalent to:
                    #   gap = self.fundamental_gaps[spin]
                    #   e_start = gap.out_state.eig
                    gap = getattr(self, gap_name)[spin]
                    e_start = getattr(gap, state_name).eig
                    scatter_opts["color"] = i / len(items)
                    scatter_opts["colorscale"] = "dense" if spin == 0 else "Burgyl"

                    for band in range(self.mband):
                        eks = self.eigens[spin, :, band]
                        where = np.where(np.abs(e_start - eks) <= max_phfreq)[0]
                        if not np.any(where): continue
                        fig.add_scatter(x=where, y=eks[where] - e0, mode='markers',
                                        marker=scatter_opts, showlegend=False, row=ply_row, col=ply_col)

        return fig


    # TODO: Is this really useful?

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
                -  None: Don't shift energies, equivalent to ``e0 = 0``.
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

    def decorate_ax(self, ax, **kwargs) -> None:
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
        ax.set_xlabel("Wave Vector")

        # Set ticks and labels.
        klabels = kwargs.pop("klabels", None)
        ticks, labels = self._make_ticks_and_labels(klabels)
        if ticks:
            # Don't show label if previous k-point is the same.
            for il in range(1, len(labels)):
                if labels[il] == labels[il-1]: labels[il] = ""
            #print("ticks", ticks, "\nlabels", labels)
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False, size=kwargs.pop("klabel_size", "large"))
            #print("ticks", len(ticks), ticks)
            ax.set_xlim(ticks[0], ticks[-1])

    def decorate_plotly(self, fig, **kwargs) -> None:
        """
        Add q-labels and unit name to figure ``fig``.
        Use units="" to add k-labels without unit name.
        Args:
            klabels:
            klabel_size:
            iax: An int, use iax=n to decorate the nth axis when the fig has subplots.
        """
        iax = kwargs.pop("iax", 1)
        xaxis = 'xaxis%u' % iax

        fig.layout[xaxis].title.text = "Wave Vector"
        fig.layout['yaxis%u' % iax].title.text = "Energy (eV)"

        # Set ticks and labels.
        klabels = kwargs.pop("klabels", None)
        ticks, labels = self._make_ticks_and_labels(klabels)
        if ticks:
            labels = plotly_klabels(labels)
            fig.layout[xaxis].tickvals = ticks
            fig.layout[xaxis].ticktext = labels
            fig.layout[xaxis].tickfont.size = kwargs.pop("klabel_size", 16)
            fig.layout[xaxis].range = (ticks[0], ticks[-1])

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
        Helper function to plot the energies for (spin, band) on the axis ax with matplotlib..

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
        with_linewidths = kwargs.pop("with_linewidths", True) and self.has_linewidths
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

    def plotly_traces(self, fig, e0, rcd=None, spin=None, band=None, showlegend=False, line_opts=None, **kwargs):
        """
        Helper function to plot the energies for (spin, band) on figure ``fig`` with plotly.

        Args:
            fig: |plotly.graph_objects.Figure|.
            e0: Option used to define the zero of energy in the band structure plot.
            rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col) of the subplot in the grid.
            spin: Spin index. If None, all spins are plotted.
            band: Band index, If None, all bands are plotted.
            showlegend: Determines whether or not an item corresponding to this trace is shown in the legend.
            line_opts: Dict of linestyle options passed to |plotly.graph_objects.scatter.Line|
            kwargs: Passed to fig.add_scatter method.
        """
        spin_range = range(self.nsppol) if spin is None else [spin]
        band_range = range(self.mband) if band is None else [band]

        label = kwargs.pop("label", '')
        # Handle linewidths
        with_linewidths = kwargs.pop("with_linewidths", True) and self.has_linewidths
        if with_linewidths:
            lw_opts = kwargs.pop("lw_opts", dict(opacity=0.6))
            lw_fact = lw_opts.pop("fact", 2.0)

        marker_opts = kwargs.pop("marker", None)
        xx = np.arange(self.nkpt)
        e0 = self.get_e0(e0)
        for spin in spin_range:
            for band in band_range:
                yy = self.eigens[spin, :, band] - e0

                # Set label only at the first iteration
                rcd = PlotlyRowColDesc.from_object(rcd)
                ply_row, ply_col = rcd.ply_row, rcd.ply_col
                if marker_opts:
                    fig.add_scatter(x=xx, y=yy, mode="lines+markers", name=label, showlegend=showlegend, line=line_opts,
                                    marker=marker_opts, legendgroup=label, **kwargs, row=ply_row, col=ply_col)
                else:
                    fig.add_scatter(x=xx, y=yy, mode="lines", name=label, showlegend=showlegend, line=line_opts,
                                    legendgroup=label, **kwargs, row=ply_row, col=ply_col)
                showlegend = False

                if with_linewidths:
                    w = self.linewidths[spin, :, band] * lw_fact / 2
                    lw_opts.update({'color': "black" if spin == 0 else "red"})
                    fig.add_scatter(x=xx, y=yy - w, mode='lines', line=lw_opts, name='',
                                    showlegend=False, row=ply_row, col=ply_col)
                    fig.add_scatter(x=xx, y=yy + w, mode='lines', line=lw_opts, name='',
                                    showlegend=False, fill='tonexty', row=ply_row, col=ply_col)


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
    def plot_with_edos(self, edos, klabels=None, ax_list=None, e0="fermie", points=None,
                       with_gaps=False, max_phfreq=None, ylims=None, width_ratios=(2, 1), **kwargs):
        r"""
        Plot the band structure and the DOS with matplotlib.

        Args:
            edos: An instance of |ElectronDos|.
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. ``klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            ax_list: The axes for the band structure plot and the DOS plot. If ax_list is None, a new figure
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

            points: Marker object with the position and the size of the marker.
                Used for plotting purpose e.g. QP energies, energy derivatives...
            with_gaps: True to add markers and arrows showing the fundamental and the direct gap.
            max_phfreq: Max phonon frequency in eV to activate scatterplot showing
                possible phonon absorption/emission processes based on energy-conservation alone.
                All final states whose energy is within +- max_phfreq of the initial state are included.
                By default, the four electronic states defining the fundamental and the direct gaps
                are considered as initial state (not available for metals).
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
        self.plot(e0=e0, ax=ax0, ylims=ylims, klabels=klabels, points=points,
                  with_gaps=with_gaps, max_phfreq=max_phfreq, show=False)

        # Plot the DOS
        if self.nsppol == 1:
            opts = {"color": "black", "linewidth": 2.0}
            edos.plot_ax(ax1, e0, exchange_xy=True, **opts)
        else:
            for spin in self.spins:
                opts = {"color": "black", "linewidth": 2.0} if spin == 0 else \
                       {"color": "red", "linewidth": 2.0}
                edos.plot_ax(ax1, e0, spin=spin, exchange_xy=True, **opts)

        ax1.grid(True)
        ax1.yaxis.set_ticks_position("right")
        ax1.yaxis.set_label_position("right")
        set_axlims(ax1, ylims, "y")

        return fig

    @add_plotly_fig_kwargs
    def plotly_with_edos(self, edos, klabels=None, fig=None, band_rcd=None, dos_rcd=None, e0="fermie", points=None,
                         with_gaps=False, max_phfreq=None, ylims=None, width_ratios=(2, 1), **kwargs):
        r"""
        Plot the band structure and the DOS with plotly.

        Args:
            edos: An instance of |ElectronDos|.
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. ``klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}``.
            fig: The |plotly.graph_objects.Figure| with two distinct plots for the band structure plot and the DOS plot.
                If fig is None, a new figure is created.
            band_rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col)
                of the band subplot in the grid.
            dos_rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col) of the dos subplot in the grid.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
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

            points: Marker object with the position and the size of the marker.
                Used for plotting purpose e.g. QP energies, energy derivatives...
            with_gaps: True to add markers and arrows showing the fundamental and the direct gap.
            max_phfreq: Max phonon frequency in eV to activate scatterplot showing
                possible phonon absorption/emission processes based on energy-conservation alone.
                All final states whose energy is within +- max_phfreq of the initial state are included.
                By default, the four electronic states defining the fundamental and the direct gaps
                are considered as initial state (not available for metals).
            width_ratios: Defines the ratio between the band structure plot and the dos plot.

        Return: |plotly.graph_objects.Figure|
        """
        if fig is None:
            # Build fig.
            fig, _ = get_figs_plotly(nrows=1, ncols=2, sharex=False, sharey=True,
                                     horizontal_spacing=0.02, column_widths=width_ratios)
            band_rcd = PlotlyRowColDesc(0, 0, 1, 2)
            dos_rcd = PlotlyRowColDesc(0, 1, 1, 2)

        # Define the zero of energy.
        e0 = self.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
        #if not kwargs: kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the band structure
        self.plotly(e0=e0, fig=fig, rcd=band_rcd, ylims=ylims, klabels=klabels, points=points,
                  with_gaps=with_gaps, max_phfreq=max_phfreq, show=False)

        # Plot the DOS
        if self.nsppol == 1:
            opts = {"color": "black", "width": 2.0}
            edos.plotly_traces(fig, e0, exchange_xy=True, rcd=dos_rcd, line_opts=opts)
        else:
            for spin in self.spins:
                opts = {"color": "black", "width": 2.0} if spin == 0 else \
                       {"color": "red", "width": 2.0}
                edos.plotly_traces(fig, e0, spin=spin, exchange_xy=True, rcd=dos_rcd, line_opts=opts)

        plotly_set_lims(fig, ylims, "y")

        return fig

    @add_fig_kwargs
    def plot_lws_vs_e0(self, ax=None, e0="fermie", function=lambda x: x, exchange_xy=False,
                       xlims=None, ylims=None, fontsize=12, **kwargs):
        r"""
        Plot electronic linewidths vs KS energy with matplotlib.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            function: Apply this function to the values before plotting
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

        xx, yy = e0mesh, tuple([function(lw) for lw in lws])
        if exchange_xy: xx, yy = yy, xx
        ax.plot(xx, yy, kw_linestyle, color=kw_color, label=kw_label, **kwargs)
        #ax.scatter(xx, yy)

        ax.grid(True)
        ylabel = "Linewidth (eV)"
        if exchange_xy: xlabel, ylabel = ylabel, xlabel
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        if kw_linestyle:
            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def to_xmgrace(self, filepath: str) -> None:
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

    @memoized_method(maxsize=5, typed=False)
    def get_ifermi_dense_bs(self, interpolation_factor, with_velocities):
        """
        Use ifermi and BoltzTraP2 to interpolate KS energies (assumes ebands in the IBZ).

        Args:
            interpolation_factor:
            with_velocities: Interpolate velocities in the full BZ.

        .. note::

            Store results in per-instance cache via memoized_method.
        """
        err_msg = self.isnot_ibz_sampling()
        if err_msg:
            raise ValueError(err_msg)

        try:
            from ifermi.interpolate import FourierInterpolator
        except ImportError:
            raise ImportError("Cannot import ifermi package.\nPlease install the package " +
                              "following the instructions given at: https://github.com/fermisurfaces/IFermi")

        # interpolate the energies onto a dense k-point mesh
        bs = self.to_pymatgen()
        interpolator = FourierInterpolator(bs)

        nworkers = 1 # Use 1 worker because it does not seem to scale well on my Mac.
        with Timer(f"BoltzTraP2 interpolation with factor {interpolation_factor} and velocities: {with_velocities}"):
            if with_velocities:
                dense_bs, velocities = interpolator.interpolate_bands(interpolation_factor=interpolation_factor,
                                                                      return_velocities=with_velocities,
                                                                      nworkers=nworkers)
            else:
                dense_bs = interpolator.interpolate_bands(interpolation_factor=interpolation_factor,
                                                          return_velocities=with_velocities,
                                                          nworkers=nworkers)
                velocities = None

        return dict2namedtuple(dense_bs=dense_bs, velocities=velocities, interpolator=interpolator)

    def get_ifermi_fs(self, interpolation_factor=8, mu=0.0, eref="fermie", wigner_seitz=True,
                      calculate_dimensionality=False, with_velocities=False):
        """
        Use ifermi package to visualize the (interpolated) Fermi surface.
        Requires netcdf file with energies in the IBZ.
        See also <https://fermisurfaces.github.io/IFermi/>

        Args:
            interpolation_factor: The factor by which the band structure will be interpolated.
            mu: Energy offset from the reference energy determing by `eref`.
            eref: Defines the energy reference. Possible values: `fermie`, `cbm`, `vbm`.
                The energy of the isosurface is given by: `eref` + `mu`.
            wigner_seitz: Controls whether the cell is the Wigner-Seitz cell
                          or the reciprocal unit cell parallelepiped.
            calculate_dimensionality:
            with_velocities: Generate the Fermi surface and calculate the group velocity
                at the center of each triangular face.

        Returns:

        .. example::

            r = ebands.get_ifermi_fs()
            r.fs_plotter.get_plot(plot_type="plotly").show()

        """
        r = self.get_ifermi_dense_bs(interpolation_factor, with_velocities)

        from ifermi.surface import FermiSurface
        from ifermi.plot import FermiSurfacePlotter #, save_plot, show_plot FermiSlicePlotter,

        eref = eref.lower()
        if eref == "fermie":
            edge_state = None
            abs_isoenergy = self.fermie

        elif eref in ("vbm", "cbm"):

            if eref == "vbm" and mu >= 0.0:
                cprint("WARNING: when eref == 'vbm', mu is expected to be < 0", color="red")
            if eref == "cbm" and mu <= 0.0:
                cprint("WARNING: when eref == 'cbm', mu is expected to be > 0", color="red")

            edge_state = self.get_edge_state(eref)
            mu = -self.fermie + edge_state.eig + mu
            abs_isoenergy = edge_state.eig + mu

        else:
            raise ValueError(f"Invalid value for eref: {eref}")

        # generate the Fermi surface
        with Timer(f"Building Fermi surface with wigner_seitz: {wigner_seitz}, eref: {eref} and mu: {mu} (eV)"):

            from ifermi.kpoints import kpoints_from_bandstructure
            dense_kpoints = kpoints_from_bandstructure(r.dense_bs) if r.velocities else None

            fs = FermiSurface.from_band_structure(
                  r.dense_bs, mu=mu, wigner_seitz=wigner_seitz,
                  calculate_dimensionality=False,
                  property_data=r.velocities,
                  property_kpoints=dense_kpoints,
            )

        fs_plotter = FermiSurfacePlotter(fs)

        return dict2namedtuple(fs=fs, fs_plotter=fs_plotter, dense_bs=r.dense_bs, velocities=r.velocities,
                               interpolator=r.interpolator, edge_state=edge_state, abs_isoenergy=abs_isoenergy)

    #def get_ifermi_slices(self, interpolation_factor=5, mu=0, eref="cbm", wigner_seitz=True,
    #                      with_velocities=False):

    #    r = self.get_ifermi_fs(interpolation_factor=interpolation_factor, mu=mu, eref=eref,
    #                           wigner_seitz=wigner_seitz, with_velocities=with_velocities)

    #    #plane_normals = [(1, 0, 0), ]
    #    plane_normals = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    #    #plane_normals += [(1, 1, 0), (0, 1, 1), (1, 1, 1)]
    #    #plane_normals = [(1, 1, 1),]
    #    print("edge_state:", r.edge_state.kpoint)
    #    distance = r.edge_state.kpoint.norm  / units.bohr_to_ang / (2 * np.pi)
    #    #distance = r.edge_state.kpoint.norm / units.bohr_to_ang

    #    print("matrix1", r.edge_state.kpoint.lattice.matrix)
    #    print("matrix2", r.fs.structure.lattice.reciprocal_lattice.matrix)
    #    print("distance:", distance)
    #    #distance = 0
    #    print("reciprocal_lattice:\n", self.structure.reciprocal_lattice.as_dict(verbosity=1))
    #    from ifermi.plot import FermiSlicePlotter

    #    expose_web = True
    #    from abipy.tools.plotting import MplExpose, PanelExpose
    #    if expose_web:
    #        e = PanelExpose(title=f"e-Bands of {self.structure.formula}")
    #    else:
    #        e = MplExpose(verbose=1)

    #    with e:
    #        for plane_normal in plane_normals:
    #            fermi_slice = r.fs.get_fermi_slice(plane_normal=plane_normal, distance=distance)
    #            slice_plotter = FermiSlicePlotter(fermi_slice)
    #            plt = slice_plotter.get_plot()
    #            fig = plt.gcf()
    #            fig.suptitle(f"Plane normal to {plane_normal} at distance {distance:.3f} from the Γ-point")
    #            e(fig)

    def to_bxsf(self, filepath):
        """
        Export the full band structure to ``filepath`` in BXSF format
        suitable for the visualization of isosurfaces with xcrysden_ (xcrysden --bxsf FILE).
        Require k-points in IBZ and gamma-centered k-mesh.
        """
        err_msg = self.isnot_ibz_sampling(require_gamma_centered=True)
        if err_msg:
            raise ValueError(err_msg)

        self.get_ebands3d().to_bxsf(filepath)

    #@memoized_method(maxsize=5, typed=False)
    def get_ebands3d(self):
        err_msg = self.isnot_ibz_sampling()
        if err_msg:
            raise ValueError(err_msg)

        return ElectronBands3D(self.structure, self.kpoints, self.has_timrev, self.eigens, self.fermie)

    def derivatives(self, spin, band, order=1, acc=4):
        """
        Compute the derivative of the eigenvalues wrt to k using finite difference.

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

    #def effective_masses(self, spin, band, acc=4):
    #    """
    #    Compute the effective masses for the given ``spin`` and ``band`` index.
    #    Use finite difference with accuracy ``acc``.

    #    Returns:
    #        |numpy-array| of size self.nkpt with effective masses.
    #    """
    #    ders2 = self.derivatives(spin, band, order=2, acc=acc) * (units.eV_to_Ha / units.bohr_to_ang**2)
    #    return 1. / ders2

    def get_effmass_analyzer(self):
        """"
        Return an instance of EffMassAnalyzer to compute effective masses with finite differences
        """
        from abipy.electrons.effmass_analyzer import EffMassAnalyzer
        return EffMassAnalyzer(self, copy=True)

    def get_effmass_line(self, spin, kpoint, band, acc=4):
        """
        Compute the effective masses along a k-line. Requires band energies on a k-path.

        Args:
            spin: Spin index.
            kpoint: integer, list of fractional coordinates or |Kpoint| object.
            band: Band index.
            acc: accuracy
        """
        warnings.warn("You may want to use `emana = ebands.get_effmass_analyzer()` for a more flexible API.")
        if not self.kpoints.is_path:
            raise ValueError("get_effmass_line requires k-points along a path. Got:\n %s" % repr(self.kpoints))

        # We have to understand if the k-point is a vertex or not.
        # If it is a vertex, we have to compute the left and right derivative.
        # If kpt is inside the line, left and right derivatives are supposed to be equal.
        from abipy.tools.derivatives import finite_diff

        for ik in self.kpoints.get_all_kindices(kpoint):
            for iline, line in enumerate(self.kpoints.lines):
                if line[-1] >= ik >= line[0]: break
            else:
                raise ValueError("Cannot find k-index `%s` in lines: `%s`" % (ik, self.kpoints.lines))

            kpos = line.index(ik)
            is_inside = kpos not in (0, len(line) - 1)
            do_right = (not is_inside) and kpos != 0 and iline != len(self.kpoints.lines) - 1

            evals_on_line, h_left, vers_left = self._eigens_hvers_iline(spin, band, iline)
            d2 = finite_diff(evals_on_line, h_left, order=2, acc=acc, index=kpos)
            em_left = 1. / (d2.value * (units.eV_to_Ha / units.bohr_to_ang ** 2))
            em_right = em_left
            h_right, vers_right = h_left, vers_left

            if do_right:
                kpos_right = self.kpoints.lines[iline + 1].index(ik)
                assert kpos_right == 0
                evals_on_line, h_right, vers_right = self._eigens_hvers_iline(spin, band, iline + 1)
                d2 = finite_diff(evals_on_line, h_right, order=2, acc=acc, index=kpos_right)
                em_right = 1. / (d2.value * (units.eV_to_Ha / units.bohr_to_ang ** 2))

            lines = []; app = lines.append
            app("For spin: %s, band: %s, k-point: %s, eig: %.3f [eV], accuracy: %s" % (
                spin, band, repr(self.kpoints[ik]), self.eigens[spin, ik, band], acc))
            #app("K-point: %s, eigenvalue: %s (eV)" % (repr(self.kpoint), self.eig))
            #app("h_left: %s, h_right %s" % (self.h_left, self.h_right))
            #app("is_inside: %s, vers_left: %s, vers_right: %s" % (self.is_inside, self.vers_left, self.vers_right))
            if em_left != em_right:
                app("emass_left: %.3f, emass_right: %.3f" % (em_left, em_right))
            else:
                app("emass: %.3f" % em_left)

            print("\n".join(lines))

    def _eigens_hvers_iline(self, spin, band, iline):
        line = self.kpoints.lines[iline]
        evals_on_line = self.eigens[spin, line, band]
        h = self.kpoints.ds[line[0]]

        if not np.allclose(h, self.kpoints.ds[line[:-1]]):
            raise ValueError("For finite difference derivatives, the path must be homogeneous!\n" +
                             str(self.kpoints.ds[line[:-1]]))

        return evals_on_line, h, self.kpoints.versors[line[0]]

    def interpolate(self, lpratio=5, knames=None, vertices_names=None, line_density=20,
                    kmesh=None, is_shift=None, bstart=0, bstop=None, filter_params=None, verbose=0):
        """
        Interpolate energies in k-space along a k-path and, optionally, in the IBZ for DOS calculations.
        Note that the interpolation will likely fail if there are symmetrical k-points in the input set of k-points
        so it's recommended to call this method with energies obtained in the IBZ.

        Args:
            lpratio: Ratio between the number of star functions and the number of ab-initio k-points.
                The default should be OK in many systems, larger values may be required for accurate derivatives.
            knames: List of strings with the k-point labels for the k-path. Has precedence over ``vertices_names``.
            vertices_names: Used to specify the k-path for the interpolated band structure
                It's a list of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path. ``line_density`` defines
                the density of the sampling. If None, the k-path is automatically generated according
                to the point group of the system.
            line_density: Number of points in the smallest segment of the k-path.
                If 0, use list of k-points given in vertices_names
            kmesh: Used to activate the interpolation on the homogeneous mesh for DOS (uses spglib_ API).
                kmesh is given by three integers and specifies mesh numbers along reciprocal primitive axis.
            is_shift: three integers (spglib_ API). When is_shift is not None, the kmesh is shifted along
                the axis in half of adjacent mesh points irrespective of the mesh numbers. None means unshifted mesh.
            bstart, bstop: Select the range of band to be used in the interpolation
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

        if self.nband > self.nelect and self.nband > 20 and bstart == 0 and bstop is None:
            cprint("Bands object contains nband %s with nelect %s. You may want to use bstart, bstop to select bands." % (
                    self.nband, self.nelect), "yellow")

        # Build interpolator.
        from abipy.core.skw import SkwInterpolator
        cell = (self.structure.lattice.matrix, self.structure.frac_coords, self.structure.atomic_numbers)

        skw = SkwInterpolator(lpratio, self.kpoints.frac_coords, self.eigens[:,:,bstart:bstop], self.fermie, self.nelect,
                              cell, fm_symrel, self.has_timrev,
                              filter_params=filter_params, verbose=verbose)

        # Generate k-points for interpolation.
        if knames is not None:
            kpath = Kpath.from_names(self.structure, knames, line_density=line_density)
        else:
            if vertices_names is None:
                vertices_names = [(k.frac_coords, k.name) for k in self.structure.hsym_kpoints]
            kpath = Kpath.from_vertices_and_names(self.structure, vertices_names, line_density=line_density)

        # Interpolate energies.
        eigens_kpath = skw.interp_kpts(kpath.frac_coords).eigens

        # Build new ebands object.
        occfacts_kpath = np.zeros_like(eigens_kpath)
        ebands_kpath = self.__class__(self.structure, kpath, eigens_kpath, self.fermie, occfacts_kpath,
                                      self.nelect, self.nspinor, self.nspden, smearing=self.smearing)
        ebands_kmesh = None
        if kmesh is not None:
            # Get kpts and weights in the IBZ.
            kdos = Ktables(self.structure, kmesh, is_shift, self.has_timrev)
            eigens_kmesh = skw.interp_kpts(kdos.ibz).eigens

            # Build new ebands object with k-mesh
            #kptopt = kptopt_from_timrev()
            ksampling = KSamplingInfo.from_mpdivs(mpdivs=kmesh, shifts=[0, 0, 0], kptopt=1)
            kpts_kmesh = IrredZone(self.structure.reciprocal_lattice, kdos.ibz, weights=kdos.weights,
                                   names=None, ksampling=ksampling)
            occfacts_kmesh = np.zeros_like(eigens_kmesh)

            ebands_kmesh = self.__class__(self.structure, kpts_kmesh, eigens_kmesh, self.fermie, occfacts_kmesh,
                                          self.nelect, self.nspinor, self.nspden, smearing=self.smearing)

        return dict2namedtuple(ebands_kpath=ebands_kpath, ebands_kmesh=ebands_kmesh, interpolator=skw)

    def get_collinear_mag(self) -> float:
        """
        Calculates the total collinear magnetization in Bohr magneton as the difference
        between the spin up and spin down densities.
        Note that we assume an IBZ sampling.

        Returns:
            float: the total magnetization.

        Raises: ValueError if total magnetization is undefined e.g. nspinor 2.
        """
        if self.nsppol == 1:
            if self.nspinor == 1 or (self.nspinor == 2 and self.nspden == 1):
                return 0.0
            else:
                raise ValueError("Cannot calculate collinear magnetization for nsppol: {}, "
                                 "nspinor {}, nspden {}".format(self.nsppol, self.nspinor, self.nspden))
        else:
            rhoup = np.sum(self.kpoints.weights[:, None] * self.occfacts[0])
            rhoudown = np.sum(self.kpoints.weights[:, None] * self.occfacts[1])
            return rhoup - rhoudown


def dataframe_from_ebands(ebands_objects, index=None, with_spglib=True) -> pd.DataFrame:
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

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: ElectronBandsPlotter
    """
    # Used in iter_lineopt to generate matplotlib linestyles.
    _LINE_COLORS = ["blue", "red", "green", "magenta", "yellow", "black"]
    _LINE_STYLES = ["-", ":", "--", "-.",]
    _LINE_STYLES_PLOTLY = ['solid', "dot", 'dash', 'dashdot',]
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

    def __len__(self):
        return len(self.ebands_dict)

    def add_plotter(self, other) -> ElectronBandsPlotter:
        """Merge two plotters, return new plotter."""
        if not isinstance(other, self.__class__):
            raise TypeError("Don't know to add %s to %s" % (other.__class__, self.__class__))

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
    def ebands_list(self) -> List[ElectronBands]:
        """"List of |ElectronBands| objects."""
        return list(self.ebands_dict.values())

    @property
    def edoses_list(self) -> List[ElectronDos]:
        """"List of |ElectronDos| objects."""
        return list(self.edoses_dict.values())

    def iter_lineopt(self):
        """Generates matplotlib linestyles."""
        for o in itertools.product(self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def iter_lineopt_plotly(self):
        """Generates plotly linestyles."""
        for o in itertools.product(self._LINE_WIDTHS,  self._LINE_STYLES_PLOTLY, self._LINE_COLORS):
            yield {"line_width": o[0], "line_dash": o[1], "line_color": o[2]}

    def add_ebands(self, label, bands, edos=None, edos_kwargs=None) -> None:
        """
        Adds a band structure and optionally an edos to the plotter.

        Args:
            label: label for the bands. Must be unique.
            bands: |ElectronBands| object.
            edos: |ElectronDos| object.
            edos_kwargs: optional dictionary with the options passed to ``get_edos`` to compute the electron DOS.
                Used only if ``edos`` is not None and it's not an |ElectronDos| instance.
        """
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
        for meth_name in ("gridplot", "boxplot"):
            yield getattr(self, meth_name)(show=False)

    def yield_plotly_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of plotly figures with minimal input from the user.
        """
        for meth_name in ("gridplot", "boxplot"):
            yield getattr(self, meth_name)(show=False)

    @add_fig_kwargs
    def combiplot(self, e0="fermie", ylims=None, width_ratios=(2, 1), fontsize=8,
                  linestyle_dict=None, **kwargs):
        """
        Plot the band structure and the DOS on the same figure with matplotlib.
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
            fontsize: fontsize for legend.
            linestyle_dict: Dictionary mapping labels to matplotlib linestyle options.

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
            if linestyle_dict is not None and label in linestyle_dict:
                my_kwargs.update(linestyle_dict[label])
            else:
                my_kwargs.update(lineopt)

            opts_label[label] = my_kwargs.copy()

            # Get energy zero.
            mye0 = self.edoses_dict[label].fermie if e0 == "edos_fermie" else ebands.get_e0(e0)

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

    # An alias for combiplot.
    plot = combiplot

    @add_plotly_fig_kwargs
    def combiplotly(self, e0="fermie", ylims=None, width_ratios=(2, 1), fontsize=12,
                    linestyle_dict=None, **kwargs):
        """
        Plot the band structure and the DOS on the same figure with plotly.
        Use ``gridplotly`` to plot band structures on different figures.

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
            width_ratios: Defines the ratio between the band structure plot and the dos plot.
                Used when there are DOS stored in the plotter.
            fontsize: fontsize for titles and legend.
            linestyle_dict: Dictionary mapping labels to linestyle options passed to |plotly.graph_objects.scatter|.

        Returns: |plotly.graph_objects.Figure|
        """
        if self.edoses_dict:
            # Bands and DOS will share the y-axis
            nrows, ncols = (1, 2)
            fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=[], sharex=False, sharey=True,
                                      horizontal_spacing=0.02, column_widths=width_ratios)
        else:
            nrows, ncols = (1, 1)
            fig, _ = get_fig_plotly()

        plotly_set_lims(fig, ylims, 'y')

        # Plot ebands.
        my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        nkpt_list = [ebands.nkpt for ebands in self.ebands_dict.values()]
        if any(nk != nkpt_list[0] for nk in nkpt_list):
            cprint("WARNING: Bands have different number of k-points:\n%s" % str(nkpt_list), "yellow")

        for (label, ebands), lineopt in zip(self.ebands_dict.items(), self.iter_lineopt_plotly()):
            i += 1
            if linestyle_dict is not None and label in linestyle_dict:
                my_kwargs.update(linestyle_dict[label])
            else:
                my_kwargs.update(lineopt)

            opts_label[label] = my_kwargs.copy()

            # Get energy zero.
            mye0 = self.edoses_dict[label].fermie if e0 == "edos_fermie" else ebands.get_e0(e0)

            # Use relative paths if label is a file.
            if os.path.isfile(label): label = os.path.relpath(label)

            rcd = PlotlyRowColDesc(0, 0, nrows, ncols)
            ebands.plotly_traces(fig, mye0, spin=None, band=None, rcd=rcd, showlegend=True, label=label, **my_kwargs)

            # Set ticks and labels, legends.
            if i == 0:
                ebands.decorate_plotly(fig, iax=rcd.iax)

        fig.layout.legend.font.size = fontsize
        fig.layout.title.font.size = fontsize

        # Add DOSes
        if self.edoses_dict:
            rcd = PlotlyRowColDesc(0, 1, nrows, ncols)
            for label, edos in self.edoses_dict.items():
                print(label)
                ebands = self.edoses_dict[label]
                mye0 = ebands.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
                edos.plotly_traces(fig, mye0, rcd=rcd, exchange_xy=True, trace_name=label, **opts_label[label])

        return fig

    plotly = combiplotly

    @add_fig_kwargs
    def gridplot(self, e0="fermie", with_dos=True, with_gaps=False, max_phfreq=None,
                 ylims=None, fontsize=8, **kwargs):
        """
        Plot multiple electron bandstructures and optionally DOSes on a grid with matplotlib.

        Args:
            eb_objects: List of objects from which the band structures are extracted.
                Each item in eb_objects is either a string with the path of the netcdf file,
                or one of the abipy object with an ``ebands`` attribute or a |ElectronBands| object.
            edos_objects: List of objects from which the electron DOSes are extracted.
                Accept filepaths or |ElectronDos| objects. If edos_objects is not None,
                each subplot in the grid contains a band structure with DOS else a simple band structure plot.
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
            with_gaps: True to add markers and arrows showing the fundamental and the direct gap.
            max_phfreq: Max phonon frequency in eV to activate scatterplot showing
                possible phonon absorptions/emission processes based on energy-conservation alone.
                All final states whose energy is within +- max_phfreq of the initial state are included.
                By default, the four electronic states defining the fundamental and the direct gaps
                are considered as initial state (not available for metals).
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ```(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: fontsize for subtitles.

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
                ebands.plot(ax=ax, e0=e0, with_gaps=with_gaps, max_phfreq=max_phfreq, fontsize=fontsize, show=False)
                set_axlims(ax, ylims, "y")
                # This to handle with_gaps = True
                title = ax.get_title()
                if not title: ax.set_title(titles[i], fontsize=fontsize)
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
                ebands.plot_with_edos(edos, e0=mye0, ax_list=(ax0, ax1), with_gaps=with_gaps,
                                      max_phfreq=max_phfreq, show=False)

                # This to handle with_gaps = True
                title = ax0.get_title()
                if not title: ax0.set_title(titles[i], fontsize=fontsize)
                if i % ncols != 0:
                    for ax in (ax0, ax1):
                        ax.set_ylabel("")

        return fig

    @add_plotly_fig_kwargs
    def gridplotly(self, e0="fermie", with_dos=True, with_gaps=False, max_phfreq=None,
                 ylims=None, fontsize=12, **kwargs):
        """
        Plot multiple electron bandstructures and optionally DOSes on a grid with plotly.

        Args:
            eb_objects: List of objects from which the band structures are extracted.
                Each item in eb_objects is either a string with the path of the netcdf file,
                or one of the abipy object with an ``ebands`` attribute or a |ElectronBands| object.
            edos_objects: List of objects from which the electron DOSes are extracted.
                Accept filepaths or |ElectronDos| objects. If edos_objects is not None,
                each subplot in the grid contains a band structure with DOS else a simple band structure plot.
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
            with_gaps: True to add markers and arrows showing the fundamental and the direct gap.
            max_phfreq: Max phonon frequency in eV to activate scatterplot showing
                possible phonon absorptions/emission processes based on energy-conservation alone.
                All final states whose energy is within +- max_phfreq of the initial state are included.
                By default, the four electronic states defining the fundamental and the direct gaps
                are considered as initial state (not available for metals).
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ```(left, right)``
            fontsize: fontsize for subtitles.

        Returns: |plotly.graph_objects.Figure|
        """
        titles = list(self.ebands_dict.keys())
        ebands_list, edos_list = self.ebands_list, self.edoses_list

        nrows, ncols = 1, 1
        numeb = len(ebands_list)

        if not edos_list or not with_dos:
            # Plot grid with bands only.
            if numeb > 1:
                ncols = 2
                nrows = numeb // ncols + numeb % ncols
            fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=list(range(1, numeb+1)),
                                     sharex=False, sharey=True)

            for i, ebands in enumerate(ebands_list):
                irow, icol = divmod(i, ncols)
                band_rcd = PlotlyRowColDesc(irow, icol, nrows, ncols)
                ebands.plotly(e0=e0, fig=fig, rcd=band_rcd, with_gaps=with_gaps, max_phfreq=max_phfreq,
                              fontsize=fontsize, show=False)
                plotly_set_lims(fig, ylims, 'y')
                # This to handle with_gaps = True
                if not with_gaps:
                    fig.layout.annotations[i].text = titles[i]
                    fig.layout.annotations[i].font.size = fontsize
                if (irow, icol) != (0, 0):
                    fig.layout['yaxis%u' % band_rcd.iax].title.text = ""
        else:
            # Special treatment required for phbands with DOS.
            numeb *= 2
            ncols = 4
            nrows = numeb // ncols + numeb % ncols
            # Plot grid with bands + DOS.
            fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=list(range(1, numeb + 1)),
                                     column_widths=[2, 1]*2, horizontal_spacing=0.02, sharex=False, sharey=True)
            # all sub_fig in the same row will share y

            for i, (ebands, edos) in enumerate(zip(ebands_list, edos_list)):
                # Align bands and DOS.
                irow, icol = divmod(i, 2)
                band_rcd = PlotlyRowColDesc(irow, icol * 2, nrows, ncols)
                dos_rcd = PlotlyRowColDesc(irow, icol * 2 + 1, nrows, ncols)
                plotly_set_lims(fig, ylims, 'y', iax=i*2+1)

                # Define the zero of energy and plot
                mye0 = ebands.get_e0(e0) if e0 != "edos_fermie" else edos.fermie
                ebands.plotly_with_edos(edos, fig=fig, band_rcd=band_rcd, dos_rcd=dos_rcd, e0=mye0, with_gaps=with_gaps,
                                        max_phfreq=max_phfreq, show=False)

                # This to handle with_gaps = True
                if not with_gaps:
                    fig.layout.annotations[i * 2].text = titles[i]
                    fig.layout.annotations[i * 2 + 1].text = ''
                    fig.layout.annotations[i * 2].font.size = fontsize
                if i % 2 != 0:
                    fig.layout['yaxis%u' % band_rcd.iax].title.text = ""

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
        df_list = []
        for label, ebands in self.ebands_dict.items():
            # Get the dataframe, select bands and add column with label
            df = ebands.get_dataframe(e0=e0, brange=brange)
            df["label"] = label
            df_list.append(df)
            if ebands.nsppol == 2: spin_polarized = True

        # Merge dataframes ignoring index (not meaningful)
        data = pd.concat(df_list, ignore_index=True)

        import seaborn as sns
        if not spin_polarized:
            ax, fig, plt = get_ax_fig_plt(ax=ax)
            ax.grid(True)
            sns.boxplot(x="band", y="eig", data=data, hue="label", ax=ax, **kwargs)
            if swarm:
                sns.swarmplot(x="band", y="eig", data=data, hue="label", color=".25", ax=ax)
        else:
            # Generate two subplots for spin-up / spin-down channels.
            import matplotlib.pyplot as plt
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

    @add_fig_kwargs
    def plot_band_edges(self, e0="fermie", epad_ev=1.0, set_fermie_to_vbm=True,
                        colormap="viridis", fontsize=8, **kwargs):
        """
        Plot the band edges for electrons and holes on two separated plots for all ebands in ebands_dict.
        Useful for comparing band structures obtained with/without SOC or bands obtained with different settings.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            epad_ev: Add this energy window in eV above VBM and below CBM.
            set_fermie_to_vbm: True if Fermi energy should be recomputed and fixed at max occupied energy level.
            colormap: matplotlib colormap.
            fontsize: legend and title fontsize.
        """
        # Two subplots for CBM and VBM
        num_plots, ncols, nrows = 2, 1, 2
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()
        cmap = plt.get_cmap(colormap)
        nb = len(self.ebands_dict.items())

        for ix, ax in enumerate(ax_list):
            for iband, (label, ebands) in enumerate(self.ebands_dict.items()):
                if set_fermie_to_vbm:
                    # This is needed when the fermi energy is computed in the GS part
                    # with a mesh that does not contain the band edges.
                    ebands.set_fermie_to_vbm()

                if ix == 0:
                    # Conduction
                    ymin = min((ebands.lumos[spin].eig for spin in ebands.spins)) - 0.1
                    ymax = ymin + epad_ev
                elif ix == 1:
                    # Valence
                    ymax = max((ebands.homos[spin].eig for spin in ebands.spins)) + 0.1
                    ymin = ymax - epad_ev
                else:
                    raise ValueError("Wrong ix: %s" % ix)

                # Define ylims and energy shift.
                this_e0 = ebands.get_e0(e0)
                ylims = (ymin - this_e0, ymax - this_e0)
                ebands.plot(ax=ax, e0=e0, color=cmap(float(iband) / nb), ylims=ylims,
                            label=label if ix == 0 else None, show=False)
            if ix == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

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
            nbv.new_code_cell("df = plotter.get_ebands_frame()\ndisplay(frame)"),
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
    def read_ebands(self) -> ElectronBands:
        """
        Returns an instance of |ElectronBands|. Main entry point for client code
        """
        ebands = ElectronBands(
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

        return ebands

    def read_nband_sk(self):
        """|numpy-array| with the number of bands indexed by [s, k]."""
        return self.read_value("number_of_states")

    def read_nspinor(self) -> int:
        """Number of spinors."""
        return self.read_dimvalue("number_of_spinor_components")

    def read_nsppol(self) -> int:
        """Number of independent spins (collinear case)."""
        return self.read_dimvalue("number_of_spins")

    def read_nspden(self) -> int:
        """Number of spin-density components"""
        # FIXME: default 1 is needed for SIGRES files (abinit8)
        return self.read_dimvalue("number_of_components", default=1)

    def read_tsmear(self) -> float:
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

    def read_nelect(self) -> float:
        """
        Number of valence electrons. Note that it's a float because we may have added
        extra fractional charge to the unit cell with a compensanting background.
        """
        return self.read_value("number_of_electrons")

    def read_smearing(self) -> Smearing:
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

    def __init__(self, mesh, spin_dos, nelect, fermie=None, spin_idos=None):
        """
        Args:
            mesh: array-like object with the mesh points in eV.
            spin_dos: array-like object with the DOS for the different spins (even if spin-unpolarized calculation).
                Shape is:
                      (1, nw) if spin-unpolarized.
                      (2, nw) if spin-polarized.
            nelect: Number of electrons in the unit cell.
            fermie: Fermi level in eV. If None, fermie is obtained from the idos integral.
            spin_idos: array-like object with the IDOS for the different spins (even if spin-unpolarized calculation).
                Shape is:
                      (1, nw) if spin-unpolarized.
                      (2, nw) if spin-polarized case.

                This argument is usually used when we have an IDOS computed with a more accurate method e.g.
                tetrahedron integration so that we can use these values instead of integrating the input DOS.

        .. note::

            mesh is given in eV, spin_dos is in states/eV.
        """
        spin_dos = np.atleast_2d(spin_dos)
        self.nsppol = len(spin_dos)
        self.nelect = nelect
        if spin_idos is not None:
            spin_idos = np.atleast_2d(spin_idos)
            assert len(spin_idos) == self.nsppol

        # Save DOS and IDOS for each spin.
        sumv = np.zeros(len(mesh))
        self.spin_dos, self.spin_idos = [], []
        for ispin, values in enumerate(spin_dos):
            sumv += values
            f = Function1D(mesh, values)
            self.spin_dos.append(f)
            # Compute IDOS or take it from spin_idos.
            if spin_idos is None:
                self.spin_idos.append(f.integral())
            else:
                self.spin_idos.append(Function1D(mesh, spin_idos[ispin]))

        # Total DOS and IDOS.
        if self.nsppol == 1: sumv = 2 * sumv
        self.tot_dos = Function1D(mesh, sumv)
        if spin_idos is None:
            # Compute IDOS from DOS
            self.tot_idos = self.tot_dos.integral()
        else:
            # Get IDOS from input (e.g. tetra)
            if self.nsppol == 1: sumv = 2 * spin_idos[0]
            if self.nsppol == 2: sumv = spin_idos[0] + spin_idos[1]
            self.tot_idos = Function1D(mesh, sumv)

        if fermie is not None:
            self.fermie = float(fermie)
        else:
            # *Compute* fermie from nelect. Note that this value may differ
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
    def as_edos(cls, obj: Any, edos_kwargs: dict) -> ElectronDos:
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
                if hasattr(abifile, "ebands"):
                    return abifile.ebands.get_edos(**edos_kwargs)
                elif hasattr(abifile, "edos"):
                    # This to handle e.g. the _EDOS file.
                    return abifile.edos
                else:
                    raise TypeError("Don't know how to extract ElectronDos object from: `%s`" % str(obj))

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
                i = len(idos) - 1
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
                except Exception:
                    raise TypeError("Wrong value for e0: %s" % str(e0))
        else:
            # Assume number
            return float(e0)

    def plot_ax(self, ax, e0, spin=None, what="dos", fact=1.0, exchange_xy=False, **kwargs):
        """
        Helper function to plot the DOS data on the matplotlib axis ``ax``.

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

    def plotly_traces(self, fig, e0, spin=None, what="dos", fact=1.0, exchange_xy=False, rcd=None, trace_name='',
                      showlegend=False, line_opts=None, **kwargs):
        """
        Helper function to plot the DOS data on the ``fig`` with plotly.

        Args:
            fig: |plotly.graph_objects.Figure|.
            e0: Option used to define the zero of energy in the band structure plot.
            spin: selects the spin component, None for total DOS, IDOS.
            what: string selecting what will be plotted. "dos" for DOS, "idos" for IDOS
            fact: Multiplication factor for DOS/IDOS. Usually +-1 for spin DOS
            exchange_xy: True to exchange x-y axis.
            rcd: PlotlyRowColDesc object used to specify the (row, col) of the subplot in the grid.
            trace_name: Name of the trace.
            showlegend: Determines whether or not an item corresponding to this trace is shown in the legend.
            line_opts: Dict of linestyle options passed to |plotly.graph_objects.scatter.Line|
            kwargs: Options passed to fig.add_scatter method.
        """
        dosf, idosf = self.dos_idos(spin=spin)
        e0 = self.get_e0(e0)

        w2f = {"dos": dosf, "idos": idosf}
        if what not in w2f:
            raise ValueError("Unknown value for what: `%s`" % str(what))
        f = w2f[what]

        xx, yy = f.mesh - e0, f.values * fact
        if exchange_xy: xx, yy = yy, xx
        rcd = PlotlyRowColDesc.from_object(rcd)
        ply_row, ply_col = rcd.ply_row, rcd.ply_col
        fig.add_scatter(x=xx, y=yy, mode='lines', name=trace_name, showlegend=showlegend, legendgroup=trace_name,
                        line=line_opts, **kwargs, row=ply_row, col=ply_col)

    @add_fig_kwargs
    def plot(self, e0="fermie", spin=None, ax=None, exchange_xy=False, xlims=None, ylims=None, **kwargs):
        """
        Plot electronic DOS with matplotlib.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                - Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                - None: Don't shift energies, equivalent to ``e0 = 0``.
            spin: Selects the spin component, None if total DOS is wanted.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            exchange_xy: True to exchange x-y axis.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            ylims: Set data limits for the y-axis.
            kwargs: options passed to ``ax.plot``.

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        e0 = self.get_e0(e0)

        for spin in range(self.nsppol):
            opts = {"color": "black", "linewidth": 1.0} if spin == 0 else \
                   {"color": "red", "linewidth": 1.0}
            opts.update(kwargs)
            spin_sign = +1 if spin == 0 else -1
            x, y = self.spin_dos[spin].mesh - e0, spin_sign * self.spin_dos[spin].values
            if exchange_xy: x, y = y, x
            ax.plot(x, y, **opts)

        ax.grid(True)
        xlabel, ylabel = 'Energy (eV)', 'DOS (states/eV)'
        set_ax_xylabels(ax, xlabel, ylabel, exchange_xy)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")

        return fig

    @add_plotly_fig_kwargs
    def plotly(self, e0="fermie", spin=None, fig=None, exchange_xy=False, xlims=None, ylims=None,
               trace_name='', showlegend=False, **kwargs):
        """
        Plot electronic DOS with plotly.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                - Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                - None: Don't shift energies, equivalent to ``e0 = 0``.
            spin: Selects the spin component, None if total DOS is wanted.
            fig: |plotly.graph_objects.Figure|.
            exchange_xy: True to exchange x-y axis.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``.
            trace_name: Name of the trace.
            showlegend: Determines whether or not an item corresponding to this trace is shown in the legend.
            kwargs: Options passed to fig.add_scatter method.

        Return: |plotly.graph_objects.Figure|
        """
        fig, _ = get_fig_plotly(fig=fig)
        e0 = self.get_e0(e0)

        for spin in range(self.nsppol):
            opts = {"color": "black", "width": 1.0} if spin == 0 else \
                   {"color": "red", "width": 1.0}
            # opts.update(kwargs)
            spin_sign = +1 if spin == 0 else -1
            x, y = self.spin_dos[spin].mesh - e0, spin_sign * self.spin_dos[spin].values
            if exchange_xy: x, y = y, x
            fig.add_scatter(x=x, y=y, mode='lines', name=trace_name, showlegend=showlegend, line=opts, **kwargs)

        xlabel, ylabel = 'Energy (eV)', 'DOS (states/eV)'
        if exchange_xy: xlabel, ylabel = ylabel, xlabel
        fig.layout.xaxis.title = xlabel
        fig.layout.yaxis.title = ylabel
        plotly_set_lims(fig, xlims, "x")
        plotly_set_lims(fig, ylims, "y")

        return fig

    @add_fig_kwargs
    def plot_dos_idos(self, e0="fermie", ax_list=None, xlims=None, height_ratios=(1, 2), **kwargs):
        """
        Plot electronic DOS and Integrated DOS on two different subplots with matplotlib.

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
            opts = {"color": "black", "linewidth": 1.0} if spin == 0 else \
                   {"color": "red", "linewidth": 1.0}
            opts.update(kwargs)
            # Plot Total dos if unpolarized.
            if self.nsppol == 1: spin = None
            self.plot_ax(ax_list[0], e0, spin=spin, what="idos", **opts)
            self.plot_ax(ax_list[1], e0, spin=spin, what="dos", **opts)

        return fig

    @add_plotly_fig_kwargs
    def plotly_dos_idos(self, e0="fermie", fig=None, xlims=None, height_ratios=(1, 2), **kwargs):
        """
        Plot electronic DOS and Integrated DOS on two different subplots with plotly.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``.
            fig: The |plotly.graph_objects.Figure| with two distinct plots for the DOS and IDOS plot.
                If fig is None, a new figure is created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``.
            height_ratios:
            kwargs: Options passed to plotly_traces method.

        Return: |plotly.graph_objects.Figure|
        """
        if fig is None:
            fig, _ = get_figs_plotly(nrows=2, ncols=1, sharex=True, sharey=False,
                                     vertical_spacing=0.05, row_heights=height_ratios)
            plotly_set_lims(fig, xlims, "x")
            fig.layout['yaxis1'].title = {'text': "TOT IDOS"}
            fig.layout['yaxis2'].title = {'text': "TOT DOS"}
            fig.layout['xaxis2'].title = {'text': 'Energy (eV)'}

        for spin in range(self.nsppol):
            opts = {"color": "black", "width": 1.0} if spin == 0 else \
                   {"color": "red", "width": 1.0}
            # Plot Total dos if unpolarized.
            if self.nsppol == 1: spin = None
            rcd = PlotlyRowColDesc(0, 0, 2, 1)
            self.plotly_traces(fig, e0, spin=spin, what="idos", rcd=rcd, line_opts=opts, **kwargs)
            rcd = PlotlyRowColDesc(1, 0, 2, 1)
            self.plotly_traces(fig, e0, spin=spin, what="dos", rcd=rcd,  line_opts=opts, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_up_minus_down(self, e0="fermie", ax=None, xlims=None, **kwargs):
        """
        Plot Dos_up - Dow_down with matplotlib.

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

    @add_plotly_fig_kwargs
    def plotly_up_minus_down(self, e0="fermie", fig=None, xlims=None, **kwargs):
        """
        Plot Dos_up - Dow_down with plotly.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            fig: |plotly.graph_objects.Figure| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``.
            kwargs: Options passed to fig.add_scatter method.

        Return: |plotly.graph_objects.Figure|
        """
        dos_diff = self.up_minus_down
        idos_diff = dos_diff.integral()

        e0 = self.get_e0(e0)
        if not kwargs:
            line_opts = {"color": "black", "width": 1.0}

        fig, _ = get_fig_plotly(fig=fig)
        fig.add_scatter(x=dos_diff.mesh - e0, y=dos_diff.values, mode='lines', name='DOS',showlegend=False,
                        line=line_opts, **kwargs)
        fig.add_scatter(x=idos_diff.mesh - e0, y=idos_diff.values, mode='lines', name='IDOS',showlegend=False,
                        line=line_opts, **kwargs)

        plotly_set_lims(fig, xlims, "x")
        fig.layout.xaxis.title = 'Energy (eV)'
        fig.layout.yaxis.title = 'Dos_up - Dos_down (states/eV)'

        return fig

    def to_pymatgen(self):
        """
        Return a pymatgen DOS object from an Abipy |ElectronDos| object.
        """

        from pymatgen.electronic_structure.dos import Dos
        den = {s: d.values for d, s in zip(self.spin_dos, [PmgSpin.up, PmgSpin.down])}
        pmg_dos = Dos(energies=self.spin_dos[0].mesh, densities=den, efermi=self.fermie)

        return pmg_dos


class ElectronDosPlotter(NotebookWriter):
    """
    Class for plotting multiple electronic DOSes.

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
    def edos_list(self) -> List[ElectronDos]:
        """List of DOSes"""
        return list(self.edoses_dict.values())

    def add_edos(self, label, edos, edos_kwargs=None) -> None:
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
    def combiplot(self, what_list="dos", spin_mode="automatic", e0="fermie",
                  ax_list=None,  xlims=None, fontsize=8, **kwargs):
        """
        Plot the the DOSes on the same figure with matplotlib.
        Use ``gridplot`` to plot DOSes on different figures with matplotlib.

        Args:
            what_list: Selects quantities to plot e.g. ["dos", "idos"] to plot DOS and integrated DOS.
                "dos" for DOS only and "idos" for IDOS only
            spin_mode: "total" for total (I)DOS, "resolved" for plotting individual contributions.
                Meaningful only if nsppol == 2.
                "automatic" to use "resolved" if at least one DOS is polarized.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize (int): fontsize for legend

        Return: |matplotlib-Figure|
        """
        if spin_mode == "automatic":
            spin_mode = "resolved" if any(edos.nsppol == 2 for edos in self.edoses_dict.values()) else "total"

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
                    lines = None
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

    @add_plotly_fig_kwargs
    def combiplotly(self, what_list="dos", spin_mode="automatic", e0="fermie",
                    ax_list=None,  xlims=None, fontsize=12, **kwargs):
        """
        Plot the the DOSes on the same figure with plotly.
        Use ``gridplotly`` to plot DOSes on different figures with plotly.

        Args:
            what_list: Selects quantities to plot e.g. ["dos", "idos"] to plot DOS and integrated DOS.
                "dos" for DOS only and "idos" for IDOS only
            spin_mode: "total" for total (I)DOS, "resolved" for plotting individual contributions.
                Meaningful only if nsppol == 2.
                "automatic" to use "resolved" if at least one DOS is polarized.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g ``e0 = 0.5``: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to ``e0 = 0``
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
            fontsize (int): fontsize for legend

        Return: |plotly.graph_objects.Figure|
        """
        if spin_mode == "automatic":
            spin_mode = "resolved" if any(edos.nsppol == 2 for edos in self.edoses_dict.values()) else "total"

        what_list = list_strings(what_list)
        nrows, ncols = len(what_list), 1

        import plotly.colors as pcolors
        # l2color = pcolors.qualitative.Light24
        # ten colors
        l2color = pcolors.DEFAULT_PLOTLY_COLORS

        fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=[], sharex=True, sharey=False)

        can_use_basename = self._can_use_basenames_as_labels()
        for i, what in enumerate(what_list):
            rcd = PlotlyRowColDesc(i, 0, nrows, ncols)
            for j, (label, edos) in enumerate(self.edoses_dict.items()):
                if can_use_basename:
                    label = os.path.basename(label)
                else:
                    # Use relative paths if label is a file.
                    if os.path.isfile(label): label = os.path.relpath(label)

                # Note to have the same color of the same label.
                opt = {"color": l2color[j]}

                # Here I handle spin and spin_mode.
                if edos.nsppol == 1 or spin_mode == "total":
                    # Plot total values
                    if i==0:
                        edos.plotly_traces(fig, e0=e0, what=what, spin=None, rcd=rcd, trace_name=label, showlegend=True,
                                           line_opts=opt)
                    else:
                        edos.plotly_traces(fig, e0=e0, what=what, spin=None, rcd=rcd, trace_name=label, line_opts=opt)
                elif spin_mode == "resolved":
                    # Plot spin resolved quantiies with sign.
                    for spin in range(edos.nsppol):
                        fact = 1 if spin == 0 else -1
                        if (i == 0) & (spin==0) :
                            edos.plotly_traces(fig, e0=e0, what=what, spin=spin, fact=fact, rcd=rcd,
                                           trace_name=label, line_opts=opt, showlegend=True)
                        else:
                            edos.plotly_traces(fig, e0=e0, what=what, spin=spin, fact=fact, rcd=rcd,
                                           trace_name=label, line_opts=opt)
                else:
                    raise ValueError("Wrong value for spin_mode: `%s`:" % str(spin_mode))

            fig.layout.legend.font.size = fontsize
            plotly_set_lims(fig, xlims, "x")
            fig.layout['yaxis'+str(rcd.iax)].title = {'text': 'DOS (states/eV)' if what == "dos" else "IDOS"}
            if i == len(what_list) - 1:
                fig.layout['xaxis' + str(rcd.iax)].title = {'text': 'Energy (eV)'}

        return fig

    # An alias for combiplotly.
    plotly = combiplotly

    @add_fig_kwargs
    def gridplot(self, what="dos", spin_mode="automatic", e0="fermie",
                 sharex=True, sharey=True, xlims=None, fontsize=8, **kwargs):
        """
        Plot multiple DOSes on a grid with matplotlib.

        Args:
            what: "dos" to plot DOS, "idos" for integrated DOS.
            spin_mode: "total" for total (I)DOS, "resolved" for plotting individual contributions.
                Meaningful only if nsppol == 2.
                "automatic" to use "resolved" if at least one DOS is polarized.
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
            fontsize: Axis_label and subtitle fontsize.

        Return: |matplotlib-Figure|
        """
        if spin_mode == "automatic":
            spin_mode = "resolved" if any(edos.nsppol == 2 for edos in self.edoses_dict.values()) else "total"

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
                lines = None
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
            if icol == 0:
                ax.set_ylabel('DOS (states/eV)' if what == "dos" else "IDOS", fontsize=fontsize)
            if irow == nrows - 1:
                ax.set_xlabel("Energy (eV)", fontsize=fontsize)

            #ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_plotly_fig_kwargs
    def gridplotly(self, what="dos", spin_mode="automatic", e0="fermie",
                  sharex=True, sharey=True, xlims=None, fontsize=12, **kwargs):
        """
        Plot multiple DOSes on a grid with plotly.

        Args:
            what: "dos" to plot DOS, "idos" for integrated DOS.
            spin_mode: "total" for total (I)DOS, "resolved" for plotting individual contributions.
                Meaningful only if nsppol == 2.
                "automatic" to use "resolved" if at least one DOS is polarized.
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
            fontsize: Axis_label and subtitle fontsize.

        Return: |plotly.graph_objects.Figure|
        """
        if spin_mode == "automatic":
            spin_mode = "resolved" if any(edos.nsppol == 2 for edos in self.edoses_dict.values()) else "total"

        titles = list(self.edoses_dict.keys())
        edos_list = self.edos_list

        nrows, ncols = 1, 1
        numeb = len(edos_list)
        if numeb > 1:
            ncols = 2
            nrows = numeb // ncols + numeb % ncols

        fig, _ = get_figs_plotly(nrows=nrows, ncols=ncols, subplot_titles=titles, sharex=sharex, sharey=sharey)

        import plotly.colors as pcolors
        # ten colors
        l2color = pcolors.DEFAULT_PLOTLY_COLORS

        for i, (label, edos) in enumerate(self.edoses_dict.items()):
            irow, icol = divmod(i, ncols)
            rcd = PlotlyRowColDesc(irow, icol, nrows, ncols)

            # Here I handle spin and spin_mode.
            if edos.nsppol == 1 or spin_mode == "total":
                opts = {"color": "black", "width": 1.0}
                edos.plotly_traces(fig, e0=e0, what=what, spin=None, rcd=rcd, trace_name=label, line_opts=opts)

            elif spin_mode == "resolved":
                # Plot spin resolved quantiies with sign.
                # Note to have the same color of the same label for both spins.
                opts = {"color": l2color[i], "width": 1.0}
                for spin in range(edos.nsppol):
                    fact = 1 if spin == 0 else -1
                    edos.plotly_traces(fig, e0=e0, what=what, spin=spin, rcd=rcd, trace_name=label, line_opts=opts)
            else:
                raise ValueError("Wrong value for spin_mode: `%s`:" % str(spin_mode))

            fig.layout.annotations[rcd.iax-1].font.size = fontsize
            plotly_set_lims(fig, xlims, "x")
            if icol == 0:
                fig.layout['yaxis'+str(rcd.iax)].title = {'text': 'DOS (states/eV)' if what == "dos" else "IDOS"
                                                          , "font": {"size" : fontsize}}
            if irow == nrows - 1:
                fig.layout['xaxis' + str(rcd.iax)].title = {'text': 'Energy (eV)', "font": {"size" : fontsize}}

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

    #def yield_plotly_figs(self, **kwargs):  # pragma: no cover
    #    """
    #    This function *generates* a predefined list of plotly figures with minimal input from the user.
    #    """
    #    yield self.combiplotly(show=False)
    #    yield self.gridplotly(show=False)

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
            raise ValueError("Wrong inshape: %s" % str(inshape))

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
        print("Producing BXSF file in:", tmp_filepath)
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
    def plot_isosurfaces(self, e0="fermie", cmap=None, verbose=0, **kwargs):
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
            from skimage.measure import marching_cubes_lewiner as marching_cubes
        except ImportError:
            try:
                from skimage.measure import marching_cubes
            except ImportError:
                raise ImportError("scikit-image not installed.\n"
                    "Please install with it with `conda install scikit-image` or `pip install scikit-image`")

        e0 = self.get_e0(e0)
        isobands = self.get_isobands(e0)
        if isobands is None: return None
        if verbose: print("Bands for isosurface:", isobands)

        #from pymatgen.electronic_structure.plotter import plot_lattice_vectors, plot_wigner_seitz
        ax, fig, plt = get_ax3d_fig_plt(ax=None)
        plot_unit_cell(self.reciprocal_lattice, ax=ax, color="k", linewidth=1)
        #plot_wigner_seitz(self.reciprocal_lattice, ax=ax, color="k", linewidth=1)

        for spin in self.spins:
            for ib, band in enumerate(isobands[spin]):
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

                if cmap is not None:
                    cmap = plt.get_cmap(cmap)
                    kwargs["color"] = cmap(float(ib) / len(isobands[spin]))

                ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], **kwargs)
                    #, cmap='Spectral', lw=1, antialiased=True)

                # mayavi package:
                #mlab.triangular_mesh([v[0] for v in verts], [v[1] for v in verts], [v[2] for v in verts], faces)
                #, color=(0, 0, 0))

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
        Contour plot with matplotlib.

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
    #def make_fermisurfer_dir(self, workdir)

#class PhononBands3D(Bands3D):
#    pass


class RobotWithEbands(object):
    """
    Mixin class for robots associated to files with |ElectronBands|.
    """
    def combiplot_ebands(self, **kwargs):
        """Wraps combiplot method of |ElectronBandsPlotter|. kwargs passed to combiplot."""
        return self.get_ebands_plotter().combiplot(**kwargs)

    def combiplotly_ebands(self, **kwargs):
        """Wraps combiplotly method of |ElectronBandsPlotter|. kwargs passed to combiplotly."""
        return self.get_ebands_plotter().combiplotly(**kwargs)

    def gridplot_ebands(self, **kwargs):
        """Wraps gridplot method of |ElectronBandsPlotter|. kwargs passed to gridplot."""
        return self.get_ebands_plotter().gridplot(**kwargs)

    def gridplotly_ebands(self, **kwargs):
        """Wraps gridplotly method of |ElectronBandsPlotter|. kwargs passed to gridplotly."""
        return self.get_ebands_plotter().gridplotly(**kwargs)

    def boxplot_ebands(self, **kwargs):
        """Wraps boxplot method of |ElectronBandsPlotter|. kwargs passed to boxplot."""
        return self.get_ebands_plotter().boxplot(**kwargs)

    def combiboxplot_ebands(self, **kwargs):
        """Wraps combiboxplot method of |ElectronDosPlotter|. kwargs passed to combiboxplot."""
        return self.get_ebands_plotter().combiboxplot(**kwargs)

    def combiplot_edos(self, **kwargs):
        """Wraps combiplot method of |ElectronDosPlotter|. kwargs passed to combiplot."""
        return self.get_edos_plotter().combiplot(**kwargs)

    #def combiplotly_edos(self, **kwargs):
    #    """Wraps combiplotly method of |ElectronDosPlotter|. kwargs passed to combiplotly."""
    #    return self.get_edos_plotter().combiplotly(**kwargs)

    def gridplot_edos(self, **kwargs):
        """Wraps gridplot method of |ElectronDosPlotter|. kwargs passed to gridplot."""
        return self.get_edos_plotter().gridplot(**kwargs)

    #def gridplotly_edos(self, **kwargs):
    #    """Wraps gridplotly method of |ElectronDosPlotter|. kwargs passed to gridplotly."""
    #    return self.get_edos_plotter().gridplotly(**kwargs)

    def get_ebands_plotter(self, kselect=None,
                           filter_abifile=None, cls=None) -> ElectronBandsPlotter:
        """
        Build and return an instance of |ElectronBandsPlotter| or a subclass if ``cls`` is not None.

        Args:
            kselect (str): Used to select particula `ebands`.
                "path" to select bands given on a k-path, "ibz" for bands with IBZ sampling.
                None has not effect
            filter_abifile: Function that receives an ``abifile`` object and returns
                True if the file should be added to the plotter.
            cls: subclass of |ElectronBandsPlotter|.
        """
        plotter = ElectronBandsPlotter() if cls is None else cls()

        for label, abifile in self.items():
            if filter_abifile is not None and not filter_abifile(abifile): continue
            if kselect is not None:
                if kselect == "path" and not abifile.ebands.kpoints.is_path: continue
                if kselect == "ibz" and not abifile.ebands.kpoints.is_ibz: continue
            plotter.add_ebands(label, abifile.ebands)

        return plotter

    def get_edos_plotter(self, cls=None, filter_abifile=None, **kwargs) -> ElectronDosPlotter:
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


from abipy.core.mixins import TextFile #, AbinitNcFile, NotebookWriter
from abipy.abio.robots import Robot


def find_yaml_section_in_lines(lines, tag):

    magic = f"--- !{tag}"
    in_doc, buf = False, []

    for line in lines:
        if line.startswith("#"):
            for i, c in enumerate(line):
                if c != "#": break
            line = line[i:]

        if line.startswith(magic):
            in_doc = True
            continue

        if in_doc and line.startswith("..."):
            in_doc = False
            break

        if in_doc:
            buf.append(line.strip())

    if not buf:
        raise ValueError(f"Cannot fine Yaml tag: `{magic}`")

    import ruamel.yaml as yaml
    return yaml.YAML(typ='safe', pure=True).load("\n".join(buf))


class EdosFile(TextFile):
    """
    This object provides an interface to the _EDOS file
    (electron DOS usually computed with the tetrahedron method).
    The EdosFile has an ElectronDos edos object.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EdosFile
    """

    def __init__(self, filepath):
        """
        Parses the EDOS file and construct self.edos object."""

        super().__init__(filepath)

        # Fortran implementation (eV units). See edos_write in m_ebands.F90.
        #
        # write(unt,"(a)")"# Energy           DOS_TOT          IDOS_TOT         DOS[spin=UP]     IDOS[spin=UP] ..."
        # do iw=1,edos%nw
        #   write(unt,'(es17.8)',advance='no')(edos%mesh(iw) - efermi) * cfact
        #   do spin=0,edos%nsppol
        #     write(unt,'(2es17.8)',advance='no')max(edos%dos(iw,spin) / cfact, tol30), max(edos%idos(iw,spin), tol30)
        #   end do
        #   write(unt,*)
        # end do

        header = []
        data = []

        for line in self:
            if line.startswith("#"):
                header.append(line)
            else:
                line = line.strip()
                if not line: continue
                data.append([float(v) for v in line.split()])

        self.header_string = "".join(header)
        self.edos_params = find_yaml_section_in_lines(header, "EDOS_PARAMS")
        #print(self.edos_params)
        nelect = float(self.edos_params["nelect"])
        data = np.array(data).T.copy()
        mesh = data[0]

        if len(data) == 5:
            # Spin unpolarized case.
            spin_dos = data[3]
            spin_idos = data[4]
        elif len(data) == 7:
            # Spin unpolarized case.
            spin_dos = data[[3, 5]]
            spin_idos = data[[4, 6]]
        else:
            raise ValueError("Don't know how to interpret %d columns in %s" % (len(data), filepath))

        #print(mesh.shape, spin_dos.shape, spin_idos.shape)
        self.edos = ElectronDos(mesh, spin_dos, nelect, fermie=None, spin_idos=spin_idos)

    def to_string(self, verbose=0):
        """String representation."""
        lines = [self.header_string]; app = lines.append
        app(self.edos.to_string(verbose=verbose))

        return "\n".join(lines)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.edos.plot(show=False)
        yield self.edos.plot_dos_idos(show=False)
        if self.edos.nsppol == 2:
            yield self.edos.plot_up_minus_down(show=False)


class EdosRobot(Robot):
    """
    This robot analyzes the results contained in multiple EDOS files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EdosRobot
    """
    EXT = "EDOS"

    #def yield_figs(self, **kwargs):  # pragma: no cover
    #    """
    #    This function *generates* a predefined list of matplotlib figures with minimal input from the user.
    #    Used in abiview.py to get a quick look at the results.
    #    """
    #    yield self.plot_lattice_convergence(show=False)
    #    yield self.plot_gsr_convergence(show=False)
    #    for fig in self.get_ebands_plotter().yield_figs(): yield fig
