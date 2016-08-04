# coding: utf-8
"""Classes for the analysis of electronic structures."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import tempfile
import copy
import itertools
import json
import numpy as np
import pymatgen.core.units as units

from collections import OrderedDict, namedtuple, Iterable
from monty.string import is_string
from monty.json import MSONable, MontyEncoder
from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.bisect import find_le, find_gt
from monty.dev import deprecated
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from pymatgen.serializers.json_coders import pmg_serialize
from abipy.core.func1d import Function1D
from abipy.core.mixins import NotebookWriter
from abipy.core.kpoints import Kpoint, KpointList, Kpath, IrredZone, KpointsReaderMixin, kmesh_from_mpdivs
from abipy.core.structure import Structure
from abipy.iotools import ETSF_Reader, Visualizer, bxsf_write
from abipy.tools import gaussian


import logging
logger = logging.getLogger(__name__)


__all__ = [
    "ElectronBands",
    "ElectronDos",
    "frame_from_ebands",
    "ElectronBandsPlotter",
    "ElectronDosPlotter",
]


class Electron(namedtuple("Electron", "spin kpoint band eig occ")):
    """
    Single-particle state.

    .. Attributes:

        spin: spin index (C convention, i.e >= 0)
        kpoint: :class:`Kpoint` object.
        band: band index. (C convention, i.e >= 0).
        eig: KS eigenvalue.
        occ: Occupation factor.

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

    #def __str__(self):

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

    def asdict(self):
        """Convert self into a dict."""
        return super(Electron, self)._asdict()

    def to_strdict(self, fmt=None):
        """Ordered dictionary mapping fields --> strings."""
        d = self.asdict()
        for (k, v) in d.items():
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
        app("Energy: %s [Ev]" % self.energy)
        app("Initial state: %s" % str(self.in_state))
        app("Final state: %s" % str(self.out_state))

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
    """Stores data and information about the smearing technique."""
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

    @property
    def has_metallic_scheme(self):
        """True if we are using a metallic scheme for occupancies."""
        return self.occopt in [3,4,5,6,7,8]


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

        fermie: Fermi level in eV. Note that, tf the band structure has been computed
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
        Initialize an instance of :class:`ElectronBands` from a netCDF file.
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
        kd = d["kpoints"].copy()
        kd.pop("@module")
        kpoints_cls = KpointList.subclass_from_name(kd.pop("@class"))
        kpoints = kpoints_cls(**kd)

        return cls(Structure.from_dict(d["structure"]), kpoints,
                   d["eigens"], d["fermie"], d["occfacts"], d["nelect"],
                   nband_sk=d["nband_sk"], smearing=d["smearing"],
                   #markers=None, widths=None
                   )

    @pmg_serialize
    def as_dict(self):
        return dict(
            structure=self.structure.as_dict(),
            kpoints=self.kpoints.as_dict(),
            eigens=self.eigens.tolist(),
            fermie=float(self.fermie),
            occfacts=self.occfacts.tolist(),
            nelect=float(self.nelect),
            nband_sk=self.nband_sk.tolist(),
            smearing=self.smearing.as_dict(),
            #markers=, widths=
        )

    @classmethod
    def as_ebands(cls, obj):
        """
        Return an instance of :class:`ElectronBands` from a generic obj.
        Supports:

            - instances of cls
            - files (string) that can be open with abiopen and that provide an `ebands` attribute.
            - objects providing an `ebands` attribute
        """
        if isinstance(obj, cls):
            return obj
        elif is_string(obj):
            # path?
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

    def __init__(self, structure, kpoints, eigens, fermie, occfacts, nelect,
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
        assert self.occfacts.shape == self.occfacts.shape
        self.nsppol, self.nkpt, self.mband = self.eigens.shape

        if nband_sk is not None:
            self.nband_sk = np.array(nband_sk)
        else:
            self.nband_sk = np.array(self.nsppol*self.nkpt*[self.mband])

        self.kpoints = kpoints
        assert self.nkpt == len(self.kpoints)

        self.smearing = {} if smearing is None else smearing
        self.nelect = float(nelect)
        self.fermie = float(fermie)

        # Use the fermi level as zero of energies
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

    # TODO
    #def __repr__(self):
    #    """String representation (short version)"""

    def __str__(self):
        """String representation"""
        return self.info

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
        #odict.update(self.smearing.as_dict())

        bws = self.bandwidths
        for spin in self.spins:
            odict["bandwidth_spin%d" % spin] = bws[spin]

        fundamental_gaps = self.fundamental_gaps
        for spin in self.spins:
            odict["fundamentalgap_spin%d" % spin] = fundamental_gaps[spin].energy
            #odict["fundamentalgap_spin%d" % spin] = fundamental_gaps[spin].energy

        direct_gaps = self.direct_gaps
        for spin in self.spins:
            odict["directgap_spin%d" % spin] = direct_gaps[spin].energy

        #print(self)
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
        from abipy.tools.plotting_utils import Marker
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

    def kindex(self, kpoint):
        """
        The index of the k-point in the internal list of k-points.
        Accepts: :class:`Kpoint` instance or integer.
        """
        if isinstance(kpoint, int):
            return kpoint
        else:
            return self.kpoints.index(kpoint)

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
            my_bands = list(band)

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
            my_bands = list(band)

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

    def to_pdframe(self, e0="fermie"):
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
            kwargs: Keywork arguments passed to seaborn boxplot.
        """
        # Get the dataframe and select bands
        frame = self.to_pdframe(e0=e0)
        if brange is not None: frame = frame[brange[0] <= frame["band"] < brange[1]]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.grid(True)

        import seaborn as sns
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

        #eigenvals is a dict of energies for spin up and spin down
        #{Spin.up:[][],Spin.down:[][]}, the first index of the array
        #[][] refers to the band and the second to the index of the
        #kpoint. The kpoints are ordered according to the order of the
        #kpoints array. If the band structure is not spin polarized, we
        #only store one data set under Spin.up

        eigenvals = {Spin.up: self.eigens[0,:,:].T.copy().tolist()}

        if self.nsppol == 2:
            eigenvals[Spin.down] = self.eigens[1,:,:].T.copy().tolist()

        # FIXME: is_path does not work since info is missing in the netcdf file.
        #if self.kpoints.is_path:
        #    labels_dict = {k.name: k.frac_coords for k in self.kpoints if k.name is not None}
        #    logger.info("calling pmg BandStructureSymmLine with labels_dict %s" % str(labels_dict))
        #    return BandStructureSymmLine(self.kpoints.frac_coords, eigenvals, self.reciprocal_lattice, self.fermie,
        #                                 labels_dict,
        #                                 coords_are_cartesian=False, structure=self.structure, projections=None)
        #else:
        logger.info("Calling pmg BandStructure")
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
                        occ=self.occfacts[spin, kidx, band]
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

    @property
    def use_metallic_scheme(self):
        """True if we are using a metallic scheme for occupancies."""
        return self.smearing.use_metallic_scheme

    #def is_metal(self, spin)
    #    """True if this spin channel is metallic."""
    #    if not self.use_metallic_scheme: return False
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

    @property
    def bandwidths(self):
        """The bandwidth for each spin channel i.e. the energy difference (homo - lomo)."""
        return [self.homos[spin].eig - self.lomos[spin].eig for spin in self.spins]

    @property
    def fundamental_gaps(self):
        """List of :class:`ElectronTransition` with info on the fundamental gaps for each spin."""
        return [ElectronTransition(self.homos[spin], self.lumos[spin]) for spin in self.spins]

    @property
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

    @property
    def info(self):
        """
        Human-readable string with useful info such as band gaps, position of HOMO, LOMO...
        """
        dir_gaps = self.direct_gaps
        fun_gaps = self.fundamental_gaps
        widths = self.bandwidths
        lomos = self.lomos
        homos = self.homos

        lines = []
        app = lines.append

        app("Electron bands of %s" % self.structure.formula)
        app("Number of electrons: %s, Ferml level: %s [eV]" % (self.nelect, self.fermie))

        def indent(s):
            return "\t" + s.replace("\n", "\n\t")

        for spin in self.spins:
            app(">>> Spin %s" % spin)
            app("Direct gap:\n %s" % indent(str(dir_gaps[spin])))
            app("Fundamental gap:\n %s" % indent(str(fun_gaps[spin])))
            app("Bandwidth: %s [eV]" % widths[spin])
            app("LOMO:\n %s" % indent(str(lomos[spin])))
            app("HOMO:\n %s" % indent(str(homos[spin])))

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
            other: :class:`BandStructure` object.
            axis:  Axis along which the statistical parameters are computed.
                   The default is to compute the parameters of the flattened array.
            numpy_op: Numpy function to apply to the difference of the eigenvalues. The
                      default computes |self.eigens - other.eigens|.

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

        def callback(method, step, width):
            edos = self.get_edos(method=method, step=step, width=width)
            edos.plot()

        return ipw.interactive(
                callback,
                method=["gaussian", "tetra"],
                step=ipw.FloatSlider(value=0.1, min=1e-6, max=1, step=0.05, description="Step of linear mesh [eV]"),
                width=ipw.FloatSlider(value=0.2, min=1e-6, max=1, step=0.05, description="Gaussian broadening [eV]"),
            )

    def get_edos(self, method="gaussian", step=0.1, width=0.2, eminmax=None):
        """
        Compute the electronic DOS on a linear mesh.

        Args:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            eminmax: Min and max energy (eV) for the frequency mesh.
                If None, boundaries are automatically computed in order to cover the
                entire energy range.

        Returns:
            :class:`ElectronDos` object.
        """
        # Weights must be normalized to one.
        wsum = self.kpoints.sum_weights()
        if abs(wsum - 1) > 1.e-6:
            err_msg = "Kpoint weights should sum up to one while sum_weights is %.3f\n" % wsum
            err_msg += "The list of kpoints does not represent a homogeneous sampling of the BZ\n"
            err_msg += str(type(self.kpoints)) + "\n" + str(self.kpoints)
            raise ValueError(err_msg)

        # Compute the linear mesh.
        if eminmax is not None:
            e_min, e_max = eminmax
        else:
            e_min = self.enemin()
            e_min -= 0.1 * abs(e_min)
            e_max = self.enemax()
            e_max += 0.1 * abs(e_max)

        nw = int(1 + (e_max - e_min) / step)
        mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)
        dos = np.zeros((self.nsppol, nw))

        if method == "gaussian":
            for spin in self.spins:
                for k, kpoint in enumerate(self.kpoints):
                    weight = kpoint.weight
                    for band in range(self.nband_sk[spin,k]):
                        e = self.eigens[spin,k,band]
                        dos[spin] += weight * gaussian(mesh, width, center=e)

        else:
            raise ValueError("Method %s is not supported" % method)

        return ElectronDos(mesh, dos, self.nelect)

    def get_ejdos(self, spin, valence, conduction, method="gaussian", step=0.1, width=0.2, mesh=None):
        """
        Compute the join density of states at q==0
            :math:`\sum_{kbv} f_{vk} (1 - f_{ck}) \delta(\omega - E_{ck} - E_{vk})`

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
    def plot_ejdosvc(self, vrange, crange, method="gaussian", step=0.1, width=0.2, cumulative=True, **kwargs):
        """
        Plot the decomposition of the joint-density of States (JDOS).

        Args:
            vrange: Int or `Iterable` with the indices of the valence bands to consider.
            crange: Int or `Iterable` with the indices of the conduction bands to consider.
            method: String defining the method.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.
            cumulative: True for cumulative plots (default).

        Returns:
            `matplotlib` figure
        """
        if not isinstance(crange, Iterable): crange = [crange]
        if not isinstance(vrange, Iterable): vrange = [vrange]

        import matplotlib.pyplot as plt
        fig = plt.figure()

        for s in self.spins:
            ax = fig.add_subplot(1, self.nsppol, s+1)

            # Get total JDOS
            tot_jdos = self.get_ejdos(s, vrange, crange, method=method, step=step, width=width)

            jdos_vc = OrderedDict()
            for v in vrange:
                for c in crange:
                    jd = self.get_ejdos(s, v, c, method=method, step=step, width=width, mesh=tot_jdos.mesh)
                    jdos_vc[(v, c)] = jd

            # Plot data for this spin.
            if cumulative:
                cmap = plt.get_cmap("jet")
                cumulative = np.zeros(len(tot_jdos))
                num_plots, i = len(jdos_vc), 0

                for (v, c), jdos in jdos_vc.items():
                    label = "val=%s --> cond=%s, s=%s" % (v, c, s)
                    color = cmap(float(i)/num_plots)
                    x, y = jdos.mesh, jdos.values
                    ax.plot(x, cumulative + y, lw=1.0, label=label, color=color)
                    ax.fill_between(x, cumulative, cumulative + y, facecolor=color, alpha=0.7)
                    cumulative += jdos.values
                    i += 1
                #tot_jdos.plot_ax(ax, color="k", lw=2, label=" Total JDOS: s=%s," % s, **kwargs)

            else:
                tot_jdos.plot_ax(ax, label="Total JDOS: s=%s," % s, **kwargs)
                for (v, c), jdos in jdos_vc.items():
                    jdos.plot_ax(ax, label="val=%s --> cond=%s, s=%s" % (v,c,s), **kwargs)

        plt.legend(loc="best")
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

        return ElectronBands(
            self.structure, self.kpoints, qp_energies, fermie, self.occfacts, self.nelect,
            nband_sk=self.nband_sk, smearing=self.smearing, markers=self.markers)

    @add_fig_kwargs
    def plot(self, ax=None, klabels=None, band_range=None, e0="fermie", marker=None, width=None, **kwargs):
        """
        Plot the band structure.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            klabels: dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels. e.g.
                klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0):"L"}.
            band_range: Tuple specifying the minimum and maximum band to plot (default: all bands are plotted)
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
        self.decorate_ax(ax, klabels=klabels) #, title=title)

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

    @add_fig_kwargs
    def plot_fatbands(self, klabels=None, e0="fermie", **kwargs):
        #colormap="jet", max_stripe_width_mev=3.0, qlabels=None, **kwargs):
        """
        Plot the electronic fatbands.

        Args:
            klabels: Dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. ~klabels = { (0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L" }`.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            band_range: Tuple specifying the minimum and maximum band to plot (default: all bands are plotted)
            width: String defining the width to plot. accepts the syntax widthname:fact where
                fact is a float used to scale the stripe size.

        Returns:
            `matplotlib` figure
        """
        # Build grid of plots.
        num_plots, ncols, nrows = len(self.widths), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()
        # don't show the last ax if numeb is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        for ax, key in zip(ax_list, self.widths):
            # Decorate the axis
            self.decorate_ax(ax, klabels=klabels, title=key)

            # Plot the energies.
            self.plot_ax(ax, e0)
            # Add width around each band.
            self.plot_width_ax(ax, key, e0)

        return fig

    def decorate_ax(self, ax, **kwargs):
        title = kwargs.pop("title", None)
        if title is not None: ax.set_title(title)

        ax.grid(True)
        #ax.set_xlabel('k-point')
        ax.set_ylabel('Energy [eV]')

        # FIXME:
        # This causes the annoying warning
        #UserWarning: No labeled objects found. Use label='...' kwarg on individual plots.
        # perhaps a method ax_finalize?
        #ax.legend(loc="best", shadow=True)

        # Set ticks and labels.
        ticks, labels = self._make_ticks_and_labels(kwargs.pop("klabels", None))
        if ticks:
            ax.set_xticks(ticks, minor=False)
            ax.set_xticklabels(labels, fontdict=None, minor=False)

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
    def plot_with_edos(self, dos, klabels=None, axlist=None, e0="fermie", **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            dos: An instance of :class:`ElectronDos`.
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. `klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}`.
            axlist: The axes for the bandstructure plot and the DOS plot. If axlist is None, a new figure
                is created and the two axes are automatically generated.
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

        # Plot the DOS
        dos.plot_ax(ax2, e0, exchange_xy=True, **kwargs)
        ax2.grid(True)
        ax2.yaxis.set_ticks_position("right")
        ax2.yaxis.set_label_position("right")

        fig = plt.gcf()
        return fig

    def export_bxsf(self, filepath):
        """
        Export the full band structure on filepath in the BXSF format
        suitable for the visualization of the Fermi surface.
        """
        # Sanity check.
        errors = []
        eapp = errors.append

        if np.any(self.nband_sk != self.nband_sk[0,0]):
            eapp("nband must be constant")

        # TODO
        #if not self.has_omesh:
        #    eapp("An homogeneous sampling is needed for the Fermi surface")

        if not np.allclose(self.kpoints.shifts, 0.0):
            eapp("shifted meshes are not supported by Xcryden")

        if errors:
            raise ValueError("\n".join(errors))

        if "." not in filepath:
            raise ValueError("Cannot detect file extension in path %s: " % filepath)

        tokens = filepath.strip().split(".")
        ext = tokens[-1]

        if not tokens[0]:
            # fname == ".ext" ==> Create temporary file.
            path = tempfile.mkstemp(suffix="."+ext, text=True)[1]

        # Xcrysden requires points in the unit cell and the mesh must include the periodic images hence use pbc=True.
        ndivs = self.kpoints.mpdivs

        ibz_arr = self.kpoints.to_array()
        #print("ibz_arr", ibz_arr)

        ebands3d = EBands3D(self.structure,
                            ibz_arr=self.kpoints.to_array(), ene_ibz=self.eigens,
                            ndivs=ndivs, shifts=self.kpoints.shifts,
                            pbc=True, order="unit_cell")

        # Symmetrize bands in the unit cell.
        emesh_sbk = ebands3d.get_emesh_sbk()

        #print(self.nband, ndivs+1)
        with open(filepath, "w") as fh:
            bxsf_write(fh, self.structure, self.nsppol, self.nband, ndivs+1, emesh_sbk, self.fermie, unit="eV")

        return Visualizer.from_file(filepath)

    def derivatives(self, spin, band, order=1, acc=4, asmarker=None):
        """
        Compute the derivative of the eigenvalues along the path.
        """
        if self.has_bzpath:
            # Extract the branch.
            branch = self.eigens[spin,:,band]

            # Simulate free-electron bands. This will produce all(effective masses == 1)
            #branch = [0.5 * Ha_eV * (k.norm * Bohr_Ang)**2 for k in self.kpoints]

            # Compute derivatives by finite differences.
            ders_onlines = self.kpoints.finite_diff(branch, order=order, acc=acc)

            if asmarker is None:
                return ders_onlines
            else:
                x, y, s = [], [], []
                for i, line in enumerate(self.kpoints.lines):
                    #print(line)
                    x.extend(line)
                    y.extend(branch[line])
                    s.extend(ders_onlines[i])
                    assert len(x) == len(y) == len(s)

                self.set_marker(asmarker, (x,y,s))
                return x, y, s

        else:
            raise ValueError("Derivatives on homogeneous k-meshes are not supported yet")

    def effective_masses(self, spin, band, acc=4):
        """Compute the effective masses."""
        ders2 = self.derivatives(spin, band, order=2, acc=acc) * units.eV_to_Ha / units.bohr_to_ang**2
        return 1.0 / ders2


def frame_from_ebands(ebands_objects, index=None, with_spglib=True):
    """
    Build a pandas dataframe with the most important results available in a list of band structures.

    Args:
        struct_objects: List of objects that can be converted to structure.
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
        plotter.add_ebands("foo.nc", label="foo bands")
        plotter.add_ebands("bar.nc", label="bar bands")
        plotter.gridplot()

    Dictionary with the mapping label --> edos.
    """
    _LINE_COLORS = ["b", "r",]
    _LINE_STYLES = ["-",":","--","-.",]
    _LINE_WIDTHS = [2,]

    def __init__(self, key_ebands=None, key_edos=None, edos_kwargs=None):
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
    def combiplot(self, e0="fermie", **kwargs):
        """
        Plot the band structure and the DOS on the same figure.
        Use `gridplot` to plot band structurs on different figures.

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

    @deprecated(message="plot method of ElectronBands has been replaced by combiplot.")
    def plot(self, *args, **kwargs):
        if "align" in kwargs or "xlim" in kwargs or "ylim" in kwargs:
            raise ValueError("align|xlim|ylim options are not supported anymore.")
        return self.combiplot(*args, **kwargs)

    @add_fig_kwargs
    def gridplot(self, e0="fermie", with_dos=True, **kwargs):
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
            titles:
                List of strings with the titles to be added to the subplots.
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
            kwargs: Keywork arguments passed to seaborn boxplot.
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
            frame = ebands.to_pdframe(e0=e0)
            if brange is not None: frame = frame[brange[0] <= frame["band"] < brange[1]]
            frame["label"] = label
            frames.append(frame)
            if ebands.nsppol == 2: spin_polarized = True

        # Merge frames ignoring index (not meaningful)
        import pandas as pd
        data = pd.concat(frames, ignore_index=True)

        import seaborn as sns
        if not spin_polarized:
            ax, fig, plt = get_ax_fig_plt(ax=ax)
            sns.boxplot(x="band", y="eig", data=data, hue="label", ax=ax, **kwargs)
            if swarm:
                sns.swarmplot(x="band", y="eig", data=data, hue="label", color=".25", ax=ax)
        else:
            if ax is not None:
                raise NotImplementedError("ax == None not implemented when nsppol==2")
            fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, squeeze=False)
            for spin, ax in zip(range(2), axes.ravel()):
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
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file is created.
        Return path to the notebook.
        """
        import io, tempfile
        if nbpath is None:
            _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        #args = [(l, f.filepath) for l, f in self.items()]
        #nb.cells.extend([
        #    nbv.new_markdown_cell("# This is a markdown cell"),
        #    nbv.new_code_cell("robot = abilab.GsrRobot(*%s)\nprint(robot)" % str(args)),
        #    nbv.new_code_cell("frame = robot.get_dataframe()\ndisplay(frame)"),
        #    nbv.new_code_cell("plotter = robot.get_ebands_plotter()"),
        #    nbv.new_code_cell("fig = plotter.plot()"),
        #])

        with io.open(nbpath, 'wt', encoding="utf8") as fh:
            nbformat.write(nb, fh)
        return nbpath


class ElectronDosPlotter(object):
    """
    Class for plotting electronic DOSes.

    Usage example:

    .. code-block:: python

        plotter = ElectronDosPlotter()
        plotter.add_edos_from_file("foo.nc", label="foo dos")
        plotter.add_edos_from_file("bar.nc", label="bar dos")
        plotter.plot()
    """
    #_LINE_COLORS = ["b", "r",]
    #_LINE_STYLES = ["-",":","--","-.",]
    #_LINE_WIDTHS = [2,]

    def __init__(self, *args):
        self._edoses_dict = OrderedDict()
        for label, dos in args:
            self.add_edos(label, dos)

    #def iter_lineopt(self):
    #    """Generates style options for lines."""
    #    for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
    #        yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    @property
    def edoses_dict(self):
        """Dictionary with the DOSes"""
        return self._edoses_dict

    @property
    def edoses(self):
        """List of DOSes"""
        return list(self._edoses_dict.values())

    def add_edos_from_file(self, filepath, label=None, method="gaussian", step=0.1, width=0.2):
        """
        Adds a dos for plotting. Reads data from a Netcdf file
        """
        ebands = ElectronBands.as_ebands(filepath)
        edos = ebands.get_edos(method=method, step=step, width=width)
        if label is None: label = filepath
        self.add_edos(label, edos)

    def add_edos(self, label, edos):
        """
        Adds a DOS for plotting.

        Args:
            label: label for the DOS. Must be unique.
            dos: :class:`ElectronDos` object.
        """
        if label in self.edoses_dict:
            raise ValueError("label %s is already in %s" % (label, self._edoses_dict.keys()))

        self._edoses_dict[label] = edos

    @add_fig_kwargs
    def plot(self, ax=None, e0="fermie", **kwargs):
        """
        Plot the the DOSes.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0

        Returns:
            `matplotlib` figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for label, dos in self.edoses_dict.items():
            # Use relative paths if label is a file.
            if os.path.isfile(label): label = os.path.relpath(label)
            dos.plot_ax(ax, e0, label=label)

        ax.grid(True)
        ax.set_xlabel("Energy [eV]")
        ax.set_ylabel('DOS [states/eV]')
        ax.legend(loc="best")

        return fig

    #def animate(self, **kwargs):
    #    animator = Animator()
    #    tmpdir = tempfile.mkdtemp()
    #    for (label, dos) in self.edoses_dict.items():
    #        savefig = os.path.join(tmpdir, label + ".png")
    #        dos.plot(show=False, savefig=savefig)
    #        animator.add_figure(label, savefig)
    #    return animator.animate(**kwargs)


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
        return self.read_dimvalue("number_of_components")

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
        occopt = self.read_value("occopt")

        try:
            scheme = "".join(c for c in self.read_value("smearing_scheme"))
            scheme = scheme.strip()
        except TypeError:
            scheme = None

        # FIXME there's a problem in smearing_scheme
        #if scheme is None:
        #    logger.warning("warning: scheme is None, occopt %s" % occopt)

        return Smearing(
            scheme=scheme,
            occopt=occopt,
            tsmear_ev=units.Energy(self.read_value("smearing_width"), "Ha").to("eV")
        )


class EBands3D(object):
    """This object symmetrizes the band energies in the full Brillouin zone."""
    def __init__(self, structure, ibz_arr, ene_ibz, ndivs, shifts, pbc=False, order="unit_cell"):
        """
        Args:
            structure: :class:`Structure` object.
            ibz_arr: `ndarray` with the k-points in the IBZ in reduced coordinates.
            ene_ibz: Energies in the IBZ. shape = (nsppol, nk_ibz, nband)
            ndivs: Number of divisions used to generate the k-mesh
            shifts: List of shifts used to generate the k-mesh
            pbc: True if periodic boundary conditions should be enfored
                during the generation of the k-mesh
            order: String defining the order of the k-point in the k-mesh.
                "unit_cell":
                    Point are located in the unit_cell, i.e kx in [0, 1]
                "bz":
                    Point are located in the Brillouin zone, i.e kx in [-1/2, 1/2].
        """
        self.ibz_arr = ibz_arr
        self.ene_ibz = np.atleast_3d(ene_ibz)
        self.nsppol, self.nkibz, self.nband = self.ene_ibz.shape

        self.ndivs, self.shifts = ndivs, np.atleast_2d(shifts)
        self.pbc, self.order = pbc, order

        # Compute the full list of k-points according to order.
        self.bz_arr = kmesh_from_mpdivs(self.ndivs, shifts, pbc=pbc, order=order)

        # Compute the mapping bz --> ibz
        from abipy.extensions.klib import map_bz2ibz
        self.bz2ibz = map_bz2ibz(structure, self.bz_arr, self.ibz_arr)

        #for ik_bz, ik_ibz in enumerate(self.bz2ibz):
        #   print(ik_bz, ">>>", ik_ibz)

        if np.any(self.bz2ibz == -1):
            raise ValueError("-1 found")

    @property
    def spins(self):
        """Used to iterate over spin indices."""
        return range(self.nsppol)

    @property
    def bands(self):
        """Used to iterate over bands"""
        return range(self.nband)

    @property
    def len_bz(self):
        """Number of point in the full bz."""
        return len(self.bz_arr)

    def get_emesh_sbk(self):
        """
        Returns a `ndarray` with shape [nsppol, nband, len_bz] with the eigevanalues in the full zone.
        """
        emesh_sbk = np.empty((self.nsppol, self.nband, self.len_bz))

        for spin in self.spins:
            for band in self.bands:
                emesh_sbk[spin,band,:] = self.get_emesh_k(spin, band)

        return emesh_sbk

    def get_emesh_k(self, spin, band):
        """
        Return a `ndarray` with shape [len_bz] with the energies in the full zone for given spin and band.
        """
        emesh_k = np.empty(self.len_bz)

        ene_ibz = self.ene_ibz[spin,:,band]
        # e_{Sk} = e_{k}
        for ik_bz in range(self.len_bz):
           ik_ibz = self.bz2ibz[ik_bz]
           #print(ik_bz, ">>>", ik_ibz)
           emesh_k[ik_bz] = self.ene_ibz[spin,ik_ibz,band]

        return emesh_k

    #def plane_cut(self, values_ibz):
    #    """
    #    Symmetrize values in the IBZ to have them on the full BZ, then
    #    select a slice along the specified plane E.g. plane = (1,1,0).
    #    """
    #    assert len(values_ibz) == len(self)
    #    #indices =
    #    z0 = 0
    #    plane = np.empty((self.nx, self.ny))
    #    kx, ky = range(self.nx), range(self.ny)
    #    for x in kx:
    #        for y in ky:
    #            ibz_idx = self.map_xyz2ibz[x, y, z0]
    #            plane[x, y] = values_ibz[ibz_idx]
    #    kx, ky = np.meshgrid(kx, ky)
    #    return kx, ky, plane


class ElectronDos(object):
    """
    This object stores the electronic density of states.
    It is usually created by calling the get_edos method of :class:`ElectronBands`.
    """

    def __init__(self, mesh, spin_dos, nelect):
        """
        Args:
            mesh: array-like object with the mesh points.
            spin_dos: array-like object with the DOS value for the different spins.
                      spin_dos[1, nw] if spin-unpolarized.
                      spin_dos[2, nw] if spin-polarized case.
            nelect: Number of electrons in the unit cell.

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

        # *Compute* fermie from nelect. Note that this value could differ
        # from the one stored in ElectronBands (coming from the SCF run)
        # The accuracy of self.fermie depends on the number of k-points used for the DOS
        # and the parameters used to call ebands.get_edos.
        self.fermie = self.find_mu(self.nelect)

    def __str__(self):
        lines = []
        app = lines.append
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

    def find_mu(self, nelect, spin=None, num=500, atol=1.e-5):
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

        # Now use spline to get a more accurate mu (useful if mesh is coarse)
        e0, e1 = idos.mesh[i-1], idos.mesh[i]
        idos_spline = idos.spline
        for mu in np.linspace(e0, e1, num=num):
            if abs(idos_spline(mu) - nelect) < atol:
                return mu
        else:
            raise RuntimeError("Cannot find mu, try to increase num and/or atol")

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

    def plot_ax(self, ax, e0, spin=None, what="d", exchange_xy=False, **kwargs):
        """
        Helper function to plot the DOS data on the axis ax.

        Args:
            ax: matplotlib axis.
            e0: Option used to define the zero of energy in the band structure plot.
            spin: selects the spin component, None for total DOS, IDOS.
            what: string selecting what will be plotted:
                  "d" for DOS, "i" for IDOS. chars can be concatenated
                  hence what="id" plots both IDOS and DOS. (default "d").
            exchange_xy: True to exchange axis.
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
            xx, yy = f.mesh - e0, f.values
            if exchange_xy: xx, yy = yy, xx
            lines.extend(ax.plot(xx, yy, **kwargs))
        return lines

    @add_fig_kwargs
    def plot(self, e0="fermie", spin=None, **kwargs):
        """
        Plot DOS and IDOS.

        Args:
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            spin: Selects the spin component, None if total DOS is wanted.

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

        ax2.set_xlabel('Energy [eV]')

        ax1.set_ylabel("TOT IDOS" if spin is None else "IDOS (spin %s)" % spin)
        ax2.set_ylabel("TOT DOS" if spin is None else "DOS (spin %s)" % spin)

        self.plot_ax(ax1, e0, spin=spin, what="i", **kwargs)
        self.plot_ax(ax2, e0, spin=spin, what="d", **kwargs)

        fig = plt.gcf()
        return fig
