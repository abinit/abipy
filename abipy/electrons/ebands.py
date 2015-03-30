# coding: utf-8
"""Classes for the analysis of electronic structures."""
from __future__ import print_function, division, unicode_literals

import sys
import os
import tempfile
import copy
import itertools
import numpy as np
import pymatgen.core.units as units

from collections import OrderedDict, namedtuple, Iterable
from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.bisect import find_le, find_gt
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from abipy.core.func1d import Function1D
from abipy.core.kpoints import Kpoint, Kpath, IrredZone, KpointsReaderMixin, kmesh_from_mpdivs
from abipy.iotools import ETSF_Reader, Visualizer, bxsf_write
from abipy.tools import gaussian
from abipy.tools.animator import FilesAnimator

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "ElectronBands",
    "ElectronBandsPlotter",
    "ElectronDosPlotter",
    "ElectronsReader",
# TODO Rename it, use camel case
    "ElectronDOS",
    "ElectronDOSPlotter",
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
    #    return self.in_state == other.in_state and 
    #           self.out_state == other.out_state

    #def __ne__(self, other):
    #    return not self == other

    #def __ge__(self, other):
    #    return self.energy >=  other.energy

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
    """This object stores the electronic band structure."""
    Error = ElectronBandsError

    @classmethod
    def from_file(cls, filepath):
        """Initialize an instance of :class:`ElectronBands` from a netCDF file."""
        if filepath.endswith(".nc"):
            with ElectronsReader(filepath) as r:
                new = r.read_ebands()
        else:
            raise NotImplementedError("ElectronBands can only be initialized from nc files")

        assert new.__class__ == cls
        return new

    def __init__(self, structure, kpoints, eigens, fermie, occfacts, nelect,
                 nband_sk=None, smearing=None, markers=None, widths=None):
        """
        Args:
            structure: pymatgen structure.
            kpoints: :class:`KpointList` instance.
            eigens: Array-like object with the eigenvalues (eV) stored as [s,k,b]
                    where s: spin , k: kpoint, b: band index
            fermie: Fermi level in eV
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

        # Fix the Fermi level and use efermi as the energy zero.
        self._fix_fermie()

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
    #    return self.info

    def __str__(self):
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

        for (band, e) in enumerate(energies[1:]):
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

    def raw_print(self, stream=sys.stdout):
        """Print k-points and energies on stream."""
        stream.write("# Band structure energies in Ev.\n")
        stream.write("# idx   kpt_red(1:3)  ene(b1) ene(b2) ...\n")

        fmt_k = lambda k: " %.6f" % k
        fmt_e = lambda e: " %.6f" % e
        for spin in self.spins:
            stream.write("# spin = " + str(spin) + "\n")
            for (k, kpoint) in enumerate(self.kpoints):
                nb = self.nband_sk[spin,k]
                ene_sk = self.eigens[spin,k,:nb]
                st = str(k+1)
                for c in kpoint: st += fmt_k(c)
                for e in ene_sk: st += fmt_e(e)
                stream.write(st+"\n")

        stream.flush()

    def to_pymatgen(self, fermie=None):
        """
        Return a pymatgen bandstructure object.

        Args:
            fermie: Fermi energy in eV. If None, self.efermi is used.
        """
        from pymatgen.electronic_structure.core import Spin
        from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine

        assert np.all(self.nband_sk == self.nband_sk[0,0])

        # TODO check this
        if fermie is None: fermie = self.fermie

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
        #    return BandStructureSymmLine(self.kpoints.frac_coords, eigenvals, self.reciprocal_lattice, fermie, labels_dict,
        #                                 coords_are_cartesian=False, structure=self.structure, projections=None)

        #else:
        logger.info("Calling pmg BandStructure")
        return BandStructure(self.kpoints.frac_coords, eigenvals, self.reciprocal_lattice, fermie,
                                 labels_dict=None, coords_are_cartesian=False, structure=self.structure, projections=None)

    def _electron_state(self, spin, kpoint, band):
        """
        Build an instance of :class:`Electron` from the spin, kpoint and band index"""
        kidx = self.kindex(kpoint)
        return Electron(spin=spin,
                        kpoint=self.kpoints[kidx],
                        band=band,
                        eig=self.eigens[spin, kidx, band],
                        occ=self.occfacts[spin, kidx, band])

    #def from_scfrun(self):
    #    return self.iscf > 0

    #def has_occupations(self):
    #    return np.any(self.occfacts != 0.0)

    def _fix_fermie(self):
        """
        Fix the value of the Fermi level in semiconductors, set it to the HOMO level.
        """
        # Use the fermi level computed by Abinit for metals or if SCF run
        # FIXME
        #if self.use_metallic_scheme or self.from_scfrun: return
    
        # FXME This won't work if ferromagnetic semi-conductor.
        try:
            occfact = 2 if self.nsppol == 1 else 1
            esb_levels = []
            for k in self.kidxs:
                esb_view = self.eigens[:,k,:].T.ravel()
                for i, esb in enumerate(esb_view):
                    if (i+1) * occfact == self.nelect:
                        esb_levels.append(esb)
                        break
                else:
                    raise ValueError("Not enough bands to compute the position of the Fermi level!")

        except ValueError:
            return

        new_fermie = max(esb_levels)

        #if abs(new_fermie - self.fermie) > 0.2:
        #    print("old_fermie %s, new fermie %s" % (self.fermie, new_fermie))

        # Use fermilevel as zero of energies.
        self.fermie = new_fermie
        # FIXME this is problematic since other values e.g. QP corrections
        # are expressed in terms of KS energies and I should always have
        # the fermi level to shift energies correctly. perhaps it's better if I provide
        # an option to shift the energies in the plot.
        #self._eigens = self._eigens - new_fermie 
        #self.fermie = 0.0

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
        b = find_le(self.eigens[spin,k,:], self.fermie)
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
        b = find_gt(self.eigens[spin,k,:], self.fermie)
        return self._electron_state(spin, k, b)

    @property
    def homos(self):
        """homo states for each spin channel as a list of nsppol :class:`Electron`."""
        homos = self.nsppol * [None]

        for spin in self.spins:
            blist, enes = [], []
            for k in self.kidxs:
                # Find rightmost value less than or equal to fermie.
                b = find_le(self.eigens[spin,k,:], self.fermie)
                blist.append(b)
                enes.append(self.eigens[spin,k,b])

            homo_kidx = np.array(enes).argmax()
            homo_band = blist[homo_kidx]

            # Build Electron instance.
            homos[spin] = self._electron_state(spin, homo_kidx, homo_band)

        return homos

    @property
    def lumos(self):
        """lumo states for each spin channel as a list of nsppol :class:`Electron`."""
        lumos = self.nsppol * [None]
                                                                     
        for spin in self.spins:
            blist, enes = [], []
            for k in self.kidxs:
                # Find leftmost value greater than x.
                b = find_gt(self.eigens[spin,k,:], self.fermie)
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
        """Human-readable string with useful info such as band gaps, position of HOMO, LOMO..."""
        dir_gaps = self.direct_gaps
        fun_gaps = self.fundamental_gaps
        widths = self.bandwidths
        lomos = self.lomos
        homos = self.homos

        lines = []
        app = lines.append

        app("Electron bands of %s" % self.structure.formula)
        app("Number of electrons %s" % self.nelect)
        app("Fermi level: %s [eV]" % self.fermie)

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

        Args:
            other: :class:`BandStructure` object.
            axis:  Axis along which the statistical parameters are computed.
                   The default is to compute the parameters of the flattened array.
            numpy_op: Numpy function to apply to the difference of the eigenvalues. The
                      default computes t|self.eigens - other.eigens|.

        Returns:
            `namedtuple` with the statistical parameters in eV
        """
        ediff = numpy_op(self.eigens - other.eigens)

        return StatParams(mean=ediff.mean(axis=axis),
                          stdev=ediff.std(axis=axis),
                          min=ediff.min(axis=axis),
                          max=ediff.max(axis=axis)
                          )

    def get_edos(self, method="gaussian", step=0.1, width=0.2):
        """
        Compute the electronic DOS on a linear mesh.

        Args:
            method: String defining the method for the computation of the DOS.
            step: Energy step (eV) of the linear mesh.
            width: Standard deviation (eV) of the gaussian.

        Returns:
            :class:`ElectronDOS` object.
        """
        # Weights must be normalized to one.
        wsum = self.kpoints.sum_weights()
        if abs(wsum - 1) > 1.e-6:
            err_msg = "Kpoint weights should sum up to one while sum_weights is %.3f\n" % wsum
            err_msg += "The list of kpoints does not represent a homogeneous sampling of the BZ\n" 
            err_msg += str(type(self.kpoints)) + "\n" + str(self.kpoints)
            raise ValueError(err_msg)

        # Compute the linear mesh.
        e_min = self.enemin()
        e_min -= 0.1 * abs(e_min)

        e_max = self.enemax()
        e_max += 0.1 * abs(e_max)

        nw = int(1 + (e_max - e_min) / step)
        mesh, step = np.linspace(e_min, e_max, num=nw, endpoint=True, retstep=True)

        dos = np.zeros((self.nsppol, nw))

        if method == "gaussian":
            for spin in self.spins:
                for (k, kpoint) in enumerate(self.kpoints):
                    weight = kpoint.weight
                    for band in range(self.nband_sk[spin,k]):
                        e = self.eigens[spin,k,band]
                        dos[spin] += weight * gaussian(mesh, width, center=e)

        else:
            raise ValueError("Method %s is not supported" % method)

        return ElectronDOS(mesh, dos)

    def get_ejdos(self, spin, valence, conduction,
                  method="gaussian", step=0.1, width=0.2, mesh=None):
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
            for (k, kpoint) in enumerate(self.kpoints):
                weight = kpoint.weight
                for c in conduction:
                    ec = self.eigens[spin,k,c]
                    fc = 1 - self.occfacts[spin,k,c] / full
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

        # Calculate Quasi-particle energies with the scissors operator.
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

        # Change the energies (NB: occupations and fermie are left unchanged).
        return ElectronBands(
            self.structure, self.kpoints, qp_energies, self.fermie, self.occfacts, self.nelect,
            nband_sk=self.nband_sk, smearing=self.smearing, markers=self.markers)

    @add_fig_kwargs
    def plot(self, ax=None, klabels=None, band_range=None, marker=None, width=None, **kwargs):
        """
        Plot the band structure.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            klabels: dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels. e.g. `klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0):"L"}`.
            band_range: Tuple specifying the minimum and maximum band to plot (default: all bands are plotted)
            marker: String defining the marker to plot. Accepts the syntax `markername:fact` where
                fact is a float used to scale the marker size.
            width: String defining the width to plot. Accepts the syntax `widthname:fact` where
                fact is a float used to scale the stripe size.

        Returns:
            `matplotlib` figure
        """
        # Select the band range.
        if band_range is None:
            band_range = range(self.mband)
        else:
            band_range = range(band_range[0], band_range[1], 1)

        ax, fig, plt = get_ax_fig_plt(ax)

        # Decorate the axis (e.g add ticks and labels).
        self.decorate_ax(ax, klabels=klabels) #, title=title)

        # Plot the band energies.
        for spin in self.spins:
            if spin == 0:
                opts = {"color": "black", "linewidth": 2.0}
            else:
                opts = {"color": "red", "linewidth": 2.0}

            for band in band_range:
                self.plot_ax(ax, spin=spin, band=band, **opts)

        # Add markers to the plot.
        if marker is not None:
            try:
                key, fact = marker.split(":")
            except ValueError:
                key = marker
                fact = 1
            fact = float(fact)

            self.plot_marker_ax(ax, key, fact=fact)

        # Plot fatbands.
        if width is not None:
            try:
                key, fact = width.split(":")
            except ValueError:
                key = width
                fact = 1

            self.plot_width_ax(ax, key, fact=fact)

        return fig

    @add_fig_kwargs
    def plot_fatbands(self, klabels=None, **kwargs):  #colormap="jet", max_stripe_width_mev=3.0, qlabels=None, **kwargs):
        """
        Plot the electronic fatbands.

        Args:
            klabels: Dictionary whose keys are tuple with the reduced coordinates of the k-points.
                     The values are the labels. e.g. ~klabels = { (0.0,0.0,0.0):"$\Gamma$", (0.5,0,0):"L" }`.
            band_range: Tuple specifying the minimum and maximum band to plot (default: all bands are plotted)
            width: String defining the width to plot. accepts the syntax `widthname:fact` where
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

        for (ax, key) in zip(ax_list, self.widths):
            # Decorate the axis
            self.decorate_ax(ax, klabels=klabels, title=key)

            # Plot the energies.
            self.plot_ax(ax)

            # Add width around each band.
            self.plot_width_ax(ax, key)

        return fig

    def decorate_ax(self, ax, **kwargs):
        title = kwargs.pop("title", None)
        if title is not None:
            ax.set_title(title)

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

    def plot_ax(self, ax, spin=None, band=None, **kwargs):
        """Helper function to plot the energies for (spin,band) on the axis ax."""
        spin_range = range(self.nsppol) if spin is None else [spin]
        band_range = range(self.mband) if band is None else [band]

        # Disable labels.
        if "label" not in kwargs:
            kwargs["label"] = "_no_legend_" # Actively suppress.

        xx, lines = range(self.nkpt), []
        for spin in spin_range:
            for band in band_range:
                yy = self.eigens[spin,:,band]
                lines.extend(ax.plot(xx, yy, **kwargs))

        return lines

    def plot_width_ax(self, ax, key, spin=None, band=None, fact=1.0, **kwargs):
        """Helper function to plot fatbands for the given (spin,band) on the axis ax."""
        spin_range = range(self.nsppol) if spin is None else [spin]
        band_range = range(self.mband) if band is None else [band]

        facecolor = kwargs.pop("facecolor", "blue")
        alpha = kwargs.pop("alpha", 0.7)

        x, width = range(self.nkpt), fact * self.widths[key]

        for spin in spin_range:
            for band in band_range:
                y, w = self.eigens[spin,:,band], width[spin,:,band] * fact
                ax.fill_between(x, y-w/2, y+w/2, facecolor=facecolor, alpha=alpha)

    def plot_marker_ax(self, ax, key, fact=1.0):
        """Helper function to plot the markers on the axis ax."""
        pos, neg = self.markers[key].posneg_marker()

        # Use different symbols depending on the value of s.
        # Cannot use negative s.
        if pos:
            ax.scatter(pos.x, pos.y, s=np.abs(pos.s)*fact, marker="^", label=key + " >0")

        if neg:
            ax.scatter(neg.x, neg.y, s=np.abs(neg.s)*fact, marker="v", label=key + " <0")

    def _make_ticks_and_labels(self, klabels):
        """Return ticks and labels from the mapping qlabels."""
        if klabels is not None:
            d = OrderedDict()
            for (kcoord, kname) in klabels.items():
                # Build Kpoint instance.
                ktick = Kpoint(kcoord, self.reciprocal_lattice)
                for (idx, kpt) in enumerate(self.kpoints):
                    if ktick == kpt: 
                        d[idx] = kname

        else:
            d = self._auto_klabels

        # Return ticks, labels
        return list(d.keys()), list(d.values())

    @add_fig_kwargs
    def plot_with_edos(self, dos, klabels=None, **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            dos: An instance of :class:`ElectronDOS`.
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. `klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}`.

        Returns:
            `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        gspec = GridSpec(1, 2, width_ratios=[2, 1])
        ax1 = plt.subplot(gspec[0])
        # Align bands and DOS.
        ax2 = plt.subplot(gspec[1], sharey=ax1)

        if not kwargs:
            kwargs = {"color": "black", "linewidth": 2.0}

        # Plot the band structure
        for spin in self.spins:
            for band in range(self.mband):
                self.plot_ax(ax1, spin=spin, band=band, **kwargs)

        self.decorate_ax(ax1, klabels=klabels)

        emin = np.min(self.eigens)
        emin -= 0.05 * abs(emin)

        emax = np.max(self.eigens)
        emax += 0.05 * abs(emax)
        ax1.yaxis.set_view_interval(emin, emax)

        # Plot the DOS
        dos.plot_ax(ax2, exchange_xy=True, **kwargs)

        ax2.grid(True)
        ax2.yaxis.set_ticks_position("right")
        ax2.yaxis.set_label_position("right")

        fig = plt.gcf()
        return fig

    def widget_plot(self):
        from IPython.html import widgets # Widget definitions
        from IPython.display import display # Used to display widgets in the notebook
 
        widget = widgets.FloatSliderWidget()

        #def on_value_change(name, value):
        def on_value_change():
            #print(value)
            #print(self)
            from IPython.display import clear_output
            clear_output()
            self.get_edos().plot() #method=method, step=step, width=width)

        widget.on_trait_change(on_value_change)
        return widget

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
        """Compute the derivative of the eigenvalues along the path."""
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
        ders2 = self.derivatives(spin, band, acc=acc) * units.eV_to_Ha / units.bohr_to_ang**2
        return 1.0/ders2


class ElectronBandsPlotter(object):
    """
    Class for plotting electronic band structure and DOSes.
    Supports plots on the same graph or separated plots.

    Usage example:

    .. code-block:: python
        
        plotter = ElectronBandsPlotter()
        plotter.add_ebands_from_file("foo.nc", label="foo bands")
        plotter.add_ebands_from_file("bar.nc", label="bar bands")
        plotter.plot()
    """
    _LINE_COLORS = ["b", "r",]
    _LINE_STYLES = ["-",":","--","-.",]
    _LINE_WIDTHS = [2,]

    def __init__(self):
        self._bands_dict = OrderedDict()
        self._edoses_dict = OrderedDict()
        self._markers = OrderedDict()

    @property
    def bands_dict(self):
        """Dictionary with the mapping label --> ebands."""
        return self._bands_dict

    @property
    def edoses_dict(self):
        """Dictionary with the mapping label --> edos."""
        return self._edoses_dict

    @property
    def ebands_list(self):
        """"List of `:class:ElectronBands`."""
        return list(self._bands_dict.values())

    @property
    def edoses_list(self):
        """"List of :class:`ElectronDos`."""
        return list(self._edoses_dict.values())

    @property
    def markers(self):
        return self._markers

    def iter_lineopt(self):
        """Generates style options for lines."""
        for o in itertools.product( self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def add_ebands_from_file(self, filepath, label=None):
        """
        Adds a band structure for plotting. Reads data from a Netcdfile
        """
        from abipy.abilab import abiopen
        with abiopen(filepath) as ncfile:
            if label is None:
                label = ncfile.filepath
            self.add_ebands(label, ncfile.ebands)

    def add_ebands(self, label, bands, dos=None):
        """
        Adds a band structure for plotting.

        Args:
            label: label for the bands. Must be unique.
            bands: :class:`ElectronBands` object.
            dos: :class:`ElectronDos` object.
        """
        if label in self._bands_dict:
            raise ValueError("label %s is already in %s" % (label, list(self._bands_dict.keys())))

        self._bands_dict[label] = bands

        if dos is not None:
            self.edoses_dict[label] = dos

    def add_ebands_list(self, labels, bands_list, dos_list=None):
        """
        Add a list of Bands and DOSes.

        Args:
            labels: List of labels.
            bands_list: List of :class:`ElectronBands` objects.
            dos_list: List of :class:`ElectronDos` objects.
        """
        assert len(labels) == len(bands_list)

        if dos_list is None:
            for label, bands in zip(labels, bands_list):
                self.add_ebands(label, bands)
        else:
            assert len(dos_list) == len(bands_list)
            for label, bands, dos in zip(labels, bands_list, dos_list):
                self.add_ebands(label, bands, dos=dos)

    def bands_statdiff(self, ref=0):
        """
        Compare the reference bands with index ref with the other bands stored in the plotter.
        """
        for i, label in enumerate(self._bands_dict.keys()):
            if i == ref:
                ref_label = label
                break
        else:
            raise ValueError("ref index %s is > number of bands" % ref)

        ref_bands = self._bands_dict[ref_label]

        text = []
        for (label, bands) in self._bands_dict.items():
            if label == ref_label: continue
            stat = ref_bands.statdiff(bands)
            text.append(str(stat))

        return "\n\n".join(text)

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

        if extend:
            if key not in self._markers:
                self._markers[key] = Marker(*xys)
            else:
                # Add xys to the previous marker set.
                self._markers[key].extend(*xys)
        
        else:
            if key in self._markers:
                raise ValueError("Cannot overwrite key %s in data" % key)

            self._markers[key] = Marker(*xys)

    @add_fig_kwargs
    def plot(self, klabels=None, **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            klabels: dictionary whose keys are tuple with the reduced
                     coordinates of the k-points. The values are the labels.
                     e.g. klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        xlim            x-axis limits. None (default) for automatic determination.
        ylim            y-axis limits. None (default) for automatic determination.
        ==============  ==============================================================

        Returns:
            matplotlib figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        # Build grid of plots.
        if self.edoses_dict:
            gspec = GridSpec(1, 2, width_ratios=[2, 1])
            ax1 = plt.subplot(gspec[0])
            # Align bands and DOS.
            ax2 = plt.subplot(gspec[1], sharey=ax1)
            ax_list = [ax1, ax2]
            fig = plt.gcf()
        else:
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax_list = [ax1]

        for ax in ax_list:
            ax.grid(True)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None:
            [ax.set_ylim(ylim) for ax in ax_list]

        # Plot bands.
        lines, legends = [], []
        my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        for (label, bands), lineopt in zip(self._bands_dict.items(), self.iter_lineopt()):
            i += 1
            my_kwargs.update(lineopt)
            opts_label[label] = my_kwargs.copy()

            l = bands.plot_ax(ax1, spin=None, band=None, **my_kwargs)
            lines.append(l[0])

            # Use relative paths if label is a file.
            if os.path.isfile(label):
                legends.append("%s" % os.path.relpath(label))
            else:
                legends.append("%s" % label)

            # Set ticks and labels, legends.
            if i == 0:
                bands.decorate_ax(ax1)

        if self.markers:
            for key, markers in self.markers.items():
                pos, neg = markers.posneg_marker()
                # Use different symbols depending on the value of s.
                # Cannot use negative s.
                fact = 1
                if pos:
                    ax1.scatter(pos.x, pos.y, s=np.abs(pos.s)*fact, marker="^", label=key + " >0")

                if neg:
                    ax1.scatter(neg.x, neg.y, s=np.abs(neg.s)*fact, marker="v", label=key + " <0")

        ax1.legend(lines, legends, loc='best', shadow=True)

        # Add DOSes
        if self.edoses_dict:
            ax = ax_list[1]
            for (label, dos) in self.edoses_dict.items():
                dos.plot_ax(ax, exchange_xy=True, **opts_label[label])

        return fig

    def animate_files(self, **kwargs):
        """
        See http://visvis.googlecode.com/hg/vvmovie/images2gif.py for a (much better) approach
        """
        animator = FilesAnimator()
        figures = OrderedDict()

        for label, bands in self.bands_dict.items():
            if self.edoses_dict:
                fig = bands.plot_with_edos(self.edoses_dict[label], show=False)
            else:
                fig = bands.plot(show=False)

            figures[label] = fig

        animator.add_figures(labels=figures.keys(), figure_list=figures.values())
        return animator.animate(**kwargs)

    def animate(self, **kwargs):
        """
        See http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
        """
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation

        fig, ax = plt.subplots()
        bands = list(self.bands_dict.values())

        plot_opts = {"color": "black", "linewidth": 2.0}

        def cbk_animate(i):
            #line.set_ydata(np.sin(x+i/10.0))  # update the data
            #print("in animate with %d" % i)
            return bands[i].plot_ax(ax, spin=None, band=None, **plot_opts)
            #lines = bands[i].plot_ax(ax, spin=None, band=None)
            #line = lines[0]
            #return line

        # initialization function: plot the background of each frame
        def init():
            return bands[0].plot_ax(ax, spin=None, band=None, **plot_opts)
            #line.set_data([], [])
            #return line,

        anim = animation.FuncAnimation(fig, cbk_animate, frames=len(bands), interval=250, blit=True, init_func=init)

        #anim.save('im.mp4', metadata={'artist':'gmatteo'})

        if kwargs.get("show", True): plt.show()

        return anim


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
        from abipy.abilab import abiopen
        with abiopen(filepath) as ncfile:
            ebands = ncfile.ebands

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
    def plot(self, ax=None, **kwargs):
        """
        Plot the band structure and the DOS.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        for (label, dos) in self.edoses_dict.items():
            # Use relative paths if label is a file.
            if os.path.isfile(label): label = os.path.relpath(label)
            dos.plot_ax(ax, label=label)

        ax.grid(True)
        ax.set_xlabel("Energy [eV]")
        ax.set_ylabel("DOS")
        ax.legend(loc="best")

        return fig

    #def animate(self, **kwargs):
    #    animator = Animator()
    #    tmpdir = tempfile.mkdtemp()
    #
    #    for (label, dos) in self.edoses_dict.items():
    #        savefig = os.path.join(tmpdir, label + ".png")
    #        dos.plot(show=False, savefig=savefig)
    #        animator.add_figure(label, savefig)
    #    return animator.animate(**kwargs)


class ElectronsReader(ETSF_Reader, KpointsReaderMixin):
    """This object reads band structure data from a netcdf file written"""
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

    #def read_xcinfo(self):
    #   """Returns a dictionary with info on the XC functional."""
    #    return XcInfo.from_ixc(self.read_value("ixc"))


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
                    Point are located in the Brilluoin zone, i.e kx in [-1/2, 1/2].
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


class ElectronDOS(object):
    """This object stores the electronic density of states."""

    def __init__(self, mesh, spin_dos):
        """
        Args:
            mesh: array-like object with the mesh points.
            spin_dos: array-like object with the DOS value for the different spins.
                      spin_dos[nw] if spin-unpolarized.
                      spin_dos[nsppol, nw] if spin-polarized case.

        .. note::
            mesh is given in eV, spin_dos is in states/eV.
        """
        spin_dos = np.atleast_2d(spin_dos)
        self.nsppol = len(spin_dos)

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

    def dos_idos(self, spin=None):
        """
        Returns DOS and IDOS for given spin. Total DOS and IDOS if spin is None.
        """
        if spin is None:
            return self.tot_dos, self.tot_idos
        else:
            return self.spin_dos[spin], self.spin_idos[spin]

    def find_mu(self, nelect, spin=None, num=500, atol=1.e-5):
        """Finds the chemical potential given the number of electrons."""
        idos = self.tot_idos if spin is None else self.spin_idos[spin]

        # Cannot use bisection because DOS might be negative due to smearing.
        # This one is safer albeit slower.
        for i, (ene, intg) in enumerate(idos):
            if intg > nelect:
                break
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

    def plot_ax(self, ax, spin=None, what="d", exchange_xy=False, *args, **kwargs):
        """
        Helper function to plot the data on the axis ax.

        Args:
            ax: matplotlib axis.
            spin: selects the spin component, None for total DOS, IDOS.
            what: string selecting what will be plotted:
                  "d" for DOS, "i" for IDOS. chars can be concatenated
                  hence what="id" plots both IDOS and DOS. (default "d").
            exchange_xy: True to exchange axis.
            kwargs:
                Options passes to matplotlib.

        Return:
            list of lines that were added.
        """
        #print("spin",spin)
        dosf, idosf = self.dos_idos(spin=spin)

        opts = [c.lower() for c in what]

        lines = []
        for c in opts:
            if c == "d": f = dosf
            if c == "i": f = idosf
            ls = f.plot_ax(ax, exchange_xy=exchange_xy, *args, **kwargs)
            lines.extend(ls)
        return lines

    @add_fig_kwargs
    def plot(self, spin=None, **kwargs):
        """
        Plot DOS and IDOS.

        Args:
            spin: Selects the spin component, None if total DOS is wanted.

        Returns:
            matplotlib figure.
        """
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        gspec = GridSpec(2, 1, height_ratios=[1,2])
        ax1 = plt.subplot(gspec[0])
        ax2 = plt.subplot(gspec[1])

        for ax in (ax1, ax2):
            ax.grid(True)

        ax2.set_xlabel('Energy [eV]')

        ax1.set_ylabel("TOT IDOS" if spin is None else "IDOS (spin %s)" % spin)
        ax2.set_ylabel("TOT DOS" if spin is None else "DOS (spin %s)" % spin)

        self.plot_ax(ax1, spin=spin, what="i", **kwargs)
        self.plot_ax(ax2, spin=spin, what="d", **kwargs)

        fig = plt.gcf()
        return fig


class ElectronDOSPlotter(object):
    """
    Class for plotting multiple electron DOSes.
    """
    def __init__(self):
        self._doses = OrderedDict()

    def add_dos(self, label, dos):
        """
        Adds a DOS for plotting.

        Args:
            label: label for the DOS. Must be unique.
            dos: :class:`ElectroDos` object.
        """
        if label in self._doses:
            raise ValueError("label %s is already in %s" % (label, self._doses.keys()))

        self._doses[label] = dos

    def add_dos_dict(self, dos_dict, key_sort_func=None):
        """
        Add a dictionary of DOSes, with an optional sorting function for the keys.

        Args:
            dos_dict: dict of {label: dos}
            key_sort_func: function used to sort the dos_dict keys.
        """
        if key_sort_func:
            keys = sorted(dos_dict.keys(), key=key_sort_func)
        else:
            keys = dos_dict.keys()

        for label in keys:
            self.add_dos(label, dos_dict[label])

    @add_fig_kwargs
    def plot(self, ax=None, *args, **kwargs):
        """
        Get a matplotlib plot showing the DOSes.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        xlim            x-axis limits. None (default) for automatic determination.
        ylim            y-axis limits.  None (default) for automatic determination.
        ==============  ==============================================================
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        ax.grid(True)

        xlim = kwargs.pop("xlim", None)
        if xlim is not None: ax.set_xlim(xlim)

        ylim = kwargs.pop("ylim", None)
        if ylim is not None: ax.set_ylim(ylim)

        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel('DOS [states/eV]')

        lines, legends = [], []
        for (label, dos) in self._doses.items():
            l = dos.plot_ax(ax, *args, **kwargs)[0]

            lines.append(l)
            legends.append("DOS: %s" % label)

        # Set legends.
        ax.legend(lines, legends, loc='best', shadow=True)

        return fig
