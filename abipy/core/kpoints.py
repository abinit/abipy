# coding: utf-8
"""This module defines objects describing the sampling of the Brillouin Zone."""
from __future__ import print_function, division, unicode_literals

import collections
import numpy as np

from monty.collections import AttrDict
from monty.functools import lazy_property
from abipy.iotools import ETSF_Reader
from abipy.tools.derivatives import finite_diff

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "issamek",
    "wrap_to_ws",
    "wrap_to_bz",
    "as_kpoints",
    "Kpoint",
    "Kpath",
    "IrredZone",
    "rc_list",
    "kmesh_from_mpdivs",
]

# Tolerance used to compare k-points.
_ATOL_KDIFF = 1e-8


def is_integer(x, atol=_ATOL_KDIFF):
    """
    True if all x is integer within the absolute tolerance atol.

    >>> assert is_integer([1., 2.])
    >>> assert is_integer(1.01, atol=0.011)
    >>> assert not is_integer([1.01, 2])
    """
    int_x = np.around(x)
    return np.allclose(int_x, x, atol=atol)

    #return (np.abs(int_x[0] - x[0]) < atol and 
    #        np.abs(int_x[1] - x[1]) < atol and 
    #        np.abs(int_x[2] - x[2]) < atol )


def issamek(k1, k2, atol=1e-08):
    """
    True if k1 and k2 are equal modulo a lattice vector.

    >>> assert issamek([1,1,1], [0,0,0])
    >>> assert issamek([1.1,1,1], [0,0,0], atol=0.1)
    >>> assert not issamek(0.00003, 1)
    """
    return is_integer(np.asarray(k1)-np.asarray(k2), atol=atol)


def wrap_to_ws(x):
    """
    Transforms x in its corresponding reduced number in the interval ]-1/2,1/2].
    """
    w = x % 1
    return np.where(w > 0.5, w-1.0, w)


def wrap_to_bz(x):
    """
    Transforms x in its corresponding reduced number in the interval [0,1[."
    """
    return x % 1


def rc_list(mp, sh, pbc=False, order="bz"):
    """
    Returns a `ndarray` with the linear mesh used to sample one dimension of the reciprocal space.
    Note that the coordinates are always ordered so that rc[i+1] > rc[i].
    so that we can easily plot quantities defined on the 3D multidimensional mesh.

    Args:
        mp:
            Number of Monkhorst-Pack divisions along this direction.
        sh:
            Shift 
        pbc:
            if pbc is True, periodic boundary conditions are enforced.
        order:
            Possible values ["bz", "unit_cell"]. 
            if "bz", the coordinates are forced to be in [-1/2, 1/2)
            if "unit_cell", the coordinates are forced to be in [0, 1).
    """
    rc = []
                                                   
    if order == "unit_cell":
        n = mp if not pbc else mp + 1
        for i in range(n):
            rc.append((i + sh) / mp)
                                                   
    elif order == "bz":
        for i in range(mp):
            x = (i + sh) / mp

            if x < 0.5:
                rc.append(x)
            else:
                # Insert xm1 in rc so that we still have a ordered list.
                xm1 = x - 1.0
                for i, c in enumerate(rc):
                    if c > xm1:
                        break
                else:
                    raise ValueError()

                rc.insert(i, xm1)
                                                   
        if pbc:
            rc.append(rc[0] + 1.0)
                                                   
    else:
        raise ValueError("Wrong order %s" % order)
                                                   
    return np.array(rc)


def kmesh_from_mpdivs(mpdivs, shifts, pbc=False, order="bz"):
    """
    Returns a `ndarray` with the reduced coordinates of the 
    k-points from the MP divisions and the shifts.

    Args:
        mpdivs: The three MP divisions
        shifts: Array-like object with the MP shift.
        pbc: If True, periodic images of the k-points will be includes i.e. closed mesh.
        order: "unit_cell" if the kpoint coordinates must be in [0,1)
               "bz" if the kpoint coordinates must be in [-1/2, +1/2)
    """
    shifts = np.reshape(shifts, (-1,3))
    assert np.all(np.abs(shifts) <= 0.5)

    # Build k-point grid.
    from itertools import product
    kbz = []
    for ish, shift in enumerate(shifts):
        rc0 = rc_list(mpdivs[0], shift[0], pbc=pbc, order=order)
        rc1 = rc_list(mpdivs[1], shift[1], pbc=pbc, order=order)
        rc2 = rc_list(mpdivs[2], shift[2], pbc=pbc, order=order)

        for kxyz in product(rc0, rc1, rc2):
            kbz.append(kxyz)

    return np.array(kbz)


class KpointsError(Exception):
    """Base error class for KpointList exceptions."""


def as_kpoints(obj, lattice, weights=None, names=None):
    """
    Convert obj into a list of k-points.

    Args:
        obj: :class:`Kpoint` or list of Kpoint objects or array-like object.
        lattice: Reciprocal lattice.
        weights: k-point weights. Ignored if obj is already a `Kpoint` instance or a list
                 of `Kpoint` items.
        name: string with the name of the k-point. Ignored if obj is already a `Kpoint`
              instance or a list of `Kpoint` items.
    """
    # K-point?
    if isinstance(obj, Kpoint):
        return [obj]

    # Iterable with K-points?
    if isinstance(obj, collections.Iterable):
        if isinstance(obj[0], Kpoint):
            assert all([isinstance(o, Kpoint) for o in obj])
            return obj

    # Assume array-like
    obj = np.reshape(np.asarray(obj), (-1, 3))
    ndim = obj.ndim

    if ndim == 1:
        return [Kpoint(obj, lattice, weight=weights, name=names)]

    elif ndim == 2:
        nk = len(obj)
        if weights is None: weights = nk * [None]
        if names is None: names = nk * [None]
        return [Kpoint(rc, lattice, weight=w, name=l) for (rc, w, l) in zip(obj, weights, names)]

    else:
        raise ValueError("ndim > 2 is not supported")


class Kpoint(object):
    """Class defining one k-point."""

    # TODO: Fix problem with pickle
    #__slots__ = [
    #    "_frac_coords",
    #    "_lattice",
    #    "_weight",
    #    "_name",
    #    "_hash",
    #]

    # Tolerance used to compare k-points.
    @property
    def ATOL_KDIFF(self):
        return _ATOL_KDIFF

    def __init__(self, frac_coords, lattice, weight=None, name=None):
        """
        Args:
            frac_coords: Reduced coordinates.
            lattice: :class:`Lattice` object describing the reciprocal lattice.
            weights: k-point weight (optional, set to zero if not given).
            name: string with the name of the k-point (optional)
        """
        self._frac_coords = np.asarray(frac_coords)
        assert len(self.frac_coords) == 3

        self._lattice = lattice
        self.set_weight(weight)
        self.set_name(name)

    #def __array__(self, **kwargs):
    #    """np.array(self)"""
    #    print(kwargs)
    #    dtype = kwargs.pop("dtype", None)
    #    if dtype is None:
    #        return self._frac_coords
    #    else:
    #        return np.array(self._frac_coords, dtype=dtype)

    def __hash__(self):
        """
        Kpoint objects can be used as keys in dictionaries.
        
        .. warning::

            The hash is computed from the fractional coordinates (floats). 
            Hence one should avoid using hashes for implementing search algorithms
            in which new Kpoints are, for example generated by means of 
            symmetry operations. This means that a dict of Kpoint objects
            is safe to use only when we are sure than we are going to access
            its entries with the *same* keys used to generate the dict!.
        """
        try:
            return self._hash
        except AttributeError:
            self._hash = hash(tuple(wrap_to_ws(self.frac_coords)))
            return self._hash

    @property
    def frac_coords(self):
        """Fractional coordinates of the k-points."""
        return self._frac_coords

    @property
    def lattice(self):
        """Reciprocal lattice."""
        return self._lattice

    @property
    def weight(self):
        """Weight of the k-point. 0.0 if the weight is not defined."""
        if self._weight is None:
            return 0.0
        else:
            return self._weight

    def set_weight(self, weight):
        """Set the weight of the k-point."""
        self._weight = weight

    @property
    def cart_coords(self):
        """Cartesian coordinates of the k-point."""
        return self.lattice.get_cartesian_coords(self.frac_coords)

    @property
    def name(self):
        """Name of the k-point. None if not available."""
        return self._name

    def set_name(self, name):
        """Set the name of the k-point."""
        # Fix typo in Latex syntax (if any).
        if name is not None and name.startswith("\\"): name = "$" + name + "$"
        self._name = name

    @property
    def on_border(self):
        """
        True if the k-point is on the border of the BZ  (lattice translations are taken into account).
        """
        kreds = wrap_to_ws(self.frac_coords)
        diff = np.abs(np.abs(kreds) - 0.5)
        return np.any(diff < self.ATOL_KDIFF)

    def __repr__(self):
        return "[%.3f, %.3f, %.3f]" % tuple(self.frac_coords)

    def __str__(self):
        s = repr(self)
        if self.name is not None: s += ", name = %s" % self.name
        if self._weight is not None: s += ", weight = %f" % self.weight

        return s

    # Kpoint algebra.
    def __add__(self, other):
        return self.__class__(self.frac_coords + other.frac_coords, self.lattice)

    def __sub__(self, other):
        return self.__class__(self.frac_coords - other.frac_coords, self.lattice)

    def __eq__(self, other):
        try:
            # Comparison between two Kpoint objects
            return issamek(self.frac_coords, other.frac_coords, atol=self.ATOL_KDIFF)

        except AttributeError:
            # Kpoint vs iterable (e.g. list)
            return issamek(self.frac_coords, other, atol=self.ATOL_KDIFF)

    def __ne__(self, other):
        return not self == other

    def __getitem__(self, slice):
        return self.frac_coords[slice]

    @classmethod
    def as_kpoint(cls, obj, lattice):
        """
        Convert obj into a Kpoint instance.

        Args:
            obj:
                :class:`Kpoint` instance or array-like with the reduced coordinates.
            lattice:
                :class:`Lattice` object defining the reciprocal lattice.
        """
        if isinstance(obj, cls):
            return obj
        else:
            return cls(obj, lattice, weight=None, name=None)

    @classmethod
    def gamma(cls, lattice, weight=None):
        """Constructor for the Gamma point."""
        return cls(np.zeros(3), lattice, weight=weight, name="$\Gamma$")

    def copy(self):
        """Deep copy."""
        return self.__class__(self.frac_coords.copy(), self.lattice.copy(),
                              weight=self.weight, name=self.name)

    @property
    def norm(self):
        """Norm of the kpoint."""
        return np.sqrt(np.dot(self.cart_coords, self.cart_coords))

    def versor(self):
        """Returns the versor i.e. ||k|| = 1"""
        cls = self.__class__
        try:
            return cls(self.frac_coords / self.norm, self.lattice, weight=self.weight)
        except ZeroDivisionError:
            return cls.gamma(self.lattice, weight=self.weight)

    def wrap_to_ws(self):
        """Returns a new `Kpoint` in the Wigner-Seitz zone."""
        return self.__class__(wrap_to_ws(self.frac_coords), self.lattice, 
                              name=self.name, weight=self.weight)

    def wrapt_to_bz(self):
        """Returns a new `Kpoint` in the first unit cell."""
        return self.__class__(wrap_to_bz(self.frac_coords), self.lattice, 
                              name=self.name, weight=self.weight)
        
    def compute_star(self, symmops, wrap_tows=True):
        """Return the star of the kpoint (tuple of `Kpoint` objects)."""
        frac_coords = [self.frac_coords]

        for sym in symmops:
            sk_coords = sym.rotate_k(self.frac_coords, wrap_tows=wrap_tows)

            # Add it only if it's not already in the list.
            for prev_coords in frac_coords:
                if issamek(sk_coords, prev_coords, atol=self.ATOL_KDIFF):
                    break
            else:
                frac_coords.append(sk_coords)

        return KpointStar(self.lattice, frac_coords, weights=None, names=len(frac_coords) * [self.name])


class KpointList(collections.Sequence):
    """
    Base class defining a sequence of :class:`Kpoint` objects. Essentially consists
    of base methods implementing the sequence protocol and helper functions.
    """
    Error = KpointsError

    def __init__(self, reciprocal_lattice, frac_coords, weights=None, names=None):
        """
        Args:
            reciprocal_lattice: :class:`Lattice` object.
            frac_coords: Array-like object with the reduced coordinates of the k-points.
            weights: List of k-point weights.
            names: List of k-point names.
        """
        self._reciprocal_lattice = reciprocal_lattice

        self._frac_coords = frac_coords = np.reshape(frac_coords, (-1, 3))

        if weights is not None:
            assert len(weights) == len(frac_coords)
            weights = np.asarray(weights)
        else:
            weights = np.zeros(len(self.frac_coords))

        if names is not None:
            assert len(names) == len(frac_coords)

        self._points = []
        for i, rcs in enumerate(frac_coords):
            name = None if names is None else names[i]
            self._points.append(Kpoint(rcs, self.reciprocal_lattice, weight=weights[i], name=name))

    @classmethod
    def from_file(cls, filepath):
        if filepath.endswith(".nc"):
            with KpointsReader(filepath) as r:
                new = r.read_kpoints()
        else:
            raise NotImplementedError("")

        new.__class__ == cls
        return new

    @property
    def reciprocal_lattice(self):
        """`Lattice` object defining the reciprocal lattice."""
        return self._reciprocal_lattice

    def __repr__(self):
        lines = ["%d) %s" % (i, repr(kpoint)) for i, kpoint in enumerate(self)]
        return "\n".join(lines)

    def __str__(self):
        lines = ["%d) %s" % (i, str(kpoint)) for i, kpoint in enumerate(self)]
        return "\n".join(lines)

    # Sequence protocol.
    def __len__(self):
        return len(self._points)

    def __iter__(self):
        return self._points.__iter__()

    def __getitem__(self, slice):
        return self._points[slice]

    def __contains__(self, kpoint):
        return kpoint in self._points

    def __reversed__(self):
        return self._points.__reversed__()

    def __add__(self, other):
        assert self.reciprocal_lattice == other.reciprocal_lattice
        return KpointList(self.reciprocal_lattice, 
                          frac_coords=[k.frac_coords for k in self] + [k.frac_coords for k in other], 
                          weights=None,
                          names=[k.name for k in self] + [k.name for k in other],
                        )

    def __eq__(self, other):
        for k1, k2 in zip(self, other):
            if k1 != k2:
                return False

        return True
        
    def __ne__(self, other):
        return not self == other

    def index(self, kpoint):
        """
        Returns: the first index of kpoint in self.

        Raises: ValueError if not found.
        """
        try:
            return self._points.index(kpoint)

        except ValueError:
            raise ValueError("\nCannot find point: %s in KpointList:\n%s" % (repr(kpoint), repr(self)))

    def find(self, kpoint):
        """
        Returns: first index of kpoint. -1 if not found
        """
        try:
            return self.index(kpoint)
        except ValueError:
            return -1

    def count(self, kpoint):
        """Return number of occurrences of kpoint"""
        return self._points.count(kpoint)

    @property
    def is_path(self):
        """True if self represents a path in the BZ."""
        return isinstance(self, Kpath)

    @property
    def is_homogeneous(self):
        """True if self represents a homogeneous sampling of the BZ."""
        return isinstance(self, IrredZone)

    @property
    def frac_coords(self):
        """Fractional coordinates of the k-point as `ndarray` of shape (len(self), 3)"""
        return self._frac_coords

    @property
    def weights(self):
        """`ndarray` with the weights of the k-points."""
        return np.array([kpoint.weight for kpoint in self])

    def sum_weights(self):
        """Returns the sum of the weights."""
        return np.sum(self.weights)

    def remove_duplicated(self):
        """Remove duplicated k-points from self. Returns new KpointList instance."""
        frac_coords, good_indices = [self[0].frac_coords], [0]

        for i, kpoint in enumerate(self[1:]):
            i += 1
            # Add it only if it's not already in the list.
            for prev_coords in frac_coords:
                if issamek(kpoint.frac_coords, prev_coords, atol=_ATOL_KDIFF):
                    break
            else:
                frac_coords.append(kpoint.frac_coords)
                good_indices.append(i)

        good_kpoints = [self[i] for i in good_indices]

        return KpointList(self.reciprocal_lattice, 
                          frac_coords=[k.frac_coords for k in good_kpoints],
                          weights=None,
                          names=[k.name for k in good_kpoints])

    def to_array(self):
        """Returns a `ndarray` [nkpy, 3] with the fractional coordinates."""
        return np.array(self.frac_coords.copy())


class KpointStar(KpointList):
    """
    Start of the kpoint. Note that the first k-point is assumed to be the base 
    of the star namely the point that is used to generate the Star.
    """
    @property
    def base_point(self):
        """The point used to generate the star."""
        return self[0]

    @property
    def name(self):
        """The name of the star (inherited from the name of the base point)."""
        return self.base_point.name


class Kpath(KpointList):
    """This object describes a path in reciprocal space."""

    def __init__(self, reciprocal_lattice, frac_coords):
        """
        Args:
            reciprocal_lattice: :class:`Lattice` object.
            frac_coords: Array-like object with the reduced coordinates of the k-points.
        """
        super(Kpath, self).__init__(reciprocal_lattice, frac_coords)

    @lazy_property
    def ds(self):
        """
        ndarray of len(self)-1 elements giving the distance between two
        consecutive k-points, i.e. ds[i] = ||k[i+1] - k[i]||.
        """
        ds = np.zeros(len(self) - 1)
        for (i, kpoint) in enumerate(self[:-1]):
            ds[i] = (self[i + 1] - kpoint).norm
        return ds

    @lazy_property
    def versors(self):
        """Tuple of len(self)-1 elements with the versors connecting k[i] to k[i+1]."""
        versors = (len(self) - 1) * [None, ]
        versors[0] = Kpoint.gamma(self.reciprocal_lattice)

        for (i, kpt) in enumerate(self[:-1]):
            versors[i] = (self[i + 1] - kpt).versor()

        return tuple(versors)

    @property
    def num_lines(self):
        """The number of lines forming the path."""
        return len(self.lines)

    @lazy_property
    def lines(self):
        """
        tuple with the list of indices of the points belonging to the same line.
        """
        lines = []
        prev, indices = self.versors[0], [0]

        for (i, v) in enumerate(self.versors[1:]):
            i += 1
            if v != prev:
                prev = v
                lines.append(indices + [i])
                indices = [i]
            else:
                indices += [i]

            lines.append(indices + [len(self) - 1])

        return tuple(lines)

    def finite_diff(self, values, order=1, acc=4):
        """
        Compute the derivatives of values by finite differences.

        Args:
            values: array-like object with the values of the path.
            order: Order of the derivative.
            acc: Accuracy: 4 corresponds to a central difference with 5 points.

        Returns:
            ndarray with the derivative.
        """
        values = np.asarray(values)
        assert len(values) == len(self)

        # Loop over the lines of the path, extract the data and
        # differenciate f(s) where s is the distance along the line.
        ders_on_lines = []

        for line in self.lines:
            vals_on_line = values[line]
            h = self.ds[line[0]]

            if not np.allclose(h, self.ds[line[:-1]]):
                raise ValueError("For finite difference derivatives, the path must be homogeneous!\n" +
                                 str(self.ds[line[:-1]]))

            der = finite_diff(vals_on_line, h, order=order, acc=acc)
            ders_on_lines.append(der)

        return np.array(ders_on_lines)


class IrredZone(KpointList):
    """
    An :class:`IrredZone` is a (immutable) sequence of points in the irreducible wedge of the BZ.
    Each point has a weight whose sum must equal 1 so that we can integrate quantities in the full Brillouin zone.
    Provides methods to symmetrize k-dependent quantities with the full symmetry of the structure. e.g.
    bands, occupation factors, phonon frequencies.
    """
    def __init__(self, reciprocal_lattice, frac_coords, weights, ksampling):
        """
        Args:
            reciprocal_lattice: :class:`Lattice` object
            frac_coords: Array-like object with the reduced coordinates of the points.
            weights: Array-like with the weights of the k-points.
            ksampling:
                TODO
        """
        super(IrredZone, self).__init__(reciprocal_lattice, frac_coords, weights=weights, names=None)

        # Weights must be normalized to one.
        wsum = self.sum_weights()

        if abs(wsum - 1) > 1.e-6:
            err_msg = "Kpoint weights should sum up to one while sum_weights is %.3f\n" % wsum
            err_msg += "The list of kpoints does not represent a homogeneous sampling of the BZ\n" 
            #err_msg += str(type(self)) + "\n" + str(self)
            #raise ValueError(err_msg)  # GA : Should not prevent a band structure from being read!
            logger.warning(err_msg)

        # FIXME
        # Quick and dirty hack to allow the reading of the k-points from WFK files
        # where info on the sampling is missing. I will regret it but at present 
        # is the only solution I found (changes in the ETSF-IO part of Abinit are needed)
        return 

        # FIXME: Check the treatment of the shifts, kptrlatt ...
        # time-reversal?
        assert ksampling.is_homogeneous
        self.kptopt = ksampling.kptopt

        shifts = ksampling.shifts 
        if shifts is None: shifts = [0.0, 0.0, 0.0]

        self._shifts = np.reshape(shifts, (-1,3))

        if ksampling.kptrlatt is not None:
            # Diagonal kptrlatt is equivalent to MP folding.
            # Non-diagonal kptrlatt is not supported.
            self.mpdivs = np.array(3, dtype=np.int)
            self.mpdivs[0] = ksampling.kptrlatt[0,0]
            self.mpdivs[1] = ksampling.kptrlatt[1,1]
            self.mpdivs[2] = ksampling.kptrlatt[2,2]

            for i in range(3):
                for j in range(3):
                    if ksampling.kptrlatt[i,j] != 0:
                        raise ValueError("Non diagonal kptrlatt is not supported")

        elif ksampling.mpdivs is not None:
            # MP folding
            self.mpdivs = ksampling.mpdivs

        else:
            raise ValueError("Either MP mesh info or kptrlatt must be present in ksampling")

        #self.nx, self.ny, self.nz = self.mpdivs
        #grids_1d = 3 * [None]
        #for i in range(3):
        #    grids_1d[i] = np.arange(0, self.mpdivs[i])
        #self.grids_1d = tuple(grids_1d)

    @property
    def shifts(self):
        """`ndarray` with the shifts."""
        return self._shifts
                                         
    @property
    def num_shifts(self):
        """Number of shifts."""
        return len(self.shifts)

    @property
    def len_bz(self):
        """Number of points in the full BZ."""
        return self.mpdivs.prod() * self.num_shifts

    #@lazy_property
    #def ktab(self):
    #Compute the mapping bz --> ibz
    #from abipy.extensions.klib import map_bz2ibz
    #return map_bz2ibz(structure=structure, bz_arr=self.bz_arr, ib_arr=self.ibz_arr)

    #def iter_bz_coords(self):
    #    """
    #    Generates the fractional coordinates of the points in the BZ.

    #    .. note:

    #        points are ordered in blocks, one block for each shift.
    #        Inside the block, points are ordered following the C convention.
    #    """
    #    for shift in self.shifts:

    #        for i in range(mpdivs[0]):
    #            x = (i + shift[0]) / mpdivs[0]
    #            for j in range(mpdivs[1]):
    #                y = (j + shift[1]) / mpdivs[1]
    #                for k in range(mpdivs[2]):
    #                    z = (k + shift[2]) / mpdivs[2]
    #                    yield np.array((x, y, z))

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


class KSamplingInfo(AttrDict):

    #KNOWN_KEYS = {
    #    "reduced_coordinates_of_kpoints": MANDATORY,
    #    "kpoint_weights": OPTIONAL,
    #    "kpoint_grid shift": OPTIONAL,
    #    "monkhorst_pack_folding": OPTIONAL,
    #    "kpoint_grid vectors": OPTIONAL,
    #    "kptopt": OPTIONAL,
    #}

    def __init__(self, *args, **kwargs):
        super(KSamplingInfo, self).__init__(*args, **kwargs)
        #print("ksampling", self)
        #for k in self:
        #   if k not in self.KNOWN_KEYS:
        #       raise ValueError("Unknow key %s" % k)

    @property
    def is_homogeneous(self):
        """True if we have a homogeneous sampling of the BZ."""
        return self.mpdivs is not None or self.kptrlatt is not None

    @property
    def is_path(self):
        """True if we have a path in the BZ."""
        return not self.is_homogeneous


def returns_None_onfail(func):
    import functools
    from numpy.ma import MaskedArray

    @functools.wraps(func)
    def wrapper(self):
        try:
            value = func(self)
            # This trick is needed because in many cases we define the netcdf variable
            # but we don't write its value. 
            return value if not isinstance(value, MaskedArray) else None
        except self.Error:
            return None

    return wrapper


class KpointsReaderMixin(object):
    """
    Mixin class that provides methods for reading k-point data from a netcdf 
    file written according to the ETSF-IO specifications.
    """
    def read_kpoints(self):
        """
        Factory function: returns an instance of [Kpath, IrredZone]
        depending on the content of the Netcdf file. Main entry point for client code.
        """
        structure = self.read_structure()

        frac_coords = self.read_kfrac_coords()
        weights = self.read_kweights()
        ksampling = self.read_ksampling_info()

        # FIXME
        # Quick and dirty hack to allow the reading of the k-points from WFK files
        # where info on the sampling is missing. I will regret it but at present 
        # is the only solution I found (changes in the ETSF-IO part of Abinit are needed)
        if ksampling.is_homogeneous or abs(sum(weights) - 1.0) < 1.e-6:
            # we have a homogeneous sampling of the BZ.
            return IrredZone(structure.reciprocal_lattice, frac_coords, weights, ksampling)

        elif ksampling.is_path:
            # we have a path in the BZ.
            return Kpath(structure.reciprocal_lattice, frac_coords)

        else:
            raise ValueError("Only homogeneous or path samplings are supported!")

    def read_ksampling_info(self):
        return KSamplingInfo(
            shifts=self.read_kshifts(),
            mpdivs=self.read_kmpdivs(),
            kptrlatt=self.read_kptrlatt(),
            kptopt=self.read_kptopt(),
        )

    def read_kfrac_coords(self):
        """Fractional coordinates of the k-points"""
        return self.read_value("reduced_coordinates_of_kpoints")

    @returns_None_onfail
    def read_kweights(self):
        """Returns the weight of the k-points. None if not found."""
        return self.read_value("kpoint_weights")

    @returns_None_onfail
    def read_kshifts(self):
        """Returns the shifts of the k-mesh in reduced coordinates. None if not found."""
        return self.read_value("kpoint_grid_shift")

    @returns_None_onfail
    def read_kmpdivs(self):
        """Returns the Monkhorst-Pack divisions defining the mesh. None if not found."""
        return self.read_value("monkhorst_pack_folding")

    @returns_None_onfail
    def read_kptrlatt(self):
        """Returns ABINIT variable kptrlatt. None if not found."""
        return self.read_value("kpoint_grid_vectors")

    @returns_None_onfail
    def read_kptopt(self):
        """Returns the ABINIT variable kptopt. None if not found."""
        return self.read_value("kptopt")


class KpointsReader(ETSF_Reader, KpointsReaderMixin):
    """This object reads k-point data from a netcdf file."""
