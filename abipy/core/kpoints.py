# coding: utf-8
"""This module defines objects describing the sampling of the Brillouin Zone."""
from __future__ import print_function, division, unicode_literals, absolute_import

import collections
import json
import sys
import numpy as np

from itertools import product
from tabulate import tabulate
from monty.json import MSONable, MontyEncoder
from monty.collections import AttrDict
from monty.functools import lazy_property
from pymatgen.core.lattice import Lattice
from pymatgen.serializers.json_coders import pmg_serialize
from pymatgen.serializers.pickle_coders import SlotPickleMixin
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
    "Ktables",
]

# Tolerance used to compare k-points.
_ATOL_KDIFF = 1e-8

# Tolerances passed to spglib.
_SPGLIB_SYMPREC = 1e-5
_SPGLIB_ANGLE_TOLERANCE = -1.0

def set_atol_kdiff(new_atol):
    """
    Change the value of the tolerance `_ATOL_KDIFF` used to compare k-points.
    Return old value

    .. warning::

        This function should be called at the beginning of the script.
    """
    global _ATOL_KDIFF
    old_atol = _ATOL_KDIFF
    _ATOL_KDIFF = new_atol
    return old_atol

def set_spglib_tols(symprec, angle_tolerance):
    """
    Change the value of the tolerances `symprec` and `angle_tolerance`
    used to call spglib. Return old values

    .. warning::

        This function should be called at the beginning of the script.
    """
    global _SPGLIB_SYMPREC, _SPGLIB_ANGLE_TOLERANCE
    old_symprec, old_angle_tolerance = _SPGLIB_SYMPREC, _SPGLIB_ANGLE_TOLERANCE
    _SPGLIB_SYMPREC, _SPGLIB_ANGLE_TOLERANCE = symprec, angle_tolerance
    return old_symprec, old_angle_tolerance


def is_integer(x, atol=None):
    """
    True if all x is integer within the absolute tolerance atol.
    Use _ATOL_KDIFF is atol is None.

    >>> assert is_integer([1., 2.])
    >>> assert is_integer(1.01, atol=0.011)
    >>> assert not is_integer([1.01, 2])
    """
    if atol is None: atol = _ATOL_KDIFF
    int_x = np.around(x)
    return np.allclose(int_x, x, atol=atol)
    #return (np.abs(int_x[0] - x[0]) < atol and
    #        np.abs(int_x[1] - x[1]) < atol and
    #        np.abs(int_x[2] - x[2]) < atol )


def issamek(k1, k2, atol=None):
    """
    True if k1 and k2 are equal modulo a lattice vector.
    Use _ATOL_KDIFF is atol is None.

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
        mpdivs: The three MP divisions.
        shifts: Array-like object with the MP shift.
        pbc: If True, periodic images of the k-points will be included i.e. closed mesh.
        order: "unit_cell" if the kpoint coordinates must be in [0,1)
               "bz" if the kpoint coordinates must be in [-1/2, +1/2)
    """
    shifts = np.reshape(shifts, (-1,3))
    assert np.all(np.abs(shifts) <= 0.5)

    # Build k-point grid.
    kbz = collections.deque()
    for ish, shift in enumerate(shifts):
        rc0 = rc_list(mpdivs[0], shift[0], pbc=pbc, order=order)
        rc1 = rc_list(mpdivs[1], shift[1], pbc=pbc, order=order)
        rc2 = rc_list(mpdivs[2], shift[2], pbc=pbc, order=order)

        for kxyz in product(rc0, rc1, rc2):
            kbz.append(kxyz)

    return np.array(kbz)


def map_bz2ibz(structure, ibz, ngkpt, has_timrev, pbc=False):
    """
    Compute the correspondence between the list of k-points in the *unit cell*
    associated to the `ngkpt` mesh and the corresponding points in the IBZ.
    This routine is mainly used to symmetrize eigenvalues in the unit cell
    e.g. to write BXSF files for electronic isosurfaces.

    Args:
        structure: Structure with symmetry operations.
        ibz: [*,3] array with reduced coordinates in the in the IBZ.
        ngkpt: Mesh divisions.
        has_timrev: True if time-reversal can be used.
        pbc: True if the mesh should contain the periodic images (closed mesh).

    Returns:
        bz2ibz: 1d array with BZ --> IBZ mapping
    """
    ngkpt = np.asarray(ngkpt, dtype=np.int)

    # Extract (FM) symmetry operations in reciprocal space.
    abispg = structure.abi_spacegroup
    symrec_fm = [s for (s, afm) in zip(abispg.symrec, abispg.symafm) if afm == 1]

    # Compute TS k_ibz.
    bzgrid2ibz = -np.ones(ngkpt, dtype=np.int)

    for ik_ibz, kibz in enumerate(ibz):
        gp_ibz = np.array(np.rint(kibz * ngkpt), dtype=np.int)
        for rot in symrec_fm:
            rot_gp = np.matmul(rot, gp_ibz)
            gp_bz = rot_gp % ngkpt
            bzgrid2ibz[gp_bz[0], gp_bz[1], gp_bz[2]] = ik_ibz
            if has_timrev:
                gp_bz = (-rot_gp) % ngkpt
                bzgrid2ibz[gp_bz[0], gp_bz[1], gp_bz[2]] = ik_ibz

    from abipy.tools.numtools import add_periodic_replicas
    if pbc:
        bzgrid2ibz = add_periodic_replicas(bzgrid2ibz)
    bz2ibz = bzgrid2ibz.flatten()

    if np.any(bz2ibz == -1):
        #for ik_bz, ik_ibz in enumerate(self.bz2ibz): print(ik_bz, ">>>", ik_ibz)
        msg = "Found %s/%s invalid entries in bz2ibz array" % ((bz2ibz == -1).sum(), len(bz2ibz))
        msg += "This can happen if there an incosistency between the input IBZ and ngkpt"
        msg += "ngkpt: %s, has_timrev: %s" % (str(ngkpt), has_timrev)
        raise ValueError(msg)

    return bz2ibz

    """
    for ik_bz, kbz in enumerate(bz):
        found = False
        for ik_ibz, kibz in enumerate(ibz):
            if found: break
            for symmop in structure.spacegroup:
                krot = symmop.rotate_k(kibz)
                if issamek(krot, kbz):
                    bz2ibz[ik_bz] = ik_ibz
                    found = True
                    break
    """


def has_timrev_from_kptopt(kptopt):
    """
    True if time-reversal symmetry can be used in the generation of the k-points in the IBZ.
    """
    kptopt = int(kptopt)
    return False if kptopt in (3, 4) else True


def map_kpoints(other_kpoints, other_lattice, ref_lattice, ref_kpoints, ref_symrecs, has_timrev):
    """
    Build mapping between a list of k-points in reduced coordinates (`other_kpoints`)
    in the reciprocal lattice `other_lattice` and a list of reference k-points given
    in the reciprocal lattice `ref_lattice` with symmetry operations `ref_symrecs`.

    Args:
        other_kpoints:
        other_lattice: matrix whose rows are the reciprocal lattice vectors in cartesian coordinates.
        ref_lattice: same meaning as other_lattice.
        ref_kpoints:
        ref_symrecs: [nsym,3,3] arrays with symmetry operations in the `ref_lattice` reciprocal space.
        has_timrev: True if time-reversal can be used.

    Returns
        (o2r_map, nmissing)

        nmissing:
            Number of k-points in ref_kpoints that cannot be mapped onto ref_kpoints.

        o2r_map[i] gives the mapping  between the i-th k-point in other_kpoints and
            ref_kpoints. Set to None if the i-th k-point does not have any image in ref.
            Each entry is a named tuple with the following attributes:

                ik_ref:
                tsign:
                isym
                g0

            kpt_other = TS kpt_ref + G0
    """
    ref_gprimd_inv = np.inv(np.asarray(ref_lattice).T)
    other_gprimd = np.asarray(other_lattice).T
    other_kpoints = np.asarray(other_kpoints).reshape((-1, 3))
    ref_kpoints = np.asarray(ref_kpoints).reshape((-1, 3))
    o2r_map = len(other_kpoints) * [None]

    tsigns = (1, -1) if has_timrev else (1,)
    kmap = collections.namedtuple("kmap", "ik_ref, tsign, isym, g0")

    for ik_oth, okpt in enumerate(other_kpoints):
        # Get other k-point in reduced coordinates in the referece lattice.
        okpt_red = np.matmul(ref_gprimd_inv, np.matmul(other_gprimd, okpt))

        # k_other = TS k_ref + G0
        found = False
        for ik_ref, kref in enumerate(ref_kpoints):
            if found: break
            for tsign in tsign:
                for isym, symrec in enumerate(ref_symrecs):
                    krot = tsign * np.matmul(symrec, kref)
                    if issamek(okpt_red, krot):
                        g0 = np.rint(okpt_red - krot)
                        o2r_map[ik_oth] = kmap(ik_ref, tsign, isym, g0)
                        found = True
                        break

        return o2r_map, o2r_map.count(None)


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
            assert all( isinstance(o, Kpoint) for o in obj)
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


class Kpoint(SlotPickleMixin):
    """
    Class defining one k-point.
    """

    __slots__ = [
        "_frac_coords",
        "_lattice",
        "_weight",
        "_name",
        "_hash",
    ]

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
        return np.any(diff < _ATOL_KDIFF)

    def __repr__(self):
        return "[%.3f, %.3f, %.3f]" % tuple(self.frac_coords)

    def __str__(self):
        s =  "[%.3f, %.3f, %.3f]" % tuple(self.frac_coords)
        if self.name is not None: s += ", name=%s" % self.name
        if self._weight is not None: s += ", weight=%.3f" % self.weight
        return s

    # Kpoint algebra.
    def __add__(self, other):
        return self.__class__(self.frac_coords + other.frac_coords, self.lattice)

    def __sub__(self, other):
        return self.__class__(self.frac_coords - other.frac_coords, self.lattice)

    def __eq__(self, other):
        try:
            # Comparison between two Kpoint objects
            return issamek(self.frac_coords, other.frac_coords)
        except AttributeError:
            # Kpoint vs iterable (e.g. list)
            return issamek(self.frac_coords, other)

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
                if issamek(sk_coords, prev_coords): break
            else:
                frac_coords.append(sk_coords)

        return KpointStar(self.lattice, frac_coords, weights=None, names=len(frac_coords) * [self.name])


class KpointList(collections.Sequence):
    """
    Base class defining a sequence of :class:`Kpoint` objects. Essentially consists
    of base methods implementing the sequence protocol and helper functions.
    The subclasses `Kpath` and `IrredZone` provide specialized methods to operate
    on k-points representing a path or list of points in the IBZ, respectively.

    .. Note:

        Algorithms usually need to know what kind of sampling we are using.
        The test can be easily implemented with:

        if kpoints.is_path:
            # code specific to k-paths.

        elif kpoints.is_ibz:
            # code specific to IBZ sampling.

    """
    Error = KpointsError

    @classmethod
    def subclass_from_name(cls, name):
        """Return the class with the given name."""
        if cls.__name__ == name: return c
        for c in cls.__subclasses__():
            if c.__name__ == name: return c

        raise ValueError("Cannot find subclass associated to name: %s" % name)

    @classmethod
    def from_dict(cls, d):
        """
        Makes Kpoints obey the general json interface used in pymatgen for easier serialization.
        """
        reciprocal_lattice = Lattice.from_dict(d["reciprocal_lattice"])
        return cls(reciprocal_lattice, d["frac_coords"],
                   weights=d["weights"], names=d["names"], ksampling=d["ksampling"])

    @pmg_serialize
    def as_dict(self):
        """
        Makes Kpoints obey the general json interface used in pymatgen for easier serialization.
        """
        if self.weights is not None: weights = self.weights.tolist()
        return dict(
            reciprocal_lattice=self.reciprocal_lattice.as_dict(),
            frac_coords=self.frac_coords.tolist(),
            weights=weights,
            names=[k.name for k in self],
            ksampling=self.ksampling,
        )

    def __init__(self, reciprocal_lattice, frac_coords, weights=None, names=None, ksampling=None):
        """
        Args:
            reciprocal_lattice: :class:`Lattice` object.
            frac_coords: Array-like object with the reduced coordinates of the k-points.
            weights: List of k-point weights. If None, weights are initialized with zeros.
            names: List of k-point names.
            ksampling: Info on the k-point sampling (used for homogeneous meshes)
        """
        self._reciprocal_lattice = reciprocal_lattice
        self._frac_coords = frac_coords = np.reshape(frac_coords, (-1, 3))
        self.ksampling = ksampling

        if weights is not None:
            if len(weights) != len(frac_coords):
                raise ValueError("len(weights) != len(frac_coords):\nweights: %s\nfrac_coords: %s" %
                    (weights, frac_coords))
            weights = np.asarray(weights)
        else:
            weights = np.zeros(len(self.frac_coords))

        if names is not None and len(names) != len(frac_coords):
            raise ValueError("len(names) != len(frac_coords):\nnames: %s\nfrac_coords: %s" %
                    (names, frac_coords))

        self._points = []
        for i, rcs in enumerate(frac_coords):
            name = None if names is None else names[i]
            self._points.append(Kpoint(rcs, self.reciprocal_lattice, weight=weights[i], name=name))

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from file."""
        if filepath.endswith(".nc"):
            with KpointsReader(filepath) as r:
                new = r.read_kpoints()
        else:
            raise NotImplementedError("Only netcdf files are supported.")

        new.__class__ = cls
        return new

    @property
    def reciprocal_lattice(self):
        """`Lattice` object defining the reciprocal lattice."""
        return self._reciprocal_lattice

    def __repr__(self):
        return "\n".join("%d) %s" % (i, repr(kpoint)) for i, kpoint in enumerate(self))

    def __str__(self):
        return "\n".join("%d) %s" % (i, str(kpoint)) for i, kpoint in enumerate(self))

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
        if self.reciprocal_lattice != other.reciprocal_lattice:
            raise ValueError("Cannot merge k-points with different reciprocal lattice.")

        return KpointList(self.reciprocal_lattice,
                          frac_coords=[k.frac_coords for k in self] + [k.frac_coords for k in other],
                          weights=None,
                          names=[k.name for k in self] + [k.name for k in other],
                        )

    def __eq__(self, other):
        if other is None or not isinstance(other, KpointList): return False
        for k1, k2 in zip(self, other):
            if k1 != k2: return False
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
            raise ValueError("Cannot find point: %s in KpointList:\n%s" % (repr(kpoint), repr(self)))

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

    def find_closest(self, obj):
        """
        Find the closest k-point in the list (not necessarily equal).

        Args:
            obj: Fractional coordinates or :class:`Kpoint` instance.

        Return:
            (ind, kpoint, dist)

            where `ind` is the index in self of the closest k-point.
            `kpoint` is the :class:`Kpoint` instance of index `ind`.
            dist is the distance between `obj` and `kpoint`.
        """
        if isinstance(obj, Kpoint):
            if obj.lattice != self.reciprocal_lattice:
                raise ValueError("Kpoint list and Kpoint object have different lattices!")
            frac_coords = obj.frac_coords
        else:
            frac_coords = np.asarray(obj)

        dist = np.empty(len(self))
        for i, kpt in enumerate(self):
            dist[i] = kpt.lattice.norm(kpt.frac_coords - frac_coords)

        ind = dist.argmin()
        return ind, self[ind], np.copy(dist[ind])

    @property
    def is_path(self):
        """True if self represents a path in the BZ."""
        return isinstance(self, Kpath)

    @property
    def is_ibz(self):
        """True if self represents a list of points in the IBZ."""
        return isinstance(self, IrredZone)

    @lazy_property
    def mpdivs_shifts(self):
        """
        The Monkhorst-Pack (MP) divisions and shifts.
        Both quantities are set to None if self is not a MP mesh.
        Use `is_mpmesh` to check whether self is a MP mesh.
        """
        ksampling = self.ksampling
        if not self.is_ibz: return (None, None)
        m = ksampling.kptrlatt.copy()
        np.fill_diagonal(m, 0)
        if np.any(m != 0): return (None, None)
        return ksampling.kptrlatt.diagonal(), ksampling.shifts

    @property
    def is_mpmesh(self):
        """
        True if self represents a Monkhorst-Pack mesh.
        i.e if the sampling has been specified in terms of divisions
        along the reciprocal lattice vectors (ngkpt)
        """
        return self.mpdivs_shifts[0] is not None

    @property
    def frac_coords(self):
        """Fractional coordinates of the k-point as `ndarray` of shape (len(self), 3)"""
        return self._frac_coords

    @property
    def names(self):
        """List with the name of the k-points."""
        return [k.name for k in self]

    @property
    def weights(self):
        """`ndarray` with the weights of the k-points."""
        return np.array([kpoint.weight for kpoint in self])

    def sum_weights(self):
        """Returns the sum of the weights."""
        return np.sum(self.weights)

    def remove_duplicated(self):
        """
        Remove duplicated k-points from self. Returns new KpointList instance.
        """
        frac_coords, good_indices = [self[0].frac_coords], [0]

        for i, kpoint in enumerate(self[1:]):
            i += 1
            # Add it only if it's not already in the list.
            for prev_coords in frac_coords:
                if issamek(kpoint.frac_coords, prev_coords): break
            else:
                frac_coords.append(kpoint.frac_coords)
                good_indices.append(i)

        good_kpoints = [self[i] for i in good_indices]

        return self.__class__(
                self.reciprocal_lattice,
                frac_coords=[k.frac_coords for k in good_kpoints],
                weights=None,
                names=[k.name for k in good_kpoints],
                ksampling=self.ksampling)

    def to_array(self):
        """Returns a `ndarray` [nkpy, 3] with the fractional coordinates."""
        return np.array(self.frac_coords.copy())

    def to_json(self):
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    def plot(self, ax=None, **kwargs):
        from pymatgen.electronic_structure.plotter import plot_brillouin_zone
        fold = False

        if self.is_path:
            labels = {k.name: k.frac_coords for k in self if k.name}
            frac_coords_lines = [self.frac_coords[line] for line in self.lines]
            return plot_brillouin_zone(self.reciprocal_lattice, lines=frac_coords_lines, labels=labels,
                                       ax=ax, fold=fold, **kwargs)
        else:
            # Not sure this works, I got points outside of the BZ in a simple with Si and Gamm-centered 8x8x8.
            # Don't know if it's a bug in matplotlib or plot_brillouin_zone.
            #print(self.frac_coords)
            return plot_brillouin_zone(self.reciprocal_lattice, kpoints=self.frac_coords,
                                       ax=ax, fold=fold, **kwargs)


class KpointStar(KpointList):
    """
    Star of the kpoint. Note that the first k-point is assumed to be the base
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
    """
    This object describes a k-path in reciprocal space.
    It provides methods to compute (line) derivatives along the path.
    The k-points do not have weights so Kpath should not be used to compute integral in the BZ.
    """

    @classmethod
    def from_vertices_and_names(cls, structure, vertices_names, line_density=20):
        """
        Generate K-path from a list of vertices and the corresponding labels.

        Args:
            structure: Structure object.
            vertices_names:  List of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path.
            line_density: Number of points used to sample the smallest segment of the path
        """
        gmet = structure.lattice.reciprocal_lattice.metric_tensor
        vnames = [str(vn[1]) for vn in vertices_names]
        vertices = np.array([vn[0] for vn in vertices_names], dtype=np.float)
        vertices.shape = (-1, 3)

        dl_vals = []
        for ik, k0 in enumerate(vertices[:-1]):
            dk = vertices[ik + 1] - k0
            dl = np.sqrt(np.dot(dk, np.matmul(gmet, dk)))
            if abs(dl) < 1e-6: dl = np.inf
            dl_vals.append(dl)

        dl_min = np.array(dl_vals).min()

        fact = dl_min / line_density
        frac_coords = collections.deque()
        knames = collections.deque()
        for ik, dl in enumerate(dl_vals):
            if dl == np.inf: continue
            numk = int(np.rint(dl / fact))
            k0 = vertices[ik]
            dk = vertices[ik + 1] - k0
            knames.append(vnames[ik])
            for ii in range(numk):
                next_k = k0 + dk * ii / numk
                frac_coords.append(next_k)
                if ii > 0: knames.append("")
        knames.append(vnames[-1])
        frac_coords.append(vertices[-1])

        return cls(structure.lattice.reciprocal_lattice,
                   frac_coords=frac_coords,
                   weights=None,
                   names=knames,
                   )

    def __str__(self):
        return self.to_string()

    def to_string(self):
        lines = []
        app = lines.append
        app("K-path contains %s lines. Number of k-points in each line: %s" % (
            len(self.lines), [len(l) for l in self.lines]))
        #for i, line in enumerate(self.lines):
        #    app("line %d: %s" % (i, line))
        header = "\n".join(lines)

        vids = {line[0] for line in self.lines}

        table = [["Idx", "Frac_coords", "Name", "ds", "Vert",]]
        for i, kpoint in enumerate(self):
            table.append([
                str(i),
                "%.5f, %.5f, %.5f" % tuple(kpoint.frac_coords),
                kpoint.name,
                self.ds[i] if i != len(self) - 1 else None,
                "*" if i in vids else " ",
            ])
        return "\n".join([header, " ", tabulate(table, headers="firstrow")])

    @lazy_property
    def ds(self):
        """
        numpy array of len(self)-1 elements giving the distance between two
        consecutive k-points, i.e. ds[i] = ||k[i+1] - k[i]|| for i=0,1,...,n-1
        """
        ds = np.zeros(len(self) - 1)
        for i, kpoint in enumerate(self[:-1]):
            ds[i] = (self[i + 1] - kpoint).norm
        return ds

    @lazy_property
    def versors(self):
        """
        Tuple of len(self)-1 elements with the versors connecting k[i] to k[i+1].
        """
        versors = (len(self) - 1) * [None, ]
        for i, kpt in enumerate(self[:-1]):
            versors[i] = (self[i + 1] - kpt).versor()
        return tuple(versors)

    @lazy_property
    def lines(self):
        """
        Nested list containing the indices of the points belonging to the same line.
        Used for extracting the eigenvalues while looping over the lines.

        Example:

            for line in self.lines:
                vals_on_line = eigens[spin, line, band]
        """
        prev = self.versors[0]
        lines = [[0]]

        for i, v in enumerate(self.versors[1:]):
            i += 1
            if v != prev:
                #print("diff", v.frac_coords - prev.frac_coords)
                prev = v
                lines[-1].append(i)
                lines.append([i])
            else:
                lines[-1].append(i)

        lines[-1].append(len(self)-1)
        return tuple(lines)

    def finite_diff(self, values, order=1, acc=4):
        """
        Compute the derivatives of values by finite differences.

        Args:
            values: array-like object with shape=(nkpt) containing the values of the path.
            order: Order of the derivative.
            acc: Accuracy: 4 corresponds to a central difference with 5 points.

        Returns:
            ragged numpy array. The i-th entry is a numpy array with the derivatives
            on the i-th line.
        """
        values = np.asarray(values)
        assert len(values) == len(self)

        # Loop over the lines of the path, extract the data on the line and
        # differentiate f(s) where s is the distance between two consecutive points along the line.
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

    .. note::

            Abinit supports different options for the specification of the BZ sampling:

                 - kptrlatt(3,3) or ngkpt(3) for the definition grid.
                 - shiftk(3, nshiftk) for the definition of multiple shifts.
                 - `kptopt` for the treatment of symmetry operations.

            All these possibilities complicate the internal implementation in particular when
            we need to recostruct the full BZ and take into account the presence of multiple shifts
            since kptrlatt may have non-zero off-diagonal components. Client code that needs to know
            how the mesh was generated can rely on the following checks:

            if not self.ibz: raise("Need an IBZ sampling")

            mpdivs, shifts = self.mpdivs_shifts
            if mpdivs is None: raise ValueError("Cannot handle kptrlatt with non-zero off-diagonal elements")
            if len(shifts) > 1: raise ValueError("Multiple shifts are not supported")
            # Code for mesh defined in terms of mpdivs and one shift.
    """
    def __init__(self, reciprocal_lattice, frac_coords, weights=None, names=None, ksampling=None):
        """
        Args:
            reciprocal_lattice: :class:`Lattice` object
            frac_coords: Array-like object with the reduced coordinates of the points.
            weights: Array-like with the weights of the k-points.
            names: List with the name of the k-points.
            ksampling: Info on the k-point sampling
        """
        super(IrredZone, self).__init__(reciprocal_lattice, frac_coords,
                                        weights=weights, names=names, ksampling=ksampling)

        # Weights must be normalized to one.
        wsum = self.sum_weights()
        if abs(wsum - 1) > 1.e-6:
            err_msg = ("The list of kpoints does not represent a homogeneous sampling of the BZ\n"
                       "Kpoint weights should sum up to one while sum_weights is %.3f\n" % wsum)
            print(err_msg)
            #raise ValueError(err_msg)

    def __str__(self):
        return self.to_string()

    def to_string(self):
        lines = []; app = lines.append

        if self.is_mpmesh:
            mpdivs, shifts = self.mpdivs_shifts
            d = "[%d, %d, %d]" % tuple(mpdivs)
            s = ", ".join("[%.1f, %.1f, %.1f]" % tuple(s) for s in shifts)
            app("K-mesh with divisions: %s, shifts: %s, kptopt: %s" % (d, s, self.ksampling.kptopt))
        else:
            for k, v in self.ksampling.items():
                app("%s: %s" % (k, v))

        return "\n".join(lines)

    #@property
    #def len_bz(self):
    #    """Number of points in the full BZ."""
    #    return self.mpdivs.prod() * self.num_shifts

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
        return self.mpdivs is not None or self.kptrlatt is not None and self.kptopt > 0

    #@property
    #def is_path(self):
    #    """True if we have a path in the BZ."""
    #    return self.kptopt < 0


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

        if ksampling.kptopt < 0:
            # We have a path in the BZ.
            kpath = Kpath(structure.reciprocal_lattice, frac_coords, ksampling=ksampling)
            for kpoint in kpath:
                kpoint.set_name(structure.findname_in_hsym_stars(kpoint))
            return kpath

        # FIXME
        # Quick and dirty hack to allow the reading of the k-points from WFK files
        # where info on the sampling is missing. I will regret it but at present
        # is the only solution I found (changes in the ETSF-IO part of Abinit are needed)
        #if ksampling.is_homogeneous or abs(sum(weights) - 1.0) < 1.e-6:
        #if np.any(ksampling.kptrlatt_orig != 0) or abs(sum(weights) - 1.0) < 1.e-6:

        #if np.any(ksampling.kptrlatt_orig != 0):
        # We have a homogeneous sampling of the BZ.
        return IrredZone(structure.reciprocal_lattice, frac_coords, weights=weights, ksampling=ksampling)

        raise ValueError(
          "Only homogeneous samplings or paths are supported!"
          "ksamping info:\n%s" % str(ksampling))

    def read_ksampling_info(self):
        # FIXME: in v8.0, the SIGRES files does not have kptopt, kptrlatt_orig and shiftk_orig
        kptrlatt = self.read_kptrlatt()
        shifts = self.read_kshifts()

        return KSamplingInfo(
            kptopt=int(self.read_value("kptopt", default=0)),
            mpdivs=self.read_kmpdivs(),
            kptrlatt=kptrlatt,
            shifts=shifts,
            kptrlatt_orig=self.read_value("kptrlatt_orig", default=kptrlatt),
            shiftk_orig=self.read_value("shiftk_orig", default=shifts),
        )

    def read_kfrac_coords(self):
        """Fractional coordinates of the k-points"""
        return self.read_value("reduced_coordinates_of_kpoints")

    #@returns_None_onfail
    def read_kweights(self):
        """Returns the weight of the k-points. None if not found."""
        return self.read_value("kpoint_weights")

    #@returns_None_onfail
    def read_kshifts(self):
        """Returns the shifts of the k-mesh in reduced coordinates. None if not found."""
        try:
            return self.read_value("shiftk")
        except self.Error:
            return self.read_value("kpoint_grid_shift")

    #@returns_None_onfail
    def read_kmpdivs(self):
        """Returns the Monkhorst-Pack divisions defining the mesh. None if not found."""
        return self.read_value("monkhorst_pack_folding")

    #@returns_None_onfail
    def read_kptrlatt(self):
        """Returns ABINIT variable kptrlatt. None if not found."""
        try:
            return self.read_value("kptrlatt")
        except self.Error:
            return self.read_value("kpoint_grid_vectors")

    #@returns_None_onfail
    #def read_kptopt(self):
    #    """Returns the ABINIT variable kptopt. None if not found."""
    #    return int(self.read_value("kptopt"))


class KpointsReader(ETSF_Reader, KpointsReaderMixin):
    """This object reads k-point data from a netcdf file."""


class Ktables(object):
    """
    Call spglib to compute the k-points in the IBZ with the corresponding weights
    as well as the mapping BZ --> IBZ.

    Args:
        mesh:
        is_shift:

    Attributes:

        mesh
        is_shift
        ibz:
        nibz
        weights:
        bz:
        nbz
        grid:
    """
    def __init__(self, structure, mesh, is_shift, has_timrev):
        """

        Args:
            structure
            mesh
            is_shift
            has_timrev
        """
        import spglib as spg
        self.mesh = np.array(mesh)
        self.is_shift = is_shift
        self.has_timrev = has_timrev
        cell = (structure.lattice.matrix, structure.frac_coords, structure.atomic_numbers)

        mapping, self.grid = spg.get_ir_reciprocal_mesh(self.mesh, cell,
            is_shift=self.is_shift, is_time_reversal=self.has_timrev, symprec=_SPGLIB_SYMPREC)

        uniq, self.weights = np.unique(mapping, return_counts=True)
        self.weights = np.asarray(self.weights, dtype=np.float) / len(self.grid)
        self.nibz = len(uniq)
        self.kshift = [0., 0., 0.] if is_shift is None else 0.5 * np.asarray(is_shift)
        self.ibz = (self.grid[uniq] + self.kshift) / self.mesh
        self.bz = (self.grid + self.kshift) / self.mesh
        self.nbz = len(self.bz)

        # All k-points and mapping to ir-grid points.
        # FIXME This is slow.
        self.bz2ibz = np.empty(len(self.bz), dtype=np.int)
        for ik_bz, ir_gp_id in enumerate(mapping):
            inds = np.where(uniq == ir_gp_id)
            assert len(inds) == 1
            self.bz2ibz[ik_bz] = inds[0]

    def __str__(self):
        return self.to_string()

    def to_string(self, **kwargs):
        """String representation"""
        lines = collections.deque()
        app = lines.append
        app("mesh %s, shift %s, time-reversal: %s, Irred points: %d" % (
            self.mesh, self.kshift, self.has_timrev, self.nibz))

        app("Irreducible k-points with number of points in star:\n")
        for ik, kpt in enumerate(self.ibz):
            app("%s: [%9.6f, %9.6f, %9.6f], nstar: %d" % (ik, kpt[0], kpt[1], kpt[2], self.weights[ik] * self.nbz))

        return "\n".join(lines)

    def print_bz2ibz(self, file=sys.stdout):
        """Print BZ --> IBZ mapping."""
        print("BZ points --> IBZ points", file=file)
        for ik_bz, ik_ibz in enumerate(self.bz2ibz):
            print("%6d) [%9.6f, %9.6f, %9.6f], ====> %6d) [%9.6f, %9.6f, %9.6f]," %
                (ik_bz, self.bz[ik_ibz][0], self.bz[ik_ibz][1], self.bz[ik_ibz][2],
                ik_ibz, self.ibz[ik_ibz][0], self.ibz[ik_ibz][1], self.ibz[ik_ibz][2]), file=file)
