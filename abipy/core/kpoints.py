# coding: utf-8
"""This module defines objects describing the sampling of the Brillouin Zone."""
from __future__ import print_function, division, unicode_literals, absolute_import

import collections
import json
import sys
import time
import numpy as np

from itertools import product
from tabulate import tabulate
from monty.json import MSONable, MontyEncoder
from monty.collections import AttrDict, dict2namedtuple
from monty.functools import lazy_property
from monty.termcolor import cprint
from monty.string import marquee
from pymatgen.core.lattice import Lattice
try:
    from pymatgen.util.serialization import pmg_serialize
except ImportError:
    from pymatgen.serializers.json_coders import pmg_serialize
try:
    from pymatgen.util.serialization import SlotPickleMixin
except ImportError:
    from pymatgen.serializers.pickle_coders import SlotPickleMixin
from abipy.iotools import ETSF_Reader
from abipy.tools.derivatives import finite_diff
from abipy.tools.numtools import add_periodic_replicas, is_diagonal

import logging
logger = logging.getLogger(__name__)

__all__ = [
    "issamek",
    "wrap_to_ws",
    "wrap_to_bz",
    "as_kpoints",
    "Kpoint",
    "KpointList",
    "KpointStar",
    "Kpath",
    "IrredZone",
    "rc_list",
    "kmesh_from_mpdivs",
    "Ktables",
    "find_points_along_path",
]

# Tolerance used to compare k-points.
_ATOL_KDIFF = 1e-8

# Tolerances passed to spglib.
_SPGLIB_SYMPREC = 1e-5
_SPGLIB_ANGLE_TOLERANCE = -1.0


def set_atol_kdiff(new_atol):
    """
    Change the value of the tolerance ``_ATOL_KDIFF`` used to compare k-points.
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
    Change the value of the tolerances ``symprec`` and ``angle_tolerance``
    used to call spglib_. Return old values

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


def issamek(k1, k2, atol=None):
    """
    True if k1 and k2 are equal modulo a lattice vector.
    Use _ATOL_KDIFF is atol is None.

    >>> assert issamek([1, 1, 1], [0, 0, 0])
    >>> assert issamek([1.1, 1, 1], [0, 0, 0], atol=0.1)
    >>> assert not issamek(0.00003, 1)
    """
    k1 = np.asarray(k1)
    k2 = np.asarray(k2)
    #if k1.shape != k2.shape:

    return is_integer(k1 - k2, atol=atol)


def wrap_to_ws(x):
    """
    Transforms x in its corresponding reduced number in the interval ]-1/2,1/2].
    """
    w = x % 1
    return np.where(w > 0.5, w - 1.0, w)


def wrap_to_bz(x):
    """
    Transforms x in its corresponding reduced number in the interval [0,1[."
    """
    return x % 1


def rc_list(mp, sh, pbc=False, order="bz"):
    """
    Returns a |numpy-array| with the linear mesh used to sample one dimension of the reciprocal space.
    Note that the coordinates are always ordered so that rc[i+1] > rc[i].
    so that we can easily plot quantities defined on the 3D multidimensional mesh.

    Args:
        mp: Number of Monkhorst-Pack divisions along this direction.
        sh: Shift
        pbc: if pbc is True, periodic boundary conditions are enforced.
        order: Possible values ["bz", "unit_cell"].
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
    Returns a |numpy-array| with the reduced coordinates of the
    k-points from the MP divisions and the shifts.

    Args:
        mpdivs: The three MP divisions.
        shifts: Array-like object with the MP shift.
        pbc: If True, periodic images of the k-points will be included i.e. closed mesh.
        order: "unit_cell" if the kpoint coordinates must be in [0,1)
               "bz" if the kpoint coordinates must be in [-1/2, +1/2)
    """
    shifts = np.reshape(shifts, (-1, 3))
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


def map_grid2ibz(structure, ibz, ngkpt, has_timrev, pbc=False):
    """
    Compute the correspondence between a *grid* of k-points in the *unit cell*
    associated to the ``ngkpt`` mesh and the corresponding points in the IBZ.
    Requires structure with Abinit symmetries.
    This routine is mainly used to symmetrize eigenvalues in the unit cell
    e.g. to write BXSF files for electronic isosurfaces.

    Args:
        structure: Structure with (Abinit) symmetry operations.
        ibz: [*, 3] array with reduced coordinates in the in the IBZ.
        ngkpt: Mesh divisions.
        has_timrev: True if time-reversal can be used.
        pbc: True if the mesh should contain the periodic images (closed mesh).

    Returns:
        bz2ibz: 1d array with BZ --> IBZ mapping
    """
    ngkpt = np.asarray(ngkpt, dtype=np.int)

    # Extract (FM) symmetry operations in reciprocal space.
    abispg = structure.abi_spacegroup
    if abispg is None:
        raise ValueError("Structure does not contain Abinit spacegroup info!")

    # Extract rotations in reciprocal space (FM part).
    symrec_fm = [o.rot_g for o in abispg.fm_symmops]

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

    if pbc:
        # Add periodic replicas.
        bzgrid2ibz = add_periodic_replicas(bzgrid2ibz)

    if np.any(bzgrid2ibz == -1):
        #for ik_bz, ik_ibz in enumerate(self.bzgrid2ibz): print(ik_bz, ">>>", ik_ibz)
        msg = "Found %s/%s invalid entries in bzgrid2ibz array" % ((bzgrid2ibz == -1).sum(), bzgrid2ibz.size)
        msg += "This can happen if there an inconsistency between the input IBZ and ngkpt"
        msg += "ngkpt: %s, has_timrev: %s" % (str(ngkpt), has_timrev)
        raise ValueError(msg)

    bz2ibz = bzgrid2ibz.flatten()
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
    True if time-reversal symmetry can be used to generate k-points in the IBZ.
    """
    # note: We assume TR if negative value i.e. band structure k-sampling.
    return int(kptopt) not in (3, 4)


def kptopt2str(kptopt, verbose=0):
    """
    Return human-readable string with meaning of kptopt.
    """
    if kptopt < 0:
        t = ("Band structure run. Use kptbounds, and ndivk (ndivsm)"
             "The absolute value of kptopt gives the number of segments of the band structure." +
             "Weights are usually irrelevant with this option")
    else:
        t = {
            0: ("Manual mode",
                "User-provided nkpt, kpt, kptnrm and wtk"),
            1: ("Use space group symmetries and TR symmetry",
                "Usual mode for GS calculations (ngkpt or kptrlatt, nshiftk and shiftk)"),
            2: ("Only TR symmetry",
                "This is to be used for DFPT at Gamma (ngkpt or kptrlatt, nshiftk and shiftk)"),
            3: ("Do not take into account any symmetry",
                "This is to be used for DFPT at non-zero q (ngkpt or kptrlatt, nshiftk and shiftk)."),
            4: ("Spatial symmetries, NO TR symmetry",
                "This has to be used for PAW calculations with SOC (pawspnorb/=0) " +
                "from ngkpt or kptrlatt, nshiftk and shiftk."),
        }[kptopt]

    return t[0] if verbose == 0 else t[0] + "\n" + t[1]


def map_kpoints(other_kpoints, other_lattice, ref_lattice, ref_kpoints, ref_symrecs, has_timrev):
    """
    Build mapping between a list of k-points in reduced coordinates (``other_kpoints``)
    in the reciprocal lattice ``other_lattice`` and a list of reference k-points given
    in the reciprocal lattice `ref_lattice` with symmetry operations ``ref_symrecs``.

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


#def find_irred_kpoints_kmesh(structure, kfrac_coords):
#    """
#    Remove k-points that are connected to each other by one of the
#    symmetry operations of the space group. Assume k-points
#    belonging to a homogeneous mesh.
#
#    Args:
#        structure: Structure object.
#        kfrac_coords: Reduced coordinates of the k-points.
#
#    Return:
#    """
#    # Wrap in [0,1[ interval.
#    uc_kcoords = np.reshape(kfrac_coords, (-1, 3)) % 1
#    numk = len(uc_kcoords)
#    nx, ny, nz = np.int(np.floor(1 / uc_kcoords.min(axis=0)))
#
#    # Compute rank and invrank
#    rank = np.array(numk, dtype=np.int)
#    invrank = {}
#    for ik, kk in enumerate(uc_kcoords):
#        rk = iz + iy * nz + ix * ny * nz
#        rank[ik] = rk
#        invrank[rank] = ik
#
#    irred_map = collections.deque()
#    irred_map.append(0)
#    kpts2irred = collections.deque()
#    kpts2irred.append((0, 0, +1))
#
#    for ik, kk in enumerate(uc_kcoords[1:]):
#        ik += 1
#        found = False
#        for ik_irr in irred_map:
#            kirr = kfrac_coords[ik_irr]
#            for isym, symmop in enumerate(structure.abi_spacegroup):
#                krot = symmop.rotate_k(kirr)
#                new_frac_coords = krot.frac_coords % 1
#                if issamek(krot, kk):
#                    #kpts2irred[ik] = ik_irr
#                    #kpts2irred[ik] = isym
#                    found = True
#                    break
#
#    #return irred_map


def find_irred_kpoints_generic(structure, kfrac_coords, verbose=1):
    """
    Remove the k-points that are connected to each other by one of the
    symmetry operations of the space group. No assumption is done
    on the initial k-point sampling, this means that one can call this
    function to treat points on a path in reciprocal space.

    Args:
        structure: |Structure| object.
        kfrac_coords: Reduced coordinates of the k-points.

    Return:
        irred_map: Index of the i-th irreducible k-point in the input kfrac_coords array.

    .. warning::

        In the worst case, the algorithm scales as nkpt ** 2 * nsym.
        hence this routine should be used only if ``kfrac_coords`` represents
        e.g. a path in the Brillouin zone or an arbitrary set of points.
    """
    start = time.time()
    print("Removing redundant k-points. This is gonna take a while... ")

    # Wrap points in [0,1[ interval.
    uc_kcoords = np.reshape(kfrac_coords, (-1, 3)) % 1

    irred_map = collections.deque()
    irred_map.append(0)
    kpts2irred = collections.deque()
    kpts2irred.append((0, 0, +1))

    for ik, kk in enumerate(uc_kcoords[1:]):
        ik += 1
        found = False
        for ik_irr in irred_map:
            kirr = kfrac_coords[ik_irr]
            for isym, symmop in enumerate(structure.abi_spacegroup):
                krot = symmop.rotate_k(kirr)
                if issamek(krot, kk):
                    found = True
                    #kpts2irred[ik] = (ik_irr, isym, symmops.time_sign)
                    break

        if not found:
            irred_map.append(ik)

    print("Completed in", time.time() - start, "[s]")
    if verbose:
        print("Entered with ", len(uc_kcoords), "k-points")
        print("Found ", len(irred_map), "irred k-points")

    return dict2namedtuple(irred_map=np.array(irred_map, dtype=np.int))


def kpath_from_bounds_and_ndivsm(bounds, ndivsm, structure):
    """
    Generate a normalized path given the extrema and the number of divisions for the smallest segment

    Args:
        bounds: (N, 3) array with the boundaries of the path in reduced coordinates.
        ndivsm: Number of divisions used for the smallest segment.

    Return:
        Array (M, 3) with fractional coordinates.
    """
    bounds = np.reshape(bounds, (-1, 3))
    nbounds = len(bounds)
    if nbounds == 1:
        raise ValueError("Need at least two points to define the k-path!")

    lens = []
    for i in range(nbounds - 1):
        v = bounds[i + 1] - bounds[i]
        lens.append(float(structure.reciprocal_lattice.norm(v)))

    # Avoid division by zero if any bounds[i+1] == bounds[i]
    minlen = np.min(lens)
    if minlen < 1e-6:
        raise ValueError("Found two equivalent consecutive points in bounds!")

    minlen = minlen / ndivsm
    ndivs = np.rint(lens / minlen).astype(np.int)
    path = []
    for i in range(nbounds - 1):
        for j in range(ndivs[i]):
            p = bounds[i] + j * (bounds[i + 1] - bounds[i]) / ndivs[i]
            path.append(p)
    path.append(bounds[-1])

    return np.array(path)


class KpointsError(Exception):
    """Base error class for KpointList exceptions."""


def as_kpoints(obj, lattice, weights=None, names=None):
    """
    Convert obj into a list of k-points.

    Args:
        obj: :class:`Kpoint` or list of Kpoint objects or array-like object.
        lattice: Reciprocal lattice.
        weights: k-point weights. Ignored if obj is already a `Kpoint` instance or a list of `Kpoint` items.
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
    Class defining one k-point. This object is immutable and can be used as key in dictionaries

    Note that we usually construct the object by passing pymatgen.reciprocal_lattice
    that is the standard reciprocal lattice used for solid state physics 
    with a factor of 2 * pi i.e. a_i . b_j = 2pi delta_ij. 
    Abinit, on the contrary, uses the crystallographic reciprocal lattice i.e. no 2pi factor.
    so pay attention when converting Abinit routines to AbiPy.
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
            lattice: |Lattice| object describing the reciprocal lattice.
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
        """|Lattice| object describing the Reciprocal lattice."""
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

    @lazy_property
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

    @lazy_property
    def on_border(self):
        """
        True if the k-point is on the border of the BZ (lattice translations are taken into account).
        """
        kreds = wrap_to_ws(self.frac_coords)
        diff = np.abs(np.abs(kreds) - 0.5)
        return np.any(diff < _ATOL_KDIFF)

    def __repr__(self):
        s = "[%+.3f, %+.3f, %+.3f]" % tuple(self.frac_coords)
        if self.name is not None:
            s += " %s" % self.name
        return s

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        s =  "[%+.3f, %+.3f, %+.3f]" % tuple(self.frac_coords)
        if self.name is not None:
            s += ", name: %s" % self.name
        if self._weight is not None: s += ", weight: %.3f" % self.weight
        return s

    # Kpoint algebra.
    def __add__(self, other):
        return self.__class__(self.frac_coords + other.frac_coords, self.lattice)

    def __sub__(self, other):
        return self.__class__(self.frac_coords - other.frac_coords, self.lattice)

    def __eq__(self, other):
        if hasattr(other, "frac_coords"):
            # Comparison between two Kpoint objects
            return issamek(self.frac_coords, other.frac_coords)
        else:
            # Kpoint vs iterable (e.g. list)
            return issamek(self.frac_coords, other)

    def __ne__(self, other):
        return not (self == other)

    def __getitem__(self, slice):
        return self.frac_coords[slice]

    @classmethod
    def as_kpoint(cls, obj, lattice):
        """
        Convert obj into a :class:`Kpoint` instance.

        Args:
            obj: :class:`Kpoint` instance or array-like with the reduced coordinates.
            lattice: |Lattice| object defining the reciprocal lattice.
        """
        if isinstance(obj, cls):
            return obj
        else:
            return cls(obj, lattice, weight=None, name=None)

    @classmethod
    def gamma(cls, lattice, weight=None):
        """Constructor for the Gamma point."""
        return cls(np.zeros(3), lattice, weight=weight, name=r"$\Gamma$")

    def copy(self):
        """Deep copy."""
        return self.__class__(self.frac_coords.copy(), self.lattice.copy(),
                              weight=self.weight, name=self.name)

    def is_gamma(self, allow_umklapp=False, atol=None):
        """
        Return True if this the gamma point.

        Args:
            allow_umklapp: True if umklapp G-vectors are allowed.
            atol: Absolute tolerance for k-point comparison (used only if umklapp).
        """
        if not allow_umklapp:
            return np.all(self.frac_coords == 0.0)
        else:
            return issamek(self.frac_coords, [0, 0, 0], atol=atol)

    @lazy_property
    def norm(self):
        """Norm of the kpoint."""
        return np.sqrt(np.dot(self.cart_coords, self.cart_coords))

    def versor(self):
        """Returns the versor i.e. math:`||k|| = 1`"""
        cls = self.__class__
        if self.norm > 1e-12:
            return cls(self.frac_coords / self.norm, self.lattice, weight=self.weight)
        else:
            return cls.gamma(self.lattice, weight=self.weight)

    def wrap_to_ws(self):
        """Returns a new |Kpoint| in the Wigner-Seitz zone."""
        return self.__class__(wrap_to_ws(self.frac_coords), self.lattice,
                              name=self.name, weight=self.weight)

    def wrap_to_bz(self):
        """Returns a new |Kpoint| in the first unit cell."""
        return self.__class__(wrap_to_bz(self.frac_coords), self.lattice,
                              name=self.name, weight=self.weight)

    def compute_star(self, symmops, wrap_tows=True):
        """Return the star of the kpoint (tuple of |Kpoint| objects)."""
        frac_coords = [self.frac_coords]

        # TODO: This part becomes a bottleneck for large nk!
        for sym in symmops:
            sk_coords = sym.rotate_k(self.frac_coords, wrap_tows=wrap_tows)

            # Add it only if it's not already in the list.
            for prev_coords in frac_coords:
                if issamek(sk_coords, prev_coords): break
            else:
                frac_coords.append(sk_coords)

        return KpointStar(self.lattice, frac_coords, weights=None, names=len(frac_coords) * [self.name])


class KpointList(collections.abc.Sequence):
    """
    Base class defining a sequence of |Kpoint| objects. Essentially consists
    of base methods implementing the sequence protocol and helper functions.
    The subclasses |Kpath| and |IrredZone| provide specialized methods to operate
    on k-points representing a path or list of points in the IBZ, respectively.
    This object is immutable.

    .. Note:

        Algorithms usually need to know what kind of sampling we are using.
        The test can be easily implemented with:

        if kpoints.is_path:
            # code specific to k-paths.

        elif kpoints.is_ibz:
            # code specific to IBZ sampling.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: KpointList
    """
    Error = KpointsError

    @classmethod
    def subclass_from_name(cls, name):
        """Return the class with the given name."""
        if cls.__name__ == name: return cls
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
            reciprocal_lattice: |Lattice| object.
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

    @property
    def reciprocal_lattice(self):
        """|Lattice| object defining the reciprocal lattice."""
        return self._reciprocal_lattice

    def __repr__(self):
        return self.to_string(func=repr)

    def __str__(self):
        return self.to_string(func=str)

    def to_string(self, func=str, title=None, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))
        lines.extend(["%d) %s" % (i, func(kpoint)) for i, kpoint in enumerate(self)])

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
        return not (self == other)

    def index(self, kpoint):
        """
        Returns: the first index of kpoint in self.

        Raises: `ValueError` if not found.
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
            obj: Fractional coordinates or |Kpoint| instance.

        Return:
            (ind, kpoint, dist)

            where ``ind`` is the index in self of the closest k-point.
            ``kpoint`` is the |Kpoint| instance of index ``ind``.
            dist is the distance between ``obj`` and ``kpoint``.
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
        if not self.is_ibz: return (None, None)
        # Test if kptrlatt is diagonal.
        if not is_diagonal(self.ksampling.kptrlatt): return (None, None)
        return self.ksampling.kptrlatt.diagonal(), self.ksampling.shifts

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
        """Fractional coordinates of the k-point as |numpy-array| of shape (len(self), 3)"""
        return self._frac_coords

    def get_cart_coords(self):
        """Cartesian coordinates of the k-point as |numpy-array| of shape (len(self), 3)"""
        return np.reshape([k.cart_coords for k in self], (-1, 3))

    @property
    def names(self):
        """List with the name of the k-points."""
        return [k.name for k in self]

    @property
    def weights(self):
        """|numpy-array| with the weights of the k-points."""
        return np.array([kpoint.weight for kpoint in self])

    def sum_weights(self):
        """Returns the sum of the weights."""
        return np.sum(self.weights)

    def check_weights(self):
        """
        Check if weights are normalized to one.
        Raises: `ValueError` if Weights are not normalized.
        """
        # Weights must be normalized to one.
        wsum = self.sum_weights()
        if abs(wsum - 1) > 1.e-6:
            err_msg = "Kpoint weights should sum up to one while sum_weights is %.3f\n" % wsum
            err_msg += "The list of kpoints does not represent a homogeneous sampling of the BZ\n"
            err_msg += "%s\n%s" % (self.__class__, self.to_string(verbose=0))
            raise ValueError(err_msg)

    def remove_duplicated(self):
        """
        Remove duplicated k-points from self. Returns new :class:`KpointList` instance.
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
        """Returns a |numpy-array| [nkpy, 3] with the fractional coordinates."""
        return np.array(self.frac_coords.copy())

    def to_json(self):
        """
        Returns a JSON_ string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)

    def plot(self, ax=None, **kwargs):
        """Plot k-points with matplotlib."""
        from pymatgen.electronic_structure.plotter import plot_brillouin_zone
        fold = False
        if self.is_path:
            labels = {k.name: k.frac_coords for k in self if k.name}
            frac_coords_lines = [self.frac_coords[line] for line in self.lines]
            return plot_brillouin_zone(self.reciprocal_lattice, lines=frac_coords_lines, labels=labels,
                                       ax=ax, fold=fold, **kwargs)
        else:
            # Not sure this works, I got points outside of the BZ in a simple with Si and Gamma-centered 8x8x8.
            # Don't know if it's a bug in matplotlib or plot_brillouin_zone.
            #print(self.frac_coords)
            return plot_brillouin_zone(self.reciprocal_lattice, kpoints=self.frac_coords,
                                       ax=ax, fold=fold, **kwargs)

    def get_k2kqg_map(self, qpt, atol_kdiff=None):
        """
        Compute mapping k_index --> (k + q)_index, g0

        Args:
            qpt: q-point in fractional coordinate or :class:`Kpoint` instance.
            atol_kdiff: Tolerance used to compare k-points.
                Use _ATOL_KDIFF is atol is None.
        """
        if atol_kdiff is None: atol_kdiff = _ATOL_KDIFF
        if isinstance(qpt, Kpoint):
            qfrac_coords = qpt.frac_coords
        else:
            qfrac_coords = np.reshape(qpt, (3,))

        k2kqg = collections.OrderedDict()
        if np.all(np.abs(qfrac_coords) <= 1e-6):
            # Gamma point, DOH!
            g0 = np.zeros(3, dtype=np.int)
            for ik, _ in enumerate(self):
                k2kqg[ik] = (ik, g0)
        else:
            # N**2 scaling but this algorithm can handle k-paths
            # Note that in principle one could have multiple k+q in k-points
            # but only the first match is considered.
            for ik, kpoint in enumerate(self):
                kpq = kpoint.frac_coords + qfrac_coords
                for ikq, ksearch in enumerate(self):
                    if issamek(kpq, ksearch.frac_coords, atol=atol_kdiff):
                        g0 = np.rint(kpq - ksearch.frac_coords)
                        k2kqg[ik] = (ikq, g0)
                        break

        return k2kqg


class KpointStar(KpointList):
    """
    Star of the kpoint. Note that the first k-point is assumed to be the base
    of the star namely the point that is used to generate the Star.

    .. inheritance-diagram:: KpointStar
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

    .. inheritance-diagram:: Kpath
    """

    @classmethod
    def from_names(cls, structure, knames, line_density=20):
        """
        Generate normalized K-path from list of k-point labels.

        Args:
            structure: |Structure| object.
            knames: List of strings with the k-point labels.
            line_density: Number of points used to sample the smallest segment of the path
        """
        kfrac_coords = structure.get_kcoords_from_names(knames)
        vertices_names = list(zip(kfrac_coords, knames))

        return cls.from_vertices_and_names(structure, vertices_names, line_density=line_density)

    @classmethod
    def from_vertices_and_names(cls, structure, vertices_names, line_density=20):
        """
        Generate normalized K-path from a list of vertices and the corresponding labels.

        Args:
            structure: |Structure| object.
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

        return cls(structure.lattice.reciprocal_lattice, frac_coords=frac_coords,
                   weights=None, names=knames)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0, title=None, **kwargs):
        """
        String representation.

        Args:
            verbose: Verbosity level. Default: 0
        """
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))
        app("K-path contains %s lines. Number of k-points in each line: %s" % (
            len(self.lines), [len(l) for l in self.lines]))
        if verbose:
            for i, line in enumerate(self.lines):
                app("line %d: %s" % (i, line))
        header = "\n".join(lines)

        vids = {line[0] for line in self.lines}

        table = [["Idx", "Frac_coords", "Name", "ds", "Vert",]]
        for i, kpoint in enumerate(self):
            tag = "*" if i in vids else " "
            if verbose == 0 and not tag: continue
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
        |numpy-array| of len(self)-1 elements giving the distance between two
        consecutive k-points, i.e. ds[i] = ||k[i+1] - k[i]|| for i=0,1,...,n-1
        """
        ds = np.zeros(len(self) - 1)
        for i, kpoint in enumerate(self[:-1]):
            ds[i] = (self[i + 1] - kpoint).norm
        return ds

    @lazy_property
    def versors(self):
        """
        Tuple of len(self) - 1 elements with the versors connecting k[i] to k[i+1].
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

    @lazy_property
    def frac_bounds(self):
        """Numpy array of shape [M, 3] with the vertexes of the path in frac coords."""
        frac_bounds = [self[line[0]].frac_coords for line in self.lines]
        frac_bounds.append(self[self.lines[-1][-1]].frac_coords)
        return np.reshape(frac_bounds, (-1, 3))

    @lazy_property
    def cart_bounds(self):
        """Numpy array of shape [M, 3] with the vertexes of the path in frac coords."""
        cart_bounds = [self[line[0]].cart_coords for line in self.lines]
        cart_bounds.append(self[self.lines[-1][-1]].cart_coords)
        return np.reshape(cart_bounds, (-1, 3))

    def find_points_along_path(self, cart_coords, dist_tol=1e-12):
        """
        Find points in ``cart_coords`` lying on the path with distance less than `dist_tol`.
        See find_points_along_path function for API.
        """
        return find_points_along_path(self.cart_bounds, cart_coords, dist_tol=dist_tol)

    def finite_diff(self, values, order=1, acc=4):
        """
        Compute the derivatives of values by finite differences.

        Args:
            values: array-like object with shape=(nkpt) containing the values of the path.
            order: Order of the derivative.
            acc: Accuracy: 4 corresponds to a central difference with 5 points.

        Returns:
            ragged numpy array. The i-th entry is a numpy array with the derivatives on the i-th line.
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

    .. inheritance-diagram:: IrredZone
    """
    #@classmethod
    #def from_ngkpt_or_kppa(cls, structure, ngkpt, shiftk, kptopt=1, verbose=0):
    #    from abipy.tools import duck
    #    if duck.is_listlike(ngkpt):
    #        return cls.from_ngkpt(structure, ngkpt, shiftk, kptopt=kptopt, verbose=verbose)
    #    else:
    #        return cls.from_kppa(structure, kppa, shiftk, kptopt=kptopt, verbose=verbose)

    @classmethod
    def from_ngkpt(cls, structure, ngkpt, shiftk, kptopt=1, verbose=0):
        """
        Build an IrredZone object from (ngkpt, shift) by calling Abinit
        to get the list of irreducible k-points.
        """
        from abipy.abio.factories import gs_input
        from abipy.data.hgh_pseudos import HGH_TABLE
        gsinp = gs_input(structure, HGH_TABLE, spin_mode="unpolarized")
        ibz = gsinp.abiget_ibz(ngkpt=ngkpt, shiftk=shiftk, kptopt=kptopt, verbose=verbose)
        ksampling = KSamplingInfo.from_mpdivs(ngkpt, shiftk, kptopt)

        return cls(structure.reciprocal_lattice, ibz.points, weights=ibz.weights,
                   names=None, ksampling=ksampling)

    @classmethod
    def from_kppa(cls, structure, kppa, shiftk, kptopt=1, verbose=0):
        """
        Build an IrredZone object from (kppa, shift) by calling Abinit
        to get the list of irreducible k-points.
        """
        from abipy.abio.factories import gs_input
        from abipy.data.hgh_pseudos import HGH_TABLE
        gsinp = gs_input(structure, HGH_TABLE, spin_mode="unpolarized", kppa=kppa)
        ibz = gsinp.abiget_ibz(ngkpt=None, shiftk=shiftk, kptopt=kptopt, verbose=verbose)
        ksampling = KSamplingInfo.from_mpdivs(gsinp["ngkpt"], shiftk, kptopt)

        return cls(structure.reciprocal_lattice, ibz.points, weights=ibz.weights,
                   names=None, ksampling=ksampling)

    def __init__(self, reciprocal_lattice, frac_coords, weights=None, names=None, ksampling=None):
        """
        Args:
            reciprocal_lattice: |Lattice| object
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

    def to_string(self, func=str, verbose=0, title=None):
        """String representation."""
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))

        if self.is_mpmesh:
            mpdivs, shifts = self.mpdivs_shifts
            d = "[%d, %d, %d]" % tuple(mpdivs)
            s = ", ".join("[%.1f, %.1f, %.1f]" % tuple(s) for s in shifts)
            app("K-mesh with divisions: %s, shifts: %s" % (d, s))
            app("kptopt: %s (%s)" % (self.ksampling.kptopt, kptopt2str(self.ksampling.kptopt)))
        else:
            app("nkpt: %d" % len(self))
            app(self.ksampling.to_string(verbose=verbose))

        app("Number of points in the IBZ: %s" % len(self))
        for i, k in enumerate(self):
            if i > 10 and verbose == 0:
                app(4 * " " + "... (More than 10 k-points)")
                break
            app("%6d) [%+.3f, %+.3f, %+.3f],  weight=%.3f" % (i, k[0], k[1], k[2], k.weight))

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
    """
    Store metadata defining the k-point sampling according to the abinit conventions.
    One should pass through one of the class methods to create the class, avoid calling __init__ directly.
    """

    KNOWN_KEYS = set([
        "mpdivs",          # Mesh divisions. Defined only if we have a sampling with diagonal kptrlatt else None.
        "kptrlatt",        # [3, 3] matrix defined only if we have a sampling else None.
        "kptrlatt_orig",   # Original set of shifts. Defined only if we have a sampling else None.
        "shifts",          # Actual shifts (Usually one). Defined only if we have a sampling else None.
        "shifts_orig",     # Original shifts specified by the user. Defined only if we have a sampling else None.
        "kptopt",          # Options for k-point generation. Negative if we have a k-path (nbounds - 1).
    ])

    @classmethod
    def as_ksampling(cls, obj):
        """"
        Convert obj into a :class:`KSamplingInfo` instance.
        Accepts: :class:`KSamplingInfo` instance, None (if info are not available) or dict-like object.
        """
        if isinstance(obj, cls): return obj
        if obj is None:
            return cls(mpdivs=None,
                       kptrlatt=None,
                       kptrlatt_orig=None,
                       shifts=None,
                       shifts_orig=None,
                       kptopt=0,
            )

        # Assume dict-like object.
        try:
            return cls(**obj)
        except Exception as exc:
            raise TypeError("Don't know how to convert `%s` into KSamplingInfo object:\n%s" % (type(obj), str(exc)))

    @classmethod
    def from_mpdivs(cls, mpdivs, shifts, kptopt):
        """
        Homogeneous sampling specified in terms of ``mpdivs`` (ngkpt in abinit),
        the set of ``shifts`` and the value of ``kptopt``.
        """
        kptrlatt = kptrlatt_orig = np.diag(mpdivs)
        shifts = shifts_orig = np.reshape(np.array(shifts), (-1, 3))

        return cls(mpdivs=mpdivs, shifts=shifts, shifts_orig=shifts_orig,
                   kptrlatt=kptrlatt, kptrlatt_orig=kptrlatt_orig, kptopt=kptopt)

    @classmethod
    def from_kptrlatt(cls, kptrlatt, shifts, kptopt):
        """
        Homogeneous sampling specified in terms of ``kptrlatt``
        the set of ``shifts`` and the value of ``kptopt``.
        """
        kptrlatt = kptrlatt_orig = np.reshape(kptrlatt, (3, 3))
        shifts = shifts_orig = np.reshape(np.array(shifts), (-1, 3))
        # Test if kptrlatt is diagonal.
        mpdivs = None if not is_diagonal(kptrlatt) else np.diag(kptrlatt)

        return cls(mpdivs=mpdivs, shifts=shifts, shifts_orig=shifts_orig,
                   kptrlatt=kptrlatt, kptrlatt_orig=kptrlatt_orig, kptopt=kptopt)

    @classmethod
    def from_kbounds(cls, kbounds):
        """
        Metadata associated to a k-path specified in terms of boundaries.
        """
        mpdivs, kptrlatt, kptrlatt_orig, shifts, shifts_orig = 5 * (None,)
        kptopt = - (len(np.reshape(kbounds, (-1, 3))) - 1)  # Note -1

        return cls(mpdivs=mpdivs, shifts=shifts, shifts_orig=shifts_orig,
                   kptrlatt=kptrlatt, kptrlatt_orig=kptrlatt_orig, kptopt=kptopt)

    def __init__(self, *args, **kwargs):
        super(KSamplingInfo, self).__init__(*args, **kwargs)
        for k in self:
           if k not in self.KNOWN_KEYS:
               raise ValueError("Unknow key %s" % k)

        # FIXME: monkhorst_pack_folding is not written in e.g. DEN.nc files
        # so we get crazy results because of netCDF4._default_fillvals
        # This part set the value of mpdivs from kptrlatt.
        if self["mpdivs"] is not None and np.any(np.abs(self["mpdivs"]) > 1e+6):
            if self.kptrlatt_orig is not None:
                # We have a sampling
                if np.all(self.kptrlatt_orig == self.kptrlatt) and is_diagonal(self.kptrlatt):
                    self["mpdivs"] = np.diag(self.kptrlatt)
                else:
                    self["mpdivs"] = None
#                    import warnings
#                    warnings.warn("""
#monkhorst_pack_folding variables has not been written to netcdf file.
#Received {mpdivs}
#Setting mpdivs to None, this may create problems in post-processing tools.
#If needed, use python netcdf to change the value of `monkhorst_pack_folding`""".format(mpdivs=self["mpdivs"]))


    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0, title=None, **kwargs):
        """String representation."""
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))

        if self.is_mesh:
            if self.has_diagonal_kptrlatt:
                app("mpdivs: %s with shifts %s and kptopt: %s" % (str(self.mpdivs), str(self.shifts), self.kptopt))
            else:
                app("kptrlatt: %s" % str(self.kptrlatt))
                app("shifts: %s" % str(self.shifts))
                app("kptrlatt_orig: %s" % str(self.kptrlatt_orig))
                app("shifts_orig: %s" % str(self.shifts_orig))
                app("kptopt: %s (%s)" % (str(self.kptopt), kptopt2str(self.kptopt)))

        elif self.is_path:
            app("Path with kptopt: %s" % self.kptopt)
        else:
            app("Neither mesh or path!")

        return "\n".join(lines)

    @property
    def is_mesh(self):
        """True if we have a mesh in the BZ."""
        return self.kptopt > 0 and (self.mpdivs is not None or self.kptrlatt is not None)

    @property
    def is_path(self):
        """True if we have a path in the BZ."""
        return self.kptopt < 0

    #@property
    #def is_homogeneous(self):
    #    """True if we have a homogeneous sampling of the BZ."""
    #    return self.kptopt > 0 and (self.mpdivs is not None or self.kptrlatt is not None)

    @property
    def has_diagonal_kptrlatt(self):
        """True if sampling with diagonal kptrlatt."""
        if self.kptrlatt is None: return False
        return is_diagonal(self.kptrlatt)


class KpointsReaderMixin(object):
    """
    Mixin class that provides methods for reading k-point data from a netcdf
    file written according to the ETSF-IO specifications.
    """
    def read_kpoints(self):
        """
        Factory function: returns an instance of :class:`Kpath` or :class:`IrredZone`
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

        #raise ValueError("Only homogeneous samplings or paths are supported!\n"
        #                 "ksampling info:\n%s" % str(ksampling))

    def read_ksampling_info(self):
        """
        Read information on the k-point sampling. Return :class:`KSamplingInfo` object.
        """
        # FIXME: in v8.0, the SIGRES files does not have kptopt, kptrlatt_orig and shiftk_orig
        kptrlatt = self.read_kptrlatt()
        shifts = self.read_kshifts()

        return KSamplingInfo(
            mpdivs=self.read_kmpdivs(),
            kptrlatt=kptrlatt,
            kptrlatt_orig=self.read_value("kptrlatt_orig", default=kptrlatt),
            shifts=shifts,
            shifts_orig=self.read_value("shiftk_orig", default=shifts),
            kptopt=int(self.read_value("kptopt", default=0)),
        )

    def read_kfrac_coords(self):
        """Fractional coordinates of the k-points"""
        return self.read_value("reduced_coordinates_of_kpoints")

    def read_kweights(self):
        """Returns the weight of the k-points. None if not found."""
        return self.read_value("kpoint_weights")

    def read_kshifts(self):
        """Returns the shifts of the k-mesh in reduced coordinates. None if not found."""
        try:
            return self.read_value("shiftk")
        except self.Error:
            return self.read_value("kpoint_grid_shift")

    def read_kmpdivs(self):
        """Returns the Monkhorst-Pack divisions defining the mesh. None if not found."""
        #return self.read_value("monkhorst_pack_folding")
        return self.none_if_masked_array(self.read_value("monkhorst_pack_folding"))

    def read_kptrlatt(self):
        """Returns ABINIT variable kptrlatt. None if not found."""
        try:
            return self.read_value("kptrlatt")
        except self.Error:
            return self.read_value("kpoint_grid_vectors")


class KpointsReader(ETSF_Reader, KpointsReaderMixin):
    """
    This object reads k-point data from a netcdf file.

    .. inheritance-diagram:: KpointsReader
    """


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

    def to_string(self, verbose=0, title=None, **kwargs):
        """String representation"""
        lines = collections.deque(); app = lines.append
        if title is not None: app(marquee(title, mark="="))

        app("mesh %s, shift %s, time-reversal: %s, Irred points: %d" % (
            self.mesh, self.kshift, self.has_timrev, self.nibz))

        app("Irreducible k-points with number of points in star:\n")
        for ik, kpt in enumerate(self.ibz):
            app("%s: [%9.6f, %9.6f, %9.6f], nstar: %d" % (ik, kpt[0], kpt[1], kpt[2], self.weights[ik] * self.nbz))

        return "\n".join(lines)

    def print_bz2ibz(self, file=sys.stdout):
        """Print BZ --> IBZ mapping."""
        print("BZ points --> IBZ points mapping", file=file)
        for ik_bz, ik_ibz in enumerate(self.bz2ibz):
            print("%6d) [%9.6f, %9.6f, %9.6f], ===> %6d) [%9.6f, %9.6f, %9.6f]," %
                (ik_bz, self.bz[ik_ibz][0], self.bz[ik_ibz][1], self.bz[ik_ibz][2],
                ik_ibz, self.ibz[ik_ibz][0], self.ibz[ik_ibz][1], self.ibz[ik_ibz][2]), file=file)


def dist_point_from_line(x0, x1, x2):
    """
    Return distance from point x0 to line x1 - x2. Cartesian coordinates are used.
    See <http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html>
    """
    denom = x2 - x1
    denomabs = np.sqrt(np.dot(denom, denom))
    numer = np.cross(x0 - x1, x0 - x2)
    numerabs = np.sqrt(np.dot(numer, numer))
    return numerabs / denomabs


def find_points_along_path(cart_bounds, cart_coords, dist_tol):
    """
    Find points in ``cart_coords`` lying on the path defined by ``cart_bounds``.

    Args:
        cart_bounds: [N, 3] array with the boundaries of the path in Cartesian coordinates.
        cart_coords: [M, 3] array with the points in Cartesian coordinate
        dist_tol: A point is considered to be on the path if its distance from the line
            is less than dist_tol.

    Return: namedtuple with the following attributes:

        (ikfound, dist_list, path_ticks)

        ikfound is a numpy array with the indices of the points lying on the path. Empty if no point is found.
        dist_list: numpy array with the distance of the points along the line.
        path_ticks: numpy array with the ticks associated to the input k-path.
    """
    ikfound, dist_list, path_ticks = [], [], [0]

    dl = 0  # cumulative length of the path
    for ibound, x0 in enumerate(cart_bounds[:-1]):
        x1 = cart_bounds[ibound + 1]
        B = x0 - x1
        #B = x1 - x0
        dk = np.sqrt(np.dot(B,B))
        #print("x0", x0, "x1", x1)
        path_ticks.append(path_ticks[ibound] + dk)
        for ik, k in enumerate(cart_coords):
            dist = dist_point_from_line(k, x0, x1)
            #print(frac_coords[ik], dist)
            if dist > dist_tol: continue
            # k-point is on the cart_bounds
            A = x0 - k
            #A = k - x0
            x = np.dot(A, B)/dk
            #print("k-x0", A, "B", B)
            #print(frac_coords[ik], x, x > 0 and x < dist_tol + dk)
            if dist_tol + dk >= x >= 0:
                # k-point is within the cart_bounds range
                # append k-point coordinate along the cart_bounds while avoing duplicate entries.
                if ikfound and ik == ikfound[-1]: continue
                ikfound.append(ik)
                dist_list.append(x + dl)

        dl = dl + dk

    # Sort dist_list and ikfound by cumulative length while removing possible duplicated entries.
    dist_list, isort = np.unique(np.array(dist_list).round(decimals=5), return_index=True)

    return dict2namedtuple(ikfound=np.array(ikfound)[isort],
                           dist_list=dist_list,
                           path_ticks=np.array(path_ticks))
