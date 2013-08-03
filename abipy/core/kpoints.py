"""This module defines objects describing the sampling of the Brillouin Zone."""
from __future__ import division, print_function

import os
import collections
import numpy as np

from numpy import allclose,  where, around, asarray 
from abipy.core.exceptions import AbipyException
from abipy.iotools import as_etsfreader
from abipy.tools.derivatives import finite_diff

__all__ = [
    "issamek",
    "wrap_to_ws",
    "wrap_to_bz",
    "askpoints",
    "Kpoint",
    "Kpath",
    "IrredZone",
    "KpointsInfo",
    "kpoints_factory",
]

# Tolerance used to compare k-points.
_ATOL_KDIFF = 1e-8

def isinteger(x, atol=_ATOL_KDIFF):
    """
    True if all x is integer within the absolute tolerance atol.

    >>> isinteger([1., 2.])
    True
    >>> isinteger(1.01, atol=0.011)
    True
    >>> isinteger([1.01, 2])
    False
    """
    int_x = around(x)
    return allclose(int_x, x, atol=atol)


def issamek(k1, k2, atol=1e-08):
    """
    True if k1 and k2 are equal modulo a lattice vector.

    Examples

    >>> issamek([1,1,1], [0,0,0])
    True
    >>> issamek([1.1,1,1], [0,0,0], atol=0.1)
    True
    >>> issamek(0.00003, 1)
    False
    """
    kdiff = asarray(k1) - asarray(k2)
    return isinteger(kdiff, atol=atol)


def wrap_to_ws(x):
    """
    Transforms x in its corresponding reduced number in the interval ]-1/2,1/2].
    """
    w = x % 1
    return where(w > 0.5, w-1.0, w)


def wrap_to_bz(x):
    """
    Transforms x in its corresponding reduced number in the interval [0,1[."
    """
    return x % 1


class KpointsError(AbipyException):
    """Base error class for KpointList exceptions."""


def askpoints(obj, lattice, weigths=None, names=None):
    """
    Convert obj into a list of k-points.

    Args:
        obj:
            Kpoint or list of Kpoint objects or array-like object.
        lattice:
            Reciprocal lattice.
        weights:
            k-point weights. Ignored if obj is already a `Kpoint` instance or a list
            of `Kpoint` items.
        name:
            string with the name of the k-point. Ignored if obj is already a `Kpoint`
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
        return [Kpoint(obj, lattice, weight=weigths, name=names)]
    elif ndim == 2:
        nk = len(obj)
        if weigths is None: weigths = nk * [None]
        if names is None: names = nk * [None]
        return [Kpoint(rc, lattice, weight=w, name=l) for (rc, w, l) in zip(obj, weigths, names)]
    else:
        raise ValueError("ndim > 2 is not supported")

##########################################################################################


class Kpoint(object):
    """Class defining one k-point."""
    __slots__ = [
        "_frac_coords",
        "_lattice",
        "_weight",
        "_name",
    ]

    # Tolerance used to compare k-points.
    @property
    def ATOL_KDIFF(self):
        return _ATOL_KDIFF

    def __init__(self, frac_coords, lattice, weight=None, name=None):
        """
        Args:
            frac_coords:
                Reduced coordinates.
            lattice:
                `Lattice` object describing the reciprocal lattice.
            weights: 
                k-point weight (optional, set to zero if not given).
            name:
                string with the name of the k-point (optional)
        """
        self._frac_coords = np.asarray(frac_coords)
        assert len(self.frac_coords) == 3

        self._lattice = lattice
        self.set_weight(weight)
        self.set_name(name)

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
        """Cartesian coordinates of the k-points."""
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

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        s = "Kpoint: [%.3f, %.3f, %.3f]" % tuple(self.frac_coords)
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
    def askpoint(cls, obj, lattice):
        """
        Convert obj into a Kpoint instance.

        Args:
            obj:
                `Kpoint` instance or array-like with the reduced coordinates
            lattice:
                `Lattice` object defining the reciprocal lattice.
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
        cart = self.lattice.get_cartesian_coords(self.frac_coords)
        return np.sqrt(np.dot(cart, cart))

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

##########################################################################################

def kpoints_factory(filepath):
    """
    Factory function: returns an instance of [Kpath, IrredZone]
    from a netcdf file written according to the ETSF-IO specifications.
    """
    file, closeit = as_etsfreader(filepath)
    structure = file.read_structure()

    kinfo = KpointsInfo.from_file(filepath)

    if kinfo.is_sampling:
        obj = IrredZone(structure.reciprocal_lattice, kinfo.frac_coords, kinfo)

    elif kinfo.is_path:
        obj = Kpath(structure.reciprocal_lattice, kinfo.frac_coords, kinfo)

    else:
        raise ValueError("Only path or sampling modes are supported!")

    if closeit:
        file.close()

    return obj

#def qpoints_factory(filepath):


class KpointList(collections.Sequence):
    """
    Base class defining a sequence of `Kpoint` objects. Essentially consists 
    of base methods implementing the sequence protocol and helper functions.
    """
    Error = KpointsError

    def __init__(self, reciprocal_lattice, frac_coords, weights=None, names=None):
        """
        Args:
            reciprocal_lattice:
                `Lattice` object.
            frac_coords:
                Array-like object with the reduced coordinates of the k-points.
            weights:
                List of k-point weights.
            names:
                List of k-point names.
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
        for (i, rcs) in enumerate(frac_coords):
            name = None if names is None else names[i]
            self._points.append(Kpoint(rcs, self.reciprocal_lattice, weight=weights[i], name=name))

    @property
    def reciprocal_lattice(self):
        return self._reciprocal_lattice

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

    #def __add__(self, other):
    #    assert self.reciprocal_lattice == other.reciprocal_lattice
    #    return KpointList(self.reciprocal_lattice, 
    #        frac_coords = [k.frac_coords for k in self] + [k.frac_coords for k in other], 
    #        weights = None,
    #        names = [k.name for k in self] + [k.name for k in other]
    #        )

    def index(self, kpoint):
        """
        Returns first index of kpoint. Raises ValueError if not found.
        """
        return self._points.index(kpoint)

    def find(self, kpoint):
        """
        Returns first index of kpoint. -1 if not found
        """
        try:
            return self.index(kpoint)
        except ValueError:
            return -1

    def count(self, kpoint):
        """Return number of occurrences of kpoint"""
        return self._points.count(kpoint)

    @property
    def frac_coords(self):
        return self._frac_coords

    @property
    def weights(self):
        return np.array([kpoint.weight for kpoint in self])

    def sum_weights(self):
        return np.sum(self.weights)

    def to_fortran_arrays(self):
        fort_arrays = collections.namedtuple("FortranKpointListArrays", "frac_coords")

        return fort_arrays(
            frac_coords=np.asfortranarray(self.frac_coords.T),
        )


class KpointStar(KpointList):
    """
    Start of the kpoint. Note that the first k-point is assumed to be the base 
    of the star namely the point that is used to generate the Star.
    """
    @property
    def base_point(self):
        return self[0]

    @property
    def name(self):
        return self.base_point.name


class Kpath(KpointList):
    """This object describes a path in reciprocal space."""

    def __init__(self, reciprocal_lattice, frac_coords, kinfo):
        #names = kinfo.pop("names", None)
        names = None

        super(Kpath, self).__init__(reciprocal_lattice, frac_coords, weights=kinfo.weights, names=names)
        # time-reversal?
        #bounds = kinfo.pop("bounds", None)
        #if bounds is None:
        #    self.bounds = np.reshape(bounds, (-1,3))
        #else:
        #    pass

    @property
    def ds(self):
        """
        ndarray of len(self)-1 elements giving the distance between two
        consecutive k-points, i.e. ds[i] = ||k[i+1] - k[i]||.
        """
        try:
            return self._ds
        except AttributeError:

            self._ds = ds = np.zeros(len(self) - 1)
            for (i, kpoint) in enumerate(self[:-1]):
                ds[i] = (self[i + 1] - kpoint).norm

            return self._ds

    @property
    def versors(self):
        """Tuple of len(self)-1 elements with the versors connecting k[i] to k[i+1]."""
        try:
            return self._versors

        except AttributeError:
            versors = (len(self) - 1) * [None, ]
            versors[0] = Kpoint.gamma(self.reciprocal_lattice)

            for (i, kpt) in enumerate(self[:-1]):
                versors[i] = (self[i + 1] - kpt).versor()
            self._versors = tuple(versors)

            return self._versors

    @property
    def num_lines(self):
        """The number of lines forming the path."""
        return len(self.lines)

    @property
    def lines(self):
        """
        tuple with the list of indices of the points belonging to the same line.
        """
        try:
            return self._lines

        except AttributeError:
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
            self._lines = tuple(lines)

            return self._lines

    def finite_diff(self, values, order=1, acc=4):
        """
        Compute the derivatives of values by finite differences.

        Args:
            values:
                array-like object with the values of the path.
            order:
                Order of the derivative.
            acc:
                Accuracy: 4 corresponds to a central difference with 5 points.

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
    An IrredZone is a (immutable) sequence of points in the irreducible wedge of the BZ.
    Each point has a weight whose sum must equal 1 so that we can integrate quantities 
    in the full Brillouin zone.
    """

    def __init__(self, reciprocal_lattice, frac_coords, kinfo):
        """
        Args:
            reciprocal_lattice:
                `Lattice`
            frac_coords:
                Array-like object with the reduced coordinates of the points.
            kinfo:
                `Kinfo` object, essentially consists of a dictionary with 
                the parameters used to generate the mesh in the full Brillouin zone. 
        """
        names = None
        #names = kinfo.pop("names", None)
        super(IrredZone, self).__init__(reciprocal_lattice, frac_coords, weights=kinfo.weights, names=names)

        # Weights must be normalized to one.
        wsum = self.sum_weights()
        if abs(wsum - 1) > 1.e-6:
            raise ValueError("K-point weights must be normalized to one but wsum = %s" % wsum)

        # time-reversal?
        self.kptopt = kinfo.kptopt
        self.shifts = kinfo.shifts

        self.mpdivs = None
        self.kptrlatt = None

        if kinfo.mpdivs is not None:
            # MP folding
            self.mpdivs = kinfo.mpdivs

        elif kinfo.kptrlatt is not None:
            # kptrlatt case.
            self.kptrlatt = kinfo.kptrlatt
            # Diagonal kptrlatt is equivalent to MP folding.
            #if self.kptrlatt  and ...
            #self.mpdivs =

        else:
            raise ValueError("Either MP mesh info or kptrlatt must be present in kinfo")

    @property
    def has_mpmesh(self):
        "True if the mesh has been defined in terms of Monkhors-Pack divisions."""
        return self.mpdivs is not None

    #@classmethod
    #def from_info(cls, reciprocal_lattice, info):
    #    """Initialize the IrredZone from reciprocal_lattice and info."""
    #    raise NotImplementedError("")


class KpointsInfo(dict):

    class MANDATORY(object):
        """Mandatory keys."""

    class OPTIONAL(object):
        """Optional keys."""

    KNOWN_KEYS = {
        "reduced_coordinates_of_kpoints": MANDATORY,
        "kpoint_weights": OPTIONAL,
        "kpoint_grid shift": OPTIONAL,
        "monkhorst_pack_folding": OPTIONAL,
        "kpoint_grid vectors": OPTIONAL,
        "kptopt": OPTIONAL,
    }

    def __init__(self, *args, **kwargs):
        super(KpointsInfo, self).__init__(*args, **kwargs)

        for k in self:
            if k not in self.KNOWN_KEYS:
                raise ValueError("Unknow key %s" % k)

        for k, v in self.KNOWN_KEYS.items():
            if v is self.MANDATORY and k not in self:
                raise ValueError("Mandatory key %s is missing" % k)

    @classmethod
    def from_file(cls, filepath):
        """
        Initialize the object from a Netcdf file with data
        saved in the ETSF-IO format.
        """
        file, closeit = as_etsfreader(filepath)

        d = {}
        for k in cls.KNOWN_KEYS:
            try:
                d[k] = file.read_value(k)
            except file.Error:
                pass

        if closeit:
            file.close()

        return cls(**d)

    @property
    def is_sampling(self):
        """True if we have a homogeneous sampling of the BZ."""
        return ("monkhorst_pack_folding" in self or
                "kpoint_grid vectors" in self)

    @property
    def is_path(self):
        """True if we have a path in the BZ."""
        return not self.is_sampling

    @property
    def frac_coords(self):
        """Fractional coordinates of the k-points"""
        return np.reshape(self.get("reduced_coordinates_of_kpoints"), (-1, 3))

    @property
    def weights(self):
        return self.get("kpoint_weights", None)

    @property
    def shifts(self):
        shifts = self.get("kpoint_grid_shift", np.zeros(3))
        return np.reshape(shifts, (-1, 3))

    @property
    def mpdivs(self):
        return self.get("monkhorst_pack_folding", None)

    @property
    def kptrlatt(self):
        return self.get("kpoint_grid vectors", None)

    @property
    def kptopt(self):
        return self.get("kptopt", None)

##########################################################################################


class Kmesh(object):
    """
    This object describes the sampling of the full Brillouin zone.

    A Kmesh object has a set of irreducible points, information on how
    to generate the mesh in the full BZ. It also provides methods
    to symmetrized and visualize quantities in k-space.
    """
    def __init__(self, structure, mpdivs, shifts, ibz):
        """
        Args:
            structure:
                `Structure` instance.
            mpdivs:
                Number of Monkhorst-Pack divisions along the reduced directions kx, ky, kz.
            shifts:
                Shifts of the mesh.
            ibz:
                k-points in the irreducible wedge.

        Provides methods to symmetrize k-dependent quantities with the full
        symmetry of the structure. e.g. bands, occupation factors, phonon frequencies.
        """
        self._structure = structure
        assert stucture.reciprocal_lattice == ibz.reciprocal_lattice

        self.mpdivs = np.asarray(mpdivs)
        assert self.mpdivs.shape == (3,)
        self.nx, self.ny, self.nz = self.mpdivs

        self._shifts = np.reshape(shifts, (-1, 3))

        self._ibz = ibz

        #grids_1d = 3 * [None]
        #for i in range(3):
        #    grids_1d[i] = np.arange(0, self.mpdivs[i])
        #self.grids_1d = tuple(grids_1d)

    @property
    def structure(self):
        """Crystalline Structure."""
        return self._structure

    @property
    def ibz(self):
        """`IrredZone` object."""
        return self._ibz

    @property
    def len_ibz(self):
        """Number of points in the IBZ."""
        return len(self.ibz)

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

    #@property
    #def tables_for_shift(self):
        #try:
        #    return self._tables_for_shift
        #except AttributeError:
        #   Build the mapping BZ --> IBZ for the different grids.
        #   map[i,j,k] --> index of the point in the IBZ
        #   self._tables_for_shift = map_mesh2ibz(self.structure, self.mpdivs, self.shifts, self.ibz)
        #   return self._tables_for_shift

    #@property
    #def flat_tables(self):

    def symmetrize_data_ibz(self, data_ibz): #, pbc=3*[False], korder="c"):
        assert len(data_ibz) == self.len_ibz 
        data_bz = np.empty(self.len_bz, dtype=data_ibz.dtype)

        bz2ibz = self.bz2ibz
        for ik_bz in range(self.len_bz):
            data_bz[ik_bz] = data_ibz[bz2ibz[ik_bz]]

        return data_bz

    def plane_cut(self, values_ibz):
        """
        Symmetrize values in the IBZ to have them on the full BZ, then
        select a slice along the specified plane E.g. plane = (1,1,0).
        """
        assert len(values_ibz) == len(self.ibz)
        #indices =
        z0 = 0
        plane = np.empty((self.nx, self.ny))

        kx, ky = range(self.nx), range(self.ny)
        for x in kx:
            for y in ky:
                ibz_idx = self.map_xyz2ibz[x, y, z0]
                plane[x, y] = values_ibz[ibz_idx]

        kx, ky = np.meshgrid(kx, ky)
        return kx, ky, plane

