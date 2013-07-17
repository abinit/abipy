"""This module defines objects describing the sampling of the Brillouin Zone."""
from __future__ import division, print_function

import os
import collections
import numpy as np

from abipy.core.exceptions import AbipyException
from abipy.iotools import as_etsfreader
from abipy.kpoints.utils import wrap_to_ws, wrap_to_bz, issamek
from abipy.tools.derivatives import finite_diff

__all__ = [
    "askpoints",
    "Kpoint",
    "Kpath",
    "IrredZone",
    "KpointsInfo",
    "kpoints_factory",
]


def askpoints(obj, lattice, weigths=None, labels=None):
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
        label:
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
        return [Kpoint(obj, lattice, weight=weigths, label=labels)]
    elif ndim == 2:
        nk = len(obj)
        if weigths is None: weigths = nk * [None]
        if labels is None: labels = nk * [None]
        return [Kpoint(rc, lattice, weight=w, label=l) for (rc, w, l) in zip(obj, weigths, labels)]
    else:
        raise ValueError("ndim > 2 is not supported")

##########################################################################################


class Kpoint(object):
    """Class defining one k-point."""
    __slots__ = [
       "_frac_coords",
       "_lattice",
       "_weight",
       "_label",
    ]
    # Tolerance used to compare k-points.
    _ATOL_KDIFF = 1e-08

    def __init__(self, frac_coords, lattice, weight=None, label=None):
        """
        Args:
            frac_coords:
                Reduced coordinates.
            lattice:
                Reciprocal lattice object.
            weights: 
                k-point weight (optional, set to zero if not given).
            label:
                string with the name of the k-point (optional)
        """
        self._frac_coords = np.asarray(frac_coords)
        assert len(self.frac_coords) == 3

        self._lattice = lattice
        self._weight = weight
        self._label = label

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
    def label(self):
        """Label of the k-point. None if not available."""
        return self._label

    def set_label(self, label):
        """Set the label of the k-point."""
        self._label = label

    def __str__(self):
        s = "(%.3f, %.3f, %.3f)" % tuple(self.frac_coords)
        if self._weight is not None:
            s += ", weight = %f" % self.weight
        return s

    # Kpoint algebra.
    def __add__(self, other):
        return self.__class__(self.frac_coords + other.frac_coords, self.lattice)

    def __sub__(self, other):
        return self.__class__(self.frac_coords - other.frac_coords, self.lattice)

    def __eq__(self, other):
        # Comparison between two Kpoints
        try:
        #if isinstance(other, Kpoint):
            #return issamek(self.frac_coords, other.frac_coords, atol=self._ATOL_KDIFF)
            kdiff = self.frac_coords - other.frac_coords
            return all(abs(np.around(kdiff) - kdiff) < self._ATOL_KDIFF)

        except AttributeError:
            # Kpoint vs iterable (e.g. list)
            return issamek(self.frac_coords, other, atol=self._ATOL_KDIFF)

    def __ne__(self, other):
        return not self == other

    def __getitem__(self, slice):
        return self.frac_coords[slice]

    #def __hash__(self):
    #    return tuple(self.frac_coords).__hash__()

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
            return cls(obj, lattice, weight=None, label=None)

    @classmethod
    def gamma(cls, lattice, weight=None):
        """Constructor for the Gamma point."""
        return cls(np.zeros(3), lattice, weight=weight, label="$\Gamma$")

    def copy(self):
        """Deep copy."""
        return self.__class__(self.frac_coords.copy(), self.lattice.copy(),
                              weight=self.weight, label=self.label)

    @property
    def norm(self):
        """Norm of the kpoint."""
        cart = self.lattice.get_cartesian_coords(self.frac_coords)
        return np.sqrt(np.dot(cart, cart))

    def versor(self):
        """Returns the versor i.e. ||k|| = 1"""
        cls = self.__class__
        try:
            return cls(self.frac_coords/self.norm, self.lattice, weight=self.weight)
        except ZeroDivisionError:
            return cls.gamma(self.lattice, weight=self.weight)

    def wraptows(self):
        """Returns a new kpoint in the Wigner-Seitz zone."""
        return self.__class__(wrap_to_ws(self.frac_coords), self.lattice, weight=self.weight)

    def wrapto1bz(self):
        """Returns a new kpoint in the first unit cell."""
        return self.__class__(wrap_to_bz(self.frac_coords), self.lattice, weight=self.weight)

    def get_star(self, symmops, wrap_tows=True):
        """Return the star of the kpoint (tuple of Kpoint object)."""
        star = []
        for sym in symmops:
            sk_coords = sym.rotate_k(self.frac_coords, wrap_tows=wrap_tows)
            star.append(Kpoint(sk_coords, self.lattice))
        return tuple(star)

##########################################################################################


class KpointsError(AbipyException):
    """Base error class for Kpoints exceptions."""


class KpointNotFoundError(KpointsError):
    """Raised when the k-point cannot be found in Kpoints."""


class Kpoints(collections.Sequence):
    """
    Base class defining a sequence of `Kpoint` objects. Essentially consists 
    of base methods implementing the sequence protocol and helper functions.
    """
    Error = KpointsError

    def __init__(self, structure, frac_coords, weights=None, labels=None):
        """
        Args:
            structure:
                `Structure` object.
            frac_coords:
                Array-like object with the reduced coordinates of the k-points.
            weights:
                List of k-point weights.
            labels:
                List of k-point labels.
        """
        self.structure = structure

        frac_coords = np.reshape(frac_coords, (-1,3))

        if weights is not None:
            assert len(weights) == len(frac_coords)

        if labels is not None:
            assert len(labels) == len(frac_coords)

        self._points = []
        for (i, rcs) in enumerate(frac_coords):
            weight = None if weights is None else weights[i]
            label = None if labels is None else labels[i]
            kpt = Kpoint(rcs, self.reciprocal_lattice, weight=weight, label=label)
            self._points.append(kpt)

    @property
    def reciprocal_lattice(self):
        return self.structure.reciprocal_lattice

    # Sequence protocol.
    def __str__(self):
        lines = ["%d) %s" % ik for ik in enumerate(self)]
        return "\n".join(lines)

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

    def index(self, kpoint):
        """
        Returns first index of kpoint.
        Raises ValueError if the value is not present..
        """
        try:
            return self._points.index(kpoint)
        except ValueError:
            raise KpointNotFoundError

    def count(self, kpoint):
        """Return number of occurrences of kpoint"""
        return self._points.count(kpoint)

    @property
    def sum_weights(self):
        """Return the sum of the weights of the k-points."""
        try:
            return self._sum_weights
        except AttributeError:
            self._sum_weights = np.sum(k.weight for k in self)
            return self._sum_weights

    def asarray(self):
        """Returns a ndarray with the fractional coordinates of the k-points."""
        coords = np.empty((len(self),3))
        for i, k in enumerate(self):
            coords[i] = k.frac_coords
        return coords

##########################################################################################


class Kpath(Kpoints):
    """This object describes a path in reciprocal space."""

    def __init__(self, structure, frac_coords, kinfo):
        #labels = kinfo.pop("labels", None)
        labels = None

        super(Kpath, self).__init__(structure, frac_coords,
                                    weights=kinfo.weights, labels=labels)

        # time-reversal?
        #bounds = kinfo.pop("bounds", None)
        #if bounds is None:
        #    self.bounds = np.reshape(bounds, (-1,3))
        #else:
        #    pass

    #@classmethod
    #def from_bounds(cls, structure, bounds, ndivsm)

    #@classmethod
    #def automatic(cls, structure, ndivsm)

    @property
    def ds(self):
        """
        ndarray of len(self)-1 elements giving the distance between two
        consecutive k-points, i.e. ds[i] = ||k[i+1] - k[i]||.
        """
        try:
            return self._ds
        except AttributeError:
            self._ds = ds = np.zeros(len(self)-1)
            for (i, kpoint) in enumerate(self[:-1]):
                ds[i] = (self[i+1] - kpoint).norm
            return self._ds

    @property
    def versors(self):
        """Tuple of len(self)-1 elements with the versors connecting k[i] to k[i+1]."""
        try:
            return self._versors
        except AttributeError:
            versors = (len(self) - 1) * [None,]
            versors[0] = Kpoint.gamma(self.reciprocal_lattice)

            for (i, kpt) in enumerate(self[:-1]):
                versors[i] = (self[i+1] - kpt).versor()
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
            for i, v in enumerate(self.versors[1:]):
                i += 1
                if v != prev:
                    prev = v
                    lines.append(indices + [i])
                    indices = [i]
                else:
                    indices += [i]
            lines.append(indices + [len(self)-1])

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

##########################################################################################


class IrredZone(Kpoints):
    """
    An IrredZone is a (immutable) sequence of points in the irreducible wedge of the BZ.
    Each point has a weight whose sum must equal 1 so that we can integrate quantities 
    in the full Brillouin zone.
    """
    def __init__(self, structure, frac_coords, kinfo):
        """
        Args:
            structure:
                pymatgen `Structure`
            frac_coords:
                Array-like object with the reduced coordinates of the points.
            kinfo:
                `Kinfo` object, essentially consists of a dictionary with 
                the parameters used to generate the mesh in the full Brillouin zone. 
        """
        labels = None
        #labels = kinfo.pop("labels", None)
        super(IrredZone, self).__init__(structure, frac_coords,
                                     weights=kinfo.weights, labels=labels)

        # Weights must be normalized to one.
        wsum = sum(kpt.weight for kpt in self)
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
    #def from_info(cls, structure, info):
    #    """Initialize the IrredZone from structure and info."""
    #    raise NotImplementedError("")

##########################################################################################


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
    def from_file(cls, file):
        """
        Initialize the object from a Netcdf file with data
        saved in the ETSF-IO format.
        """
        file, closeit = as_etsfreader(file)

        d = {}
        for k in cls.KNOWN_KEYS:
            try:
                d[k] = file.read_value(k)
            except file.Error:
                pass

        if closeit: file.close()
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
        return np.reshape(shifts, (-1,3))

    @property
    def mpdivs(self):
        return self.get("monkhorst_pack_folding", None)

    @property
    def kptrlatt(self):
        return self.get("kpoint_grid vectors", None)

    @property
    def kptopt(self):
        return self.get("kptopt", None)


def kpoints_factory(file):
    """
    Factory function: returns an instance of [Kpath, IrredZone]
    from a netcdf file written according to the ETSF-IO specifications.
    """
    file, closeit = as_etsfreader(file)
    structure = file.read_structure()

    kinfo = KpointsInfo.from_file(file)

    if kinfo.is_sampling:
        obj = IrredZone(structure, kinfo.frac_coords, kinfo)

    elif kinfo.is_path:
        obj = Kpath(structure, kinfo.frac_coords, kinfo)

    else:
        raise ValueError("Only path or sampling modes are supported!")

    if closeit:
        file.close()

    return obj

##########################################################################################


def map_mesh2ibz(structure, mpdivs, shifts, ibz):
    """
    This function computes the mapping between the 
    points in the full BZ and the points in the IBZ.

    Args:
        structure:
            pymatgen `Structure` instance.
        mpdivs
            The three MP divisions
        shifts:
            Array-like object with the MP shift.
        ibz:
            array-like object with the reduced coordinates of the
            points in the IBZ.

    Returns
        ndarray array that maps the full mesh onto ibz i.e. 

            bzmap[idx_bz] = id_ibz

        where idx_bz is the index the BZ and idx_ibz the index in the IBZ.
    """
    # Build bz grid.
    mpdivs = np.asarray(mpdivs)
    shifts = np.asarray(shifts)

    if len(shifts) != 1:
        raise ValueError("Too many shifts!")
    shift = shifts[0]

    bz = np.empty((mpdivs.prod(),3))

    count = 0
    for i in range(mpdivs[0]):
        x = (i + shift[0]) / mpdivs[0]
        for j in range(mpdivs[1]):
            y = (j + shift[1]) / mpdivs[1]
            for k in range(mpdivs[2]):
                z = (k + shift[2]) / mpdivs[2]
                bz[count] = (x,y,z)
                count += 1

    # Build mapping.
    kmap = map_bz2ibz(structure, bz, ibz)
    return np.reshape(kmap, mpdivs)


def map_bz2ibz(structure, bz, ibz):
    """
    This functon computes  the mapping between the Brillouin zone and the Irreducible wedge.
    Args:
        structure:
            pymatgen `Structure` instance.
        bz:
            Reduced coordinates of the points in the BZ
        ibz:
            Reduced coordinates of the points in the IBZ
    """
    #from pymatgen.core.finder import SymmetryFinder
    #finder = SymmetryFinder(structure)

    #Returns:
    #            Numbering of reducible kpoints. Equivalent kpoints will have the
    #            same number. The number of unique values is the number of the
    #            irreducible kpoints.

    #finder.get_ir_kpoints_mapping(kpoints, is_time_reversal=True)

    #Returns:
    #        A list of irreducible kpoints and their weights as a list of
    #        tuples [(ir_kpoint, weight)], with ir_kpoint given
    #        in fractional coordinates

    #finder.get_ir_reciprocal_mesh(mesh=(10, 10, 10), shift=(0, 0, 0), is_time_reversal=True)

    ibz = askpoints(ibz, structure.reciprocal_lattice)
    bz = askpoints(bz, structure.reciprocal_lattice)
    bz2ibz = np.empty(len(bz), np.int)
    bz2ibz.fill(-1)

    # Only Ferromagnetic symmetries are used.
    fm_symmops = structure.fm_symops

    #from .shells import Shells
    #func = lambda k: k.wraptows().norm
    #bz_shells = Shells(bz, func=func)
    #ibz_shells = Shells(ibz, func=func)

    #miss = []
    #for fsh in bz_shells:
    #    try:
    #        ish = ibz_shells.get_from_value(fsh.value)
    #    except ValueError:
    #        miss.append(fsh)
    #        continue

    #    for bz_idx, kbz in fsh.indexitem():
    #        for ibz_idx, kibz in ish.indexitem():
    #            if kbz in kibz.get_star(fm_symmops, wrap_tows=False):
    #                bz2ibz[bz_idx] = ibz_idx
    #                break
    #        #else:
    #            #raise ValueError("Not found")
    #if miss:
    #    print("len(miss) = ",len(miss))
    #    for fsh in miss:
    #        for bz_idx, kfull in fsh.indexitem():
    #            for ibz_idx, kirred in enumerate(ibz):
    #                if kfull in kirred.iter_star(fm_symmops):
    #                    bz2ibz[bz_idx] = ibz_idx
    #                    #break
    #            #else:
    #                #raise ValueError("Full k-point: %s not found" % kfull)
    #if all(bzibz != -1):
    #    return bz2ibz
    ATOL = Kpoint._ATOL_KDIFF

    stars = len(ibz) * [None]
    coord_stars = len(ibz) * [None]
    norm2star = {}
    for i, kirr in enumerate(ibz):
        stars[i] = kirr.get_star(fm_symmops)
        norm2star[kirr.norm] = i
        fc = [k.frac_coords for k in stars[i]]
        coord_stars[i] = np.reshape(fc, (-1,3))

    #fullbz_coords = np.empty((len(bz), 3))
    #ibz_coords = np.empty((len(ibz), 3))

    #for i, k in enumerate(ibz):
    #    ibz_coords[i][:] = k.frac_coords

    #for i, kfull in enumerate(bz):
    #    fullbz_coords[i][:] = k.frac_coords

    #for ik_full in range(len(bz)):
    #    found = False
    #    fullk = fullbz_coords[ik_full]
    #    for ik_irr in range(len(ibz)):
    #        kdiff = coord_stars[ik_irr] - fullk
    #        #print(kdiff.shape)
    #        mask = abs(np.around(kdiff) - kdiff) < ATOL
    #        for m in mask:
    #            if m.all():
    #                found = True
    #                bz2ibz[ik_full] = ik_irr
    #                break
    #        #if mask.any():
    #        #    found = True
    #        #    bz2ibz[ik_full] = ik_irr
    #        #    break

    #    if not found:
    #        raise ValueError("Full k-point: %s not found" % kfull)

    #return bz2ibz
    #print("done")
    #import sys
    #sys.exit(1)

    miss = 0
    for bz_idx, kfull in enumerate(bz):
        if bz2ibz[bz_idx] != -1:
            continue
        found = False

        kfull_coords = kfull.frac_coords

        #kf_norm = kfull.norm
        #if kf_norm in norm2star:
        #    ibz_idx = norm2star[kf_norm]
        #    star = stars[ibz_idx]
        #    #    if kfull in star:
        #    #        found = True
        #    #        bz2ibz[bz_idx] = ibz_idx
        #    for k in star:
        #        kdiff = k.frac_coords - kfull_coords
        #        if all(abs(np.around(kdiff) - kdiff) < ATOL):
        #            found = True
        #            bz2ibz[bz_idx] = ibz_idx
        #            break
        #if found:
        #    continue

        for ibz_idx, kirr in enumerate(ibz):
            #if kfull in kirr.iter_star(fm_symmops):
            #    bz2ibz[bz_idx] = ibz_idx
            #    break

            #kdiffs = np.reshape([k.frac_coords - kfull_coords for k in stars[ibz_idx]], (-1,3))

            for k in stars[ibz_idx]:
                kdiff = k.frac_coords - kfull_coords
                #kdiff = np.ascontiguousarray(kdiff)
                #assert kdiff.iscontiguous
                #int_kdiff = np.array([round(c) for c in kdiff], dtype=np.float)
                #int_kdiff = np.array([c for c in kdiff], dtype=np.float)
                #if all(abs(int_kdiff - kdiff) < ATOL):
                if all(abs(np.round(kdiff) - kdiff) < ATOL):
                #if all(abs(kdiff) < ATOL):
                    found = True
                    bz2ibz[bz_idx] = ibz_idx
                    break

            #if kfull in stars[ibz_idx]:
            #    found = True
            #    bz2ibz[bz_idx] = ibz_idx
            #    break

        if not found:
            #print("Full k-point: %s not found" % kfull))
            miss += 1

    if miss:
        err_msg = "Cannot map BZ %d/%d missing" % (miss, len(bz))
        print(err_msg)
        #raise ValueError(err_msg)

    print("Done")
    return bz2ibz

#########################################################################################


class Kmesh(object):
    """
    This object describes the sampling of the full Brillouin zone.

    A Kmesh object has a set of irreduble point, information on how
    to generate the mesh in the full BZ. It also provides methods
    to symmetrized and visualize quantities in k-space.
    """
    def __init__(self, structure, mpdivs, shifts, ibz):
        """
        Args:
            structure:
                `Structure` instance.
            mpdivs:
                number of Monkhorst-Pack divisions.
            shifts:
                shifts of the mesh.
            ibz:
                k-points in the irreducible wedge.

        Provides methods to  symmetrize k-dependent quantities with the full
        symmetry of the structure. e.g. bands, occupation factors, phonon frequencies.
        """
        self.structure = structure
        self.mpdivs = np.asarray(mpdivs)
        self.nx, self.ny, self.nz = self.mpdivs

        self.shifts = np.reshape(shifts, (-1,3))
        assert len(self.shifts) == 1
        self.ibz = ibz

        grids_1d = 3 * [None]
        for i in range(3):
            grids_1d[i] = np.arange(0, self.mpdivs[i])
        self.grids_1d = tuple(grids_1d)

        # Hack needed because the python code in map_mesh2ibz is very slow!
        # I need a C version, possibly included in spglib.
        bkp_file = "bzmap.arr.npy"
        if os.path.exists(bkp_file):
            self.bzmap = np.load(bkp_file)

        else:
            # Build the mapping BZ --> IBZ for the different grids.
            # map[i,j,k] --> index of the point in the IBZ
            self.bzmap = map_mesh2ibz(structure, self.mpdivs, self.shifts, ibz)

            np.save(bkp_file, self.bzmap)

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
                ibz_idx = self.bzmap[x,y,z0]
                plane[x,y] = values_ibz[ibz_idx]

        kx, ky = np.meshgrid(kx, ky)
        return kx, ky, plane

