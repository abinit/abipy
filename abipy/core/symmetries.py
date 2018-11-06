# coding: utf-8
"""Objects used to deal with symmetry operations in crystals."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import abc
import warnings
import collections
import six
import numpy as np
import spglib

from six.moves import cStringIO
from monty.string import is_string
from monty.itertools import iuptri
from monty.functools import lazy_property
from monty.collections import dict2namedtuple
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
try:
    from pymatgen.util.serialization import SlotPickleMixin
except:
    from pymatgen.serializers.pickle_coders import SlotPickleMixin
from abipy.core.kpoints import wrap_to_ws, issamek, has_timrev_from_kptopt
from abipy.iotools import as_etsfreader


__all__ = [
    "LatticeRotation",
    "AbinitSpaceGroup",
]

def wrap_in_ucell(x):
    """
    Transforms x in its corresponding reduced number in the interval [0,1[."
    """
    return x % 1


def is_integer(x, atol=1e-08):
    """
    True if all x is integer within the absolute tolerance atol.

    >>> is_integer([1., 2.])
    True
    >>> is_integer(1.01, atol=0.011)
    True
    >>> is_integer([1.01, 2])
    False
    """
    int_x = np.around(x)
    return np.allclose(int_x, x, atol=atol)


def mati3inv(mat3, trans=True):
    """
    Invert and transpose orthogonal 3x3 matrix of INTEGER elements.

    Args:
        mat3: (3, 3) matrix-like object with integer elements

    Returns:
        |numpy-array| with the TRANSPOSE of the inverse of mat3 if trans==True.
        If trans==False, the inverse of mat3 is returned.

    .. note::

       Used for symmetry operations. This function applies to *ORTHOGONAL* matrices only.
       Since these form a group, inverses are also integer arrays.
    """
    mat3 = np.reshape(np.array(mat3, dtype=np.int), (3, 3))

    mit = np.empty((3, 3), dtype=np.int)
    mit[0,0] = mat3[1,1] * mat3[2,2] - mat3[2,1] * mat3[1,2]
    mit[1,0] = mat3[2,1] * mat3[0,2] - mat3[0,1] * mat3[2,2]
    mit[2,0] = mat3[0,1] * mat3[1,2] - mat3[1,1] * mat3[0,2]
    mit[0,1] = mat3[2,0] * mat3[1,2] - mat3[1,0] * mat3[2,2]
    mit[1,1] = mat3[0,0] * mat3[2,2] - mat3[2,0] * mat3[0,2]
    mit[2,1] = mat3[1,0] * mat3[0,2] - mat3[0,0] * mat3[1,2]
    mit[0,2] = mat3[1,0] * mat3[2,1] - mat3[2,0] * mat3[1,1]
    mit[1,2] = mat3[2,0] * mat3[0,1] - mat3[0,0] * mat3[2,1]
    mit[2,2] = mat3[0,0] * mat3[1,1] - mat3[1,0] * mat3[0,1]

    dd = mat3[0,0] * mit[0,0] + mat3[1,0] * mit[1,0] + mat3[2,0] * mit[2,0]

    # Make sure matrix is not singular
    if dd == 0:
        raise ValueError("Attempting to invert integer array: %s\n ==> determinant is zero." % str(mat3))

    mit = mit // dd
    if trans:
        return mit
    else:
        return mit.T.copy()


def _get_det(mat):
    """
    Return the determinant of a 3x3 rotation matrix mat.

    raises:
        ValueError if abs(det) != 1.
    """
    det = mat[0,0]* (mat[1,1]*mat[2,2] - mat[1,2]*mat[2,1])\
        - mat[0,1]* (mat[1,0]*mat[2,2] - mat[1,2]*mat[2,0])\
        + mat[0,2]* (mat[1,0]*mat[2,1] - mat[1,1]*mat[2,0])

    if abs(det) != 1:
        raise ValueError("Determinant must be +-1 while it is %s" % det)

    return det


def indsym_from_symrel(symrel, tnons, structure, tolsym=1e-8):
    r"""
    For each symmetry operation, find the number of the position to
    which each atom is sent in the unit cell by the INVERSE of the
    symmetry operation inv(symrel); i.e. this is the atom which, when acted
    upon by the given symmetry element isym, gets transformed into atom iatom.
    indirect indexing array for atoms, see symatm.F90.

    $ R^{-1} (xred(:,iat) - \tau) = xred(:,iat_sym) + R_0 $
    * indsym(4,  isym,iat) gives iat_sym in the original unit cell.
    * indsym(1:3,isym,iat) gives the lattice vector $R_0$.

    Args:
        symrel: int (nsym,3,3) array with real space symmetries expressed in reduced coordinates.
        tnons: float (nsym, 3) array with nonsymmorphic translations for each symmetry.
        structure: |Structure| object.
        tolsym: tolerance for the symmetries

    Returns:
    """
    natom = len(structure)
    nsym = len(symrel)
    xred = np.array([site.frac_coords for site in structure], dtype=float)
    typat = {i: site.specie.symbol for i, site in enumerate(structure)}

    rm1_list = np.empty_like(symrel)
    for isym in range(nsym):
        rm1_list[isym] = mati3inv(symrel[isym], trans=False)

    # Start testmn out at large value
    testmn = 1000000
    err = 0.0
    indsym = np.empty((natom, nsym, 4))

    # Implementation is similar to Abinit routine (including the order of the loops)
    for isym in range(nsym):
        for iatom in range(natom):
            tratm = np.matmul(rm1_list[isym], xred[iatom] - tnons[isym])
            # Loop through atoms, when types agree, check for agreement after primitive translation
            for jatm in range(natom):
                if typat[jatm] != typat[iatom]: continue
                test_vec = tratm - xred[jatm]
                # Find nearest integer part of difference
                trans = np.rint(test_vec)
                # Check whether, after translation, they agree
                test_vec = test_vec - trans
                diff = np.abs(test_vec).sum()
                # Abinit uses 1e-10 but python seems to require a slightly larger value.
                #if diff < 1e-10:
                if diff < 1e-9:
                    difmin = test_vec
                    indsym[iatom, isym, :3] = trans
                    indsym[iatom, isym, 3] = jatm
                    # Break out of loop when agreement is within tolerance
                    break
                else:
                    # Keep track of smallest difference if greater than tol10
                    if diff < testmn:
                        testmn = diff
                        # Note that abs() is not taken here
                        difmin = test_vec
                        indsym[iatom, isym, :3] = trans
                        indsym[iatom, isym, 3] = jatm

        # Keep track of maximum difference between transformed coordinates and nearest "target" coordinate
        difmax = np.abs(difmin).max()
        err = max(err, difmax)
        if difmax > tolsym:
            cprint("""
Trouble finding symmetrically equivalent atoms.
Applying inverse of symm number {isym} to atom number {iatom} of typat',typat(iatom) gives tratom=',tratom(1:3)
This is further away from every atom in crystal than the allowed tolerance.
The inverse symmetry matrix is',symrec(1,1:3,isym),ch10,&
                               ',symrec(2,1:3,isym),ch10,&
                               ',symrec(3,1:3,isym)
and the nonsymmorphic transl. tnons =',(tnons(mu,isym),mu=1,3)
The nearest coordinate differs by',difmin(1:3) for indsym(nearest atom)=',indsym(4,isym,iatom)

This indicates that when symatm attempts to find atoms symmetrically
related to a given atom, the nearest candidate is further away than some tolerance.
Should check atomic coordinates and symmetry group input data.
""".format(), "red")

    if err > tolsym:
        raise ValueError("maximum err %s is larger than tolsym: %s" % (err, tolsym))

    return indsym


@six.add_metaclass(abc.ABCMeta)
class Operation(object):
    """
    Abstract base class that defines the methods that must be
    implemented by the concrete class representing some sort of operation
    """
    @abc.abstractmethod
    def __eq__(self, other):
        """O1 == O2"""

    def __ne__(self, other):
        return not (self == other)

    @abc.abstractmethod
    def __mul__(self, other):
        """O1 * O2"""

    @abc.abstractmethod
    def __hash__(self):
        """Operation can be used as dictionary keys."""

    @abc.abstractmethod
    def inverse(self):
        """Returns the inverse of self."""

    def opconj(self, other):
        """Returns X^-1 S X where X is the other symmetry operation."""
        return other.inverse() * self * other

    @abc.abstractproperty
    def isE(self):
        """True if self is the identity operator"""

    #def commute(self, other)
    #    return self * other == other * self

    #def commutator(self, other)
    #    return self * other - other * self

    #def anticommute(self, other)
    #    return self * other == - other * self

    #def direct_product(self, other)


class SymmOp(Operation, SlotPickleMixin):
    """
    Crystalline symmetry.
    """
    _ATOL_TAU =  1e-8

    __slots__ = [
        "rot_r",
        "rotm1_r",
        "tau",
        "time_sign",
        "afm_sign",
        "rot_g",
        "_det",
        "_trace",
    ]

    # TODO: Add lattice?
    def __init__(self, rot_r, tau, time_sign, afm_sign, rot_g=None):
        """
        This object represents a space group symmetry i.e. a symmetry of the crystal.

        Args:
            rot_r: (3,3) integer matrix with the rotational part in real space in reduced coordinates (C order).
            tau: fractional translation in reduced coordinates.
            time_sign: -1 if time reversal can be used, +1 otherwise.
            afm_sign: anti-ferromagnetic part [+1, -1].
        """
        rot_r = np.asarray(rot_r)

        # Store R and R^{-1} in real space.
        self.rot_r, self.rotm1_r = rot_r, mati3inv(rot_r, trans=False)
        self.tau = np.asarray(tau)

        self.afm_sign, self.time_sign = afm_sign, time_sign
        assert afm_sign in [-1, 1] and time_sign in [-1, 1]

        # Compute symmetry matrix in reciprocal space: S = R^{-1t}
        if rot_g is None:
            self.rot_g = mati3inv(rot_r, trans=True)
        else:
            assert np.all(rot_g == mati3inv(rot_r, trans=True))
            self.rot_g = rot_g

    # operator protocol.
    def __eq__(self, other):
        # Note the two fractional traslations are equivalent if they differ by a lattice vector.
        return (np.all(self.rot_r == other.rot_r) and
                is_integer(self.tau - other.tau, atol=self._ATOL_TAU) and
                self.afm_sign == other.afm_sign and
                self.time_sign == other.time_sign)

    def __mul__(self, other):
        """
        Returns a new :class:`SymmOp` which is equivalent to apply the "other" :class:`SymmOp`
        followed by this one i.e:

        {R,t} {S,u} = {RS, Ru + t}
        """
        return self.__class__(rot_r=np.dot(self.rot_r, other.rot_r),
                              tau=self.tau + np.dot(self.rot_r, other.tau),
                              time_sign=self.time_sign * other.time_sign,
                              afm_sign=self.afm_sign * other.afm_sign)

    def __hash__(self):
        """
        :class:`Symmop` can be used as keys in dictionaries.
        Note that the hash is computed from integer values.
        """
        return int(8 * self.trace + 4 * self.det + 2 * self.time_sign)

    def inverse(self):
        """Returns inverse of transformation i.e. {R^{-1}, -R^{-1} tau}."""
        return self.__class__(rot_r=self.rotm1_r,
                              tau=-np.dot(self.rotm1_r, self.tau),
                              time_sign=self.time_sign,
                              afm_sign=self.afm_sign)

    @lazy_property
    def isE(self):
        """True if identity operator."""
        return (np.all(self.rot_r == np.eye(3, dtype=np.int)) and
                is_integer(self.tau, atol=self._ATOL_TAU) and
                self.time_sign == 1 and
                self.afm_sign == 1)
    # end operator protocol.

    #@lazy_property
    #def order(self):
    #    """Order of the operation."""
    #    n = 0
    #    o = self
    #    while m < 1000:
    #        if o.isE: return n
    #        n += 1
    #        o = self * o
    #    else:
    #        raise ValueError("Cannot find order")

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        def vec2str(vec):
            return "%2d,%2d,%2d" % tuple(v for v in vec)

        s = ""
        for i in range(3):
            s += "[" + vec2str(self.rot_r[i]) + ", %.3f]  " % self.tau[i] + "[" + vec2str(self.rot_g[i]) + "] "
            if i == 2:
                s += ", time_sign = %+1d, afm_sign = %+1d, det = %+1d" % (self.time_sign, self.afm_sign, self.det)
            s += "\n"

        return s

    @lazy_property
    def is_symmorphic(self):
        """True if the fractional translation is non-zero."""
        return np.any(np.abs(self.tau) > 0.0)

    @lazy_property
    def det(self):
        """Determinant of the rotation matrix [-1, +1]."""
        return _get_det(self.rot_r)

    @lazy_property
    def trace(self):
        """Trace of the rotation matrix."""
        return self.rot_r.trace()

    @lazy_property
    def is_proper(self):
        """True if the rotational part has determinant == 1."""
        return self.det == +1

    @lazy_property
    def has_timerev(self):
        """True if symmetry contains the time-reversal operator."""
        return self.time_sign == -1

    @lazy_property
    def is_fm(self):
        """True if self if ferromagnetic symmetry."""
        return self.afm_sign == +1

    @lazy_property
    def is_afm(self):
        """True if self if anti-ferromagnetic symmetry."""
        return self.afm_sign == -1

    def rotate_k(self, frac_coords, wrap_tows=False):
        """
        Apply the symmetry operation to the k-point given in reduced coordinates.

        Sk is wrapped to the first Brillouin zone if wrap_tows is True.
        """
        sk = np.dot(self.rot_g, frac_coords) * self.time_sign

        return wrap_to_ws(sk) if wrap_tows else sk

    def preserve_k(self, frac_coords, ret_g0=True):
        """
        Check if the operation preserves the k-point modulo a reciprocal lattice vector.

        Args:
            frac_coords: Fractional coordinates of the k-point
            ret_g0: False if only the boolean result is wanted.

        Returns:
            bool, g0 = S(k) - k

            bool is True is self preserves k and g0 is an integer vector.
        """
        sk = self.rotate_k(frac_coords, wrap_tows=False)

        if ret_g0:
            return issamek(sk, frac_coords), np.array(np.round(sk - frac_coords), dtype=np.int)
        else:
            return issamek(sk, frac_coords)

    def rotate_r(self, frac_coords, in_ucell=False):
        """
        Apply the symmetry operation to a point in real space given in reduced coordinates.

        .. NOTE::

            We use the convention: symmop(r) = R^{-1] (r - tau)
        """
        rotm1_rmt = np.dot(self.rotm1_r, frac_coords - self.tau)

        return wrap_in_ucell(rotm1_rmt) if in_ucell else rotm1_rmt


class OpSequence(collections.Sequence):
    """
    Mixin class providing the basic method that are common to containers of operations.
    """
    def __len__(self):
        return len(self._ops)

    def __iter__(self):
        return self._ops.__iter__()

    def __getitem__(self, slice):
        return self._ops[slice]

    def __contains__(self, op):
        return op in self._ops

    def __eq__(self, other):
        """
        Equality test.

        .. warning::

            The order of the operations in self and  in other is not relevant.
        """
        if other is None: return False
        if len(self) != len(other):
            return False

        # Check if each operation in self is also present in other.
        # The order is irrelevant.
        founds = []
        for i, op in enumerate(self):
            if op not in other: return False
            founds.append(i)

        if len(set(founds)) == len(founds):
            return True

        warnings.warn("self contains duplicated ops! Likely a bug!")
        return False

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        lines = [str(op) for op in self]
        return "\n".join(lines)

    def show_ops(self, stream=sys.stdout):
        lines = [str(op) for op in self]
        stream.writelines("\n".join(lines))

    def count(self, op):
        """Returns the number of occurences of operation op in self."""
        return self._ops.count(op)

    def index(self, op):
        """
        Return the (first) index of operation op in self.

        Raises:
            ValueError if not found.
        """
        return self._ops.index(op)

    def find(self, op):
        """Return the (first) index of op in self. -1 if not found."""
        try:
            return self.index(op)
        except ValueError:
            return -1

    def is_group(self):
        """True if this set of operations represent a group."""
        check = 0

        # Identity must be present.
        if [op.isE for op in self].count(True) != 1:
            check += 1

        # The inverse must be in the set.
        if [op.inverse() in self for op in self].count(True) != len(self):
            check += 2

        # The product of two members must be in the set.
        op_prods = [op1 * op2 for op1 in self for op2 in self]

        d = self.asdict()
        for op12 in op_prods:
            if op12 not in d:
                print("op12 not in group\n %s" % str(op12))
                check += 1

        return check == 0

    def is_commutative(self):
        """True if all operations commute with each other."""
        for op1, op2 in iuptri(self, diago=False):
            if op1 * op2 != op2 * op1: return False
        return True

    def is_abelian_group(self):
        """True if commutative group."""
        return self.is_commutative() and self.is_group()

    def asdict(self):
        """
        Returns a dictionary where the keys are the symmetry operations and
        the values are the indices of the operations in the iterable.
        """
        return {op: idx for idx, op in enumerate(self)}

    #def is_subset(self, other)
    #    indmap = {}
    #    for i, op in self:
    #        j = other.find(op)
    #        if j != -1: indmap[i] = j
    #    return indmap

    #def is_superset(self, other)

    @lazy_property
    def mult_table(self):
        """
        Given a set of nsym 3x3 operations which are supposed to form a group,
        this routine constructs the multiplication table of the group.
        mtable[i,j] gives the index of the product S_i * S_j.
        """
        mtable = np.empty((len(self), len(self)), dtype=np.int)

        d = self.asdict()
        for i, op1 in enumerate(self):
            for j, op2 in enumerate(self):
                op12 = op1 * op2
                # Save the index of op12 in self
                try:
                    index = d[op12]
                except KeyError:
                    index = None
                mtable[i, j] = index

        return mtable

    @property
    def num_classes(self):
        """Number of classes."""
        return len(self.class_indices)

    @lazy_property
    def class_indices(self):
        """
        A class is defined as the set of distinct elements obtained by
        considering for each element, S, of the group all its conjugate
        elements X^-1 S X where X ranges over all the elements of the group.

        Returns:
            Nested list l = [cls0_indices, cls1_indices, ...] where each sublist
            contains the indices of the class. len(l) equals the number of classes.
        """
        found, class_indices = len(self) * [False], [[] for i in range(len(self))]

        num_classes = -1
        for ii, op1 in enumerate(self):
            if found[ii]: continue
            num_classes += 1

            for jj, op2 in enumerate(self):
                # Form conjugate and search it among the operations
                # that have not been found yet.
                op1_conj = op1.opconj(op2)

                for kk, op3 in enumerate(self):
                    if not found[kk] and op1_conj == op3:
                        found[kk] = True
                        class_indices[num_classes].append(kk)

        class_indices = class_indices[:num_classes + 1]
        assert sum(len(c) for c in class_indices) == len(self)
        return class_indices

    def groupby_class(self, with_inds=False):
        """
        Iterate over the operations grouped in symmetry classes.

        Args:
            with_inds: If True, [op0, op1, ...], [ind_op0, ind_op1, ...] is returned.
        """
        if with_inds:
            for indices in self.class_indices:
                yield [self[i] for i in indices], indices
        else:
            for indices in self.class_indices:
                yield [self[i] for i in indices]


class AbinitSpaceGroup(OpSequence):
    """
    Container storing the space group symmetries.
    """

    def __init__(self, spgid, symrel, tnons, symafm, has_timerev, inord="C"):
        """
        Args:
            spgid (int): space group number (from 1 to 232, 0 if cannot be specified).
            symrel: (nsym,3,3) array with the rotational part of the symmetries in real
                space (reduced coordinates are assumed, see also `inord` for the order.
            tnons: (nsym,3) array with fractional translation in reduced coordinates.
            symafm: (nsym) array with +1 for Ferromagnetic symmetry and -1 for AFM
            has_timerev: True if time-reversal symmetry is included.
            inord: storage order of mat in symrel[:]. If inord == "F", mat.T is stored
                as matrices are always stored in C-order. Use inord == "F" if you have
                read symrel from an external file produced by abinit.

        .. note::

            All the arrays are stored in C-order. Use as_fortran_arrays to extract data
            that can be passes to Fortran routines.
        """
        self.spgid = spgid
        assert 233 > self.spgid >= 0

        # Time reversal symmetry.
        self._has_timerev = has_timerev
        self._time_signs = [+1, -1] if self.has_timerev else [+1]

        self._symrel, self._tnons, self._symafm = list(map(np.asarray, (symrel, tnons, symafm)))

        if len(self.symrel) != len(self.tnons) or len(self.symrel) != len(self.symafm):
            raise ValueError("symrel, tnons and symafm must have equal shape[0]")

        inord = inord.upper()
        assert inord in ["F", "C"]
        if inord == "F":
            # Fortran to C.
            for isym in range(len(self.symrel)):
                self._symrel[isym] = self._symrel[isym].T

        self._symrec = self._symrel.copy()
        for isym in range(len(self.symrel)):
            self._symrec[isym] = mati3inv(self.symrel[isym], trans=True)

        all_syms = []
        for time_sign in self._time_signs:
            for isym in range(len(self.symrel)):
                all_syms.append(SymmOp(rot_r=self.symrel[isym],
                                       tau=self.tnons[isym],
                                       time_sign=time_sign,
                                       afm_sign=self.symafm[isym],
                                       rot_g=self.symrec[isym]))
        self._ops = tuple(all_syms)

    @classmethod
    def from_ncreader(cls, r, inord="F"):
        """
        Builds the object from a netcdf reader
        """
        kptopt = int(r.read_value("kptopt", default=1))
        symrel = r.read_value("reduced_symmetry_matrices")

        return cls(spgid=r.read_value("space_group"),
                   symrel=symrel,
                   tnons=r.read_value("reduced_symmetry_translations"),
                   symafm=r.read_value("symafm"),
                   has_timerev=has_timrev_from_kptopt(kptopt),
                   inord=inord)

    @classmethod
    def from_file(cls, ncfile, inord="F"):
        """
        Initialize the object from a Netcdf file.
        """
        r, closeit = as_etsfreader(ncfile)
        new = cls.from_ncreader(r)
        if closeit:
            file.close()

        return new

    @classmethod
    def from_structure(cls, structure, has_timerev=True, symprec=1e-5, angle_tolerance=5):
        """
        Takes a |Structure| object. Uses spglib to perform various symmetry finding operations.

        Args:
            structure: |Structure| object.
            has_timerev: True is time-reversal symmetry is included.
            symprec: Tolerance for symmetry finding.
            angle_tolerance: Angle tolerance for symmetry finding.

        .. warning::

            AFM symmetries are not supported.
        """
        # Call spglib to get the list of symmetry operations.
        spga = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
        data = spga.get_symmetry_dataset()
        symrel = data["rotations"]

        return cls(spgid=data["number"],
                   symrel=symrel,
                   tnons=data["translations"],
                   symafm=len(symrel) * [1],
                   has_timerev=has_timerev,
                   inord="C")

    def __repr__(self):
        return "spgid: %d, num_spatial_symmetries: %d, has_timerev: %s, symmorphic: %s" % (
            self.spgid, self.num_spatial_symmetries, self.has_timerev, self.is_symmorphic)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = ["spgid: %d, num_spatial_symmetries: %d, has_timerev: %s, symmorphic: %s" % (
            self.spgid, self.num_spatial_symmetries, self.has_timerev, self.is_symmorphic)]
        app = lines.append

        if verbose > 1:
            for op in self.symmops(time_sign=+1):
                app(str(op))

        return "\n".join(lines)

    @lazy_property
    def is_symmorphic(self):
        """True if there's at least one operation with non-zero fractional translation."""
        return any(op.is_symmorphic for op in self)

    @property
    def has_timerev(self):
        """True if time-reversal symmetry is present."""
        return self._has_timerev

    @property
    def symrel(self):
        """
        [nsym, 3, 3] int array with symmetries in reduced coordinates of the direct lattice.
        """
        return self._symrel

    @property
    def tnons(self):
        """
        [nsym, 3] float array with fractional translations in reduced coordinates of the direct lattice.
        """
        return self._tnons

    @property
    def symrec(self):
        """
        [nsym, 3, 3] int array with symmetries in reduced coordinates of the reciprocal lattice.
        """
        return self._symrec

    @property
    def symafm(self):
        """[nsym] int array with +1 if FM or -1 if AFM symmetry."""
        return self._symafm

    @property
    def num_spatial_symmetries(self):
        fact = 2 if self.has_timerev else 1
        return int(len(self) / fact)

    @property
    def afm_symmops(self):
        """Tuple with antiferromagnetic symmetries."""
        return self.symmops(time_sign=None, afm_sign=-1)

    @property
    def fm_symmops(self):
        """Tuple of ferromagnetic symmetries."""
        return self.symmops(time_sign=None, afm_sign=+1)

    def symmops(self, time_sign=None, afm_sign=None):
        """
        Args:
            time_sign: If specified, only symmetries with time-reversal sign time_sign are returned.
            afm_sign: If specified, only symmetries with anti-ferromagnetic part afm_sign are returned.

        returns:
            tuple of :class:`SymmOp` instances.
        """
        symmops = []
        for sym in self._ops:
            gotit = True

            if time_sign:
                assert time_sign in self._time_signs
                gotit = gotit and sym.time_sign == time_sign

            if afm_sign:
                assert afm_sign in [-1,+1]
                gotit = gotit and sym.afm_sign == afm_sign

            if gotit:
                symmops.append(sym)

        return tuple(symmops)

    def symeq(self, k1_frac_coords, k2_frac_coords, atol=None):
        """
        Test whether two k-points in fractional coordinates are symmetry equivalent
        i.e. if there's a symmetry operations TO (including time-reversal T, if present)
	such that::

            TO(k1) = k2 + G0

	Return: namedtuple with::

            isym: The index of the symmetry operation such that TS(k1) = k2 + G0
                Set to -1 if k1 and k2 are not related by symmetry.
            op: Symmetry operation.
            g0: numpy vector.
        """
        for isym, sym in enumerate(self):
            sk_coords = sym.rotate_k(k1_frac_coords, wrap_tows=False)
            if issamek(sk_coords, k2_frac_coords, atol=atol):
                g0 = sym.rotate_k(k1_frac_coords) - k2_frac_coords
                return dict2namedtuple(isym=isym, op=self[isym], g0=g0)

        return dict2namedtuple(isym=-1, op=None, g0=None)

    def find_little_group(self, kpoint):
        """
        Find the little group of the kpoint.

        Args:
            kpoint: Accept vector with the reduced coordinates or :class:`Kpoint` object.

        Returns:
            :class:`LittleGroup` object.
        """
        if hasattr(kpoint, "frac_coords"):
            frac_coords = kpoint.frac_coords
        else:
            frac_coords = np.reshape(kpoint, (3))

        to_spgrp, g0vecs = [], []

        # Exclude AFM operations.
        for isym, symmop in enumerate(self.fm_symmops):
            is_same, g0 = symmop.preserve_k(frac_coords)
            if is_same:
                to_spgrp.append(isym)
                g0vecs.append(g0)

        # List with the symmetry operation that preserve the kpoint.
        k_symmops = [self[i] for i in to_spgrp]
        return LittleGroup(kpoint, k_symmops, g0vecs)

# FIXME To maintain backward compatibility.
SpaceGroup = AbinitSpaceGroup


class LittleGroup(OpSequence):

    def __init__(self, kpoint, symmops, g0vecs):
        """
        k_symmops, g0vecs, indices

        k_symmops is a tuple with the symmetry operations that preserve the k-point i.e. Sk = k + G0
        g0vecs is the tuple for G0 vectors for each operation in k_symmops
        """
        self.kpoint = kpoint
        self._ops = symmops
        self.g0vecs = np.reshape(g0vecs, (-1, 3))
        assert len(self.symmops) == len(self.g0vecs)

        # Find the point group of k so that we know how to access the Bilbao database.
        # (note that operations are in reciprocal space, afm and time_reversal are taken out
        krots = np.array([o.rot_g for o in symmops if not o.has_timerev])
        self.kgroup = LatticePointGroup(krots)

    @lazy_property
    def is_symmorphic(self):
        """True if there's at least one operation with non-zero fractional translation."""
        return any(op.is_symmorphic for op in self)

    @property
    def symmops(self):
        return self._ops

    @lazy_property
    def on_bz_border(self):
        """
        True if the k-point is on the border of the BZ.
        """
        frac_coords = np.array(self.kpoint)
        kreds = wrap_to_ws(frac_coords)
        diff = np.abs(np.abs(kreds) - 0.5)
        return np.any(diff < 1e-8)

    def iter_symmop_g0(self):
        for symmop, g0 in zip(self.symmops, self.g0vecs):
            yield symmop, g0

    def __repr__(self):
        return "Kpoint Group: %s, Kpoint: %s" % (self.kgroup, self.kpoint)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation of little group."""
        lines = ["Kpoint-group: %s, Kpoint: %s, Symmorphic: %s" % (self.kgroup, self.kpoint, self.is_symmorphic)]
        app = lines.append
        app(" ")

        # Add character_table from Bilbao database.
        bilbao_ptgrp = bilbao_ptgroup(self.kgroup.sch_symbol)
        app(bilbao_ptgrp.to_string(verbose=verbose))
        app("")

        # Write warning if non-symmorphic little group with k-point at zone border.
        if self.is_symmorphic and self.on_bz_border:
            app("WARNING: non-symmorphic little group with k at zone-border.")
            app("Electronic states cannot be classified with this character table.")

        return "\n".join(lines)

    #def iter_symmop_g0_byclass(self):


class LatticePointGroup(OpSequence):

    def __init__(self, rotations):
        rotations = np.reshape(rotations, (-1, 3, 3))
        self._ops = [LatticeRotation(rot) for rot in rotations]

        # Call spglib to get the Herm symbol.
        # (symbol, pointgroup_number, transformation_matrix)
        herm_symbol, ptg_num, trans_mat = spglib.get_pointgroup(rotations)
        # Remove blanks from C string.
        self.herm_symbol = herm_symbol.strip()
        #print(self.herm_symbol, ptg_num, trans_mat)

        if self.sch_symbol is None:
            raise ValueError("Cannot detect point group symbol! Got sch_symbol = %s" % self.sch_symbol)

    #@classmethod
    #def from_vectors(cls, vectors)

    def __repr__(self):
        return "%s: %s, %s (%d)" % (self.__class__.__name__, self.herm_symbol, self.sch_symbol, self.spgid)

    def __str__(self):
        return "%s, %s (%d)" % (self.herm_symbol, self.sch_symbol, self.spgid)

    @property
    def sch_symbol(self):
        """Schoenflies symbol"""
        return herm2sch(self.herm_symbol)

    @property
    def spgid(self):
        """ID in the space group table."""
        return sch2spgid(self.sch_symbol)


class LatticeRotation(Operation):
    """
    This object defines a pure rotation of the lattice (proper, improper, mirror symmetry)
    that is a rotation which is compatible with a lattice. The rotation matrix is
    expressed in reduced coordinates, therefore its elements are integers.

    See:
        http://xrayweb2.chem.ou.edu/notes/symmetry.html#rotation

    .. note::

        This object is immutable and therefore we do not inherit from |numpy-array|.
    """
    _E3D = np.identity(3,  np.int)

    def __init__(self, mat):
        self.mat = np.asarray(mat, dtype=np.int)
        self.mat.shape = (3, 3)

    def _find_order_and_rootinv(self):
        """
        Returns the order of the rotation and if self is a root of the inverse.
        """
        order, root_inv = None, 0
        for ior in range(1, 7):
            rn = self ** ior

            if rn.isE:
                order = ior
                break

            if rn.isI:
                root_inv = ior

        if order is None:
            raise ValueError("LatticeRotation is not a root of unit!")

        return order, root_inv

    def __repr__(self):
        return self.name

    #def __str__(self):
    #    lines = "Rotation: " + str(self.order) + ", versor: " + str(self.versor) + ",
    #    lines.append(str(self.mat))
    #    return "\n".join(lines)

    # operator protocol.
    def __eq__(self, other):
        return np.allclose(self.mat, other.mat)

    def __mul__(self, other):
        return self.__class__(np.matmul(self.mat, other.mat))

    def __hash__(self):
        return int(8 * self.trace + 4 * self.det)

    def inverse(self):
        """
        Invert an orthogonal 3x3 matrix of INTEGER elements.
        Note use of integer arithmetic. Raise ValueError if not invertible.
        """
        return self.__class__(mati3inv(self.mat, trans=False))

    @lazy_property
    def isE(self):
        """True if it is the identity"""
        return np.allclose(self.mat, self._E3D)
    # end operator protocol.

    # Implement the unary arithmetic operations (+, -)
    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.mat)

    def __pow__(self, intexp, modulo=1):
        if intexp ==  0: return self.__class__(self._E3D)
        if intexp  >  0: return self.__class__(np.linalg.matrix_power(self.mat, intexp))
        if intexp == -1: return self.inverse()
        if intexp  <  0: return self.__pow__(-intexp).inverse()
        raise TypeError("type %s is not supported in __pow__" % type(intexp))

    @property
    def order(self):
        """Order of the rotation."""
        try:
            return self._order
        except AttributeError:
            self._order, self._root_inv = self._find_order_and_rootinv()
            return self._order

    @property
    def root_inv(self):
        try:
            return self._root_inv
        except AttributeError:
            self._order, self._root_inv = self._find_order_and_rootinv()
            return self._root_inv

    @lazy_property
    def det(self):
        """Return the determinant of a symmetry matrix mat[3,3]. It must be +-1"""
        return _get_det(self.mat)

    @lazy_property
    def trace(self):
        """The trace of the rotation matrix"""
        return self.mat.trace()

    @lazy_property
    def is_proper(self):
        """True if proper rotation"""
        return self.det == 1

    @lazy_property
    def isI(self):
        """True if self is the inversion operation."""
        return np.allclose(self.mat, -self._E3D)

    @lazy_property
    def name(self):
        # Sign of the determinant (only if improper)
        name = "-" if self.det == -1 else ""
        name += str(self.order)
        # Root of inverse?
        name += "-" if self.root_inv != 0 else "+"

        return name

    #@property
    #def rottype(self):
    #    """
    #    Receive a 3x3 orthogonal matrix and reports its type:
    #        1 Identity
    #        2 Inversion
    #        3 Proper rotation of an angle <> 180 degrees
    #        4 Proper rotation of 180 degrees
    #        5 Mirror symmetry
    #        6 Improper rotation
    #    """
    #    # Treat identity and inversion first
    #    if self.isE: return 1
    #    if self.isI: return 2
    #
    #    if self.isproper: # Proper rotation
    #        t = 3 # try angle != 180
    #        #det180 = get_sym_det(rot + self._E3D)
    #        if (self + identity).det == 0: t = 4 # 180 rotation
    #    else:
    #        # Mirror symmetry or Improper rotation
    #        t = 6
    #        #detmirror = get_sym_det(rot - self._E3D)
    #        if (self - identity).det == 0:
    #            t = 5 # Mirror symmetry if an eigenvalue is 1

    #    return t


# TODO: Need to find an easy way to map classes in internal database
# onto classes computed by client code when calculation has been done
# with non-conventional settings (spglib?)
class Irrep(object):
    """
    This object represents an irreducible representation.

    .. attributes::

        traces: all_traces[nsym]. The trace of each irrep.
        character: character[num_classes]
    """
    def __init__(self, name, dim, mats, class_range):
        """
        Args:
            name:  Name of the irreducible representation.
            dim: Dimension of the irreducible representation.
            mats: Array of shape [nsym,dim,dim] with the irreducible
                representations of the group. mats are packed in classes.
            class_range: List of tuples, each tuple gives the start and stop index for the class.
                e.g. [(0, 2), (2,4), (4,n)]
        """
        self.name = name
        self._mats = np.reshape(np.array(mats), (-1, dim, dim))

        self.traces = [m.trace() for m in self.mats]

        self.class_range = class_range
        self.nclass = len(class_range)

        # Compute character table.
        character = self.nclass * [None]
        for icls, (start, stop) in enumerate(self.class_range):
            t0 = self.traces[start]
            isok = all(t0 == self.traces[i] for i in range(start, stop))
            character[icls] = t0

        self._character = character

    @property
    def mats(self):
        return self._mats

    @property
    def character(self):
        return self._character

    #@lazy_property
    #def dataframe(self):


def bilbao_ptgroup(sch_symbol):
    """
    Returns an instance of :class:`BilbaoPointGroup` from a string with the point group symbol
    or a number with the spacegroup ID.
    """
    sch_symbol = any2sch(sch_symbol)

    from abipy.core.irrepsdb import _PTG_IRREPS_DB
    entry = _PTG_IRREPS_DB[sch_symbol].copy()
    entry.pop("nclass")
    entry["sch_symbol"] = sch_symbol

    return BilbaoPointGroup(**entry)


class BilbaoPointGroup(object):
    """
    A :class:`BilbaoPointGroup` is a :class:`Pointgroup` with irreducible representations
    """
    def __init__(self, sch_symbol, rotations, class_names, class_range, irreps):
        # Rotations are grouped in classes.
        self.sch_symbol = sch_symbol
        self.rotations = np.reshape(rotations, (-1, 3, 3))
        self.class_names = class_names
        self.nclass = len(class_names)

        # List of tuples, each tuple gives the start and stop index for the class.
        # e.g. [(0, 2), (2,4), (4,n)]
        self.class_range = class_range
        self.class_len = [stop - start for start, stop in class_range]

        # The number of irreps must equal the number of classes.
        assert len(irreps) == self.nclass
        self.irreps, self.irreps_by_name = [], {}
        for name, d in irreps.items():
            mats = d["matrices"]
            assert len(mats) == self.num_rots
            irrep = Irrep(name, d["dim"], mats, class_range=self.class_range)
            self.irreps.append(irrep)
            self.irreps_by_name[name] = irrep

    @property
    def herm_symbol(self):
        """Hermann-Mauguin symbol."""
        return herm2sch(self.sch_symbol)

    @property
    def spgid(self):
        """ID in the space group table."""
        return sch2spgid(self.sch_symbol)

    @property
    def num_rots(self):
        """Number of rotations."""
        return len(self.rotations)

    @property
    def num_irreps(self):
        """Number of irreducible representations."""
        return len(self.irreps)

    @property
    def irrep_names(self):
        """List with the names of the irreps."""
        return list(self.irreps_by_name.keys())

    @lazy_property
    def character_table(self):
        """
        Dataframe with irreps.
        """
        # Caveat: class names are not necessarly unique --> use np.stack
        import pandas as pd
        name_mult = [name + " [" + str(mult) + "]" for (name, mult) in zip(self.class_names, self.class_len)]
        columns = ["name"] + name_mult

        stack = np.stack([irrep.character for irrep in self.irreps])
        index = [irrep.name for irrep in self.irreps]
        df = pd.DataFrame(stack, columns=name_mult, index=index)
        df.index.name = "Irrep"
        df.columns.name = self.sch_symbol

	# TODO
        #print(df)
        # Convert complex --> real if all entries in a colums are real.
        #for k in name_mult:
        #    if np.all(np.isreal(df[k].values)):
        #        #df[k] = df[k].values.real
        #        df[k] = df[k].astype(float)

        return df

    def to_string(self, verbose=0):
        """
        Return string with the character_table
        """
        return self.character_table.to_string()

    #def decompose(self, character):
    #   od = collections.OrderedDict()
    #   for irrep in self.irreps:
    #       irrep.name
    #       irrep.character
    #   return od

    #def show_irrep(self, irrep_name):
    #    """Show the mapping rotation --> irrep mat."""
    #    irrep = self.irreps_by_name[irrep_name]

    #def irrep_from_character(self, character, rotations, tol=None):
    #    """
    #    Main entry point for client code.
    #    This routine receives a character computed from the user and finds the
    #    irreducible representation.
    #    """

    #def map_rotclasses(self, rotations_in_classes)
    #def map_rotation(self, rotations_in_classes)

    def auto_test(self):
        """
        Perform internal consistency check. Return 0 if success
        """
        #print("rotations\n", self.rotations)
        rot_group = LatticePointGroup(self.rotations)
        if not rot_group.is_group():
            print("rotations do not form a group!")
            return 1

        # Symmetries should be ordered in classes.
        # Here we recompute the classes by calling rot_group.class_indices.
        # We then sort the indices and we compare the results with the ref data stored in the Bilbao database.
        calc_class_inds = [sorted(l) for l in rot_group.class_indices]
        #print(calc_class_inds)
        assert len(calc_class_inds) == len(self.class_range)

        for calc_inds, ref_range in zip(calc_class_inds, self.class_range):
            ref_inds = list(range(ref_range[0], ref_range[1]))
            if calc_inds != ref_inds:
                print("Rotations are not ordered in classes.", calc_inds, ref_inds)
                return 2

        # Do we have a representation of the Group?
        mult_table = rot_group.mult_table
        max_err = 0.0

        for idx1, rot1 in enumerate(rot_group):
            for idx2, rot2 in enumerate(rot_group):
                idx_prod = mult_table[idx1, idx2]
                for irrep in self.irreps:
                    mat_prod = np.dot(irrep.mats[idx1], irrep.mats[idx2])
                    err = (mat_prod - irrep.mats[idx_prod]).max()
                    max_err = max(max_err, abs(err))

        if max_err > 1e-5:
            print("Irreps do not form a representation of the group, max_err: ", max_err)
            return 3

        # TODO
        # Test orthogonality theorem

        # Test the orthogonality relation of traces.
        max_err = 0.0
        for (ii, jj), (irp1, irp2) in iuptri(self.irreps, with_inds=True):
            trac1, trac2 = irp1.traces, irp2.traces
            err = np.vdot(trac1, trac2) / self.num_rots
            if ii == jj: err -= 1.0
            max_err = max(max_err, abs(err))

        if max_err > 1e-5:
            print("Error in orthogonality relation of traces: ", max_err)
            return 4

        # Success.
        return 0


# Schoenflies, Hermann-Mauguin, spgid
_PTG_IDS = [
    ("C1" , "1",     1),
    ("Ci" , "-1",    2),
    ("C2" , "2",     3),
    ("Cs" , "m",     6),
    ("C2h", "2/m",   10),
    ("D2" , "222",   16),
    ("C2v", "mm2",   25),
    ("D2h", "mmm",   47),
    ("C4" , "4",     75),
    ("S4" , "-4",    81),
    ("C4h", "4/m",   83),
    ("D4" , "422",   89),
    ("C4v", "4mm",   99),
    ("D2d", "-42m",  111),
    ("D4h", "4/mmm", 123),
    ("C3" , "3",     143),
    ("C3i", "-3",    147),
    ("D3" , "32",    149),
    ("C3v", "3m",    156),
    ("D3d", "-3m",   162),
    ("C6" , "6",     168),
    ("C3h", "-6",    174),
    ("C6h", "6/m",   175),
    ("D6" , "622",   177),
    ("C6v", "6mm",   183),
    ("D3h", "-6m2",  189),
    ("D6h", "6/mmm", 191),
    ("T"  , "23",    195),
    ("Th" , "m-3",   200),
    ("O"  , "432",   207),
    ("Td" , "-43m",  215),
    ("Oh" , "m-3m",  221),
]

_SCH2HERM = {t[0]: t[1] for t in _PTG_IDS}
_HERM2SCH = {t[1]: t[0] for t in _PTG_IDS}
_SPGID2SCH = {t[2]: t[0] for t in _PTG_IDS}
_SCH2SPGID = {t[0]: t[2] for t in _PTG_IDS}

sch_symbols = list(_SCH2HERM.keys())


def sch2herm(sch_symbol):
    """Convert from Schoenflies to Hermann-Mauguin."""
    return _SCH2HERM.get(sch_symbol, None)


def sch2spgid(sch_symbol):
    """Convert from Schoenflies to the space group id."""
    return _SCH2SPGID.get(sch_symbol, None)


def herm2sch(herm_symbol):
    """Convert from Hermann-Mauguin to Schoenflies."""
    return _HERM2SCH.get(herm_symbol, None)


def spgid2sch(spgid):
    """Return the Schoenflies symbol from the space group identifier."""
    return _SPGID2SCH.get(spgid, None)


def any2sch(obj):
    """Convert string or int to Schoenflies symbol. Returns None if invalid input"""
    if is_string(obj):
        if obj in sch_symbols:
            return obj
        else:
            # Try Hermann-Mauguin
            return herm2sch(obj)
    else:
        # Spacegroup ID?
        return spgid2sch(obj)
