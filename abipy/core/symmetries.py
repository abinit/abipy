"""Objects used to deal with symmetry operations in crystals."""
from __future__ import division, print_function

import numpy as np
import collections

from abipy.core.kpoints import wrap_to_ws, issamek
from abipy.iotools import as_etsfreader


__all__ = [
    "SymmOp",
    "SpaceGroup",
]

def wrap_in_ucell(x):
    """
    Transforms x in its corresponding reduced number in the interval [0,1[."
    """
    return x % 1

def isinteger(x, atol=1e-08):
    """
    True if all x is integer within the absolute tolerance atol.

    >>> isinteger([1., 2.])
    True
    >>> isinteger(1.01, atol=0.011)
    True
    >>> isinteger([1.01, 2])
    False
    """
    int_x = np.around(x)
    return np.allclose(int_x, x, atol=atol)


def mati3inv(mm, trans=True):
    """
    Invert and transpose orthogonal 3x3 matrix of INTEGER elements.

    Args:
        mm:
            3x3 matrix-like object with integer elements

    Returns:
        ndarray with the TRANSPOSE of the inverse of mm if trans==True.
        If trans==False, the inverse of mm is returned.

    .. note::

       Used for symmetry operations. This function applies to *ORTHOGONAL* matrices only.
       Since these form a group, inverses are also integer arrays.
    """
    mm = np.array(mm)
    assert mm.dtype in [np.int, np.int8, np.int16, np.int32, np.int64]

    mit = np.zeros((3,3), dtype=np.int)
    mit[0,0] = mm[1,1] * mm[2,2] - mm[2,1] * mm[1,2]
    mit[1,0] = mm[2,1] * mm[0,2] - mm[0,1] * mm[2,2]
    mit[2,0] = mm[0,1] * mm[1,2] - mm[1,1] * mm[0,2]
    mit[0,1] = mm[2,0] * mm[1,2] - mm[1,0] * mm[2,2]
    mit[1,1] = mm[0,0] * mm[2,2] - mm[2,0] * mm[0,2]
    mit[2,1] = mm[1,0] * mm[0,2] - mm[0,0] * mm[1,2]
    mit[0,2] = mm[1,0] * mm[2,1] - mm[2,0] * mm[1,1]
    mit[1,2] = mm[2,0] * mm[0,1] - mm[0,0] * mm[2,1]
    mit[2,2] = mm[0,0] * mm[1,1] - mm[1,0] * mm[0,1]

    dd = mm[0,0] * mit[0,0] + mm[1,0] * mit[1,0] + mm[2,0] * mit[2,0]

    # Make sure matrix is not singular
    if dd == 0:
        raise ValueError("Attempting to invert integer array: %s\n ==> determinant is zero." % mm)

    mit = mit // dd
    if trans:
        return mit
    else:
        return mit.T


def _get_det(mat):
    """
    Return the determinant of a 3x3 rotation matrix mat.

    raises:
        ValueError if det not in [+1,-1]
    """
    det = mat[0,0]* (mat[1,1]*mat[2,2] - mat[1,2]*mat[2,1])\
        - mat[0,1]* (mat[1,0]*mat[2,2] - mat[1,2]*mat[2,0])\
        + mat[0,2]* (mat[1,0]*mat[2,1] - mat[1,1]*mat[2,0])

    if abs(det) != 1:
        raise ValueError("determinant must be \pm 1 while it is %s" % det)
    else:
        return det


class SymmOp(object):
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

    def __init__(self, rot_r, tau, time_sign, afm_sign, rot_g=None):
        """
        This object represents a space group symmetry i.e. a symmetry of the crystal.

        Args:
            rot_r:
                3x3 integer matrix with the rotational part in real space in reduced coordinates (C order).
            tau:
                fractional translation in reduced coordinates.
            time_sign:
                -1 if time reversal can be used, otherwise +1.
            afm_sign:
                anti-ferromagnetic part [+1,-1].
        """
        rot_r = np.asarray(rot_r)

        # Store R and R^{-1} in real space.
        self.rot_r, self.rotm1_r = rot_r, mati3inv(rot_r, trans=False) 
        self.tau = np.asarray(tau)

        self.afm_sign, self.time_sign = afm_sign, time_sign
        assert afm_sign in  [-1, 1] and time_sign in [-1, 1]

        # Compute symmetry matrix in reciprocal space: S = R^{-1t}
        if rot_g is None:
            self.rot_g = mati3inv(rot_r, trans=True)
        else:
            assert np.all(rot_g == mati3inv(rot_r, trans=True))
            self.rot_g = rot_g

    @staticmethod
    def _vec2str(vec):
        return "%2d,%2d,%2d" % tuple(v for v in vec)

    def __str__(self):
        s = ""
        for i in range(3):
            s +=  "[" + self._vec2str(self.rot_r[i]) + ", %.3f]  " % self.tau[i] + "[" + self._vec2str(self.rot_g[i]) + "] "
            if i == 0:
                s += " time_sign=%2d, afm_sign=%2d, det=%2d" % (self.time_sign, self.afm_sign, self.det)
            s += "\n"

        return s

    def __eq__(self, other):
        # Note the two fractional traslations are equivalent if they differ by a lattice vector.
        return (np.all(self.rot_r == other.rot_r) and
                isinteger(self.tau-other.tau, atol=self._ATOL_TAU) and
                self.afm_sign == other.afm_sign and
                self.time_sign == other.time_sign
                )

    def __ne__(self, other):
        return not self == other

    def __mul__(self, other):
        """
        Returns a new SymmOp which is equivalent to apply  the "other" `SymmOp` followed by this one.
        {R,t} {S,u} = {RS, Ru + t}
        """
        return SymmOp(rot_r=np.dot(self.rot_r, other.rot_r),
                      tau=self.tau + np.dot(self.rot_r, other.tau),
                      time_sign=self.time_sign * other.time_sign,
                      afm_sign=self.afm_sign * other.afm_sign
                      )

    def __hash__(self):
        return 8 * self.trace + 4 * self.det + 2 * self.time_sign

    @property
    def is_identity(self):
        return (np.all(self.rot_r == np.eye(3, dtype=np.int)) and
                isinteger(self.tau, atol=self._ATOL_TAU) and
                self.time_sign == 1 and
                self.afm_sign == 1
               )

    def inverse(self):
        """
        Returns inverse of transformation i.e. {R^{-1}, -R^{-1} tau}.
        """
        return SymmOp(rot_r=self.rotm1_r,
                      tau=-np.dot(self.rotm1_r, self.tau),
                      time_sign=-self.time_sign,
                      afm_sign=self.afm_sign
                      )

    def conjugate(self, other):
        """Returns X^-1 S X where X is the other symmetry operation."""
        return other.inverse() * self * other

    @property
    def det(self):
        """Determinant of the rotation matrix [-1, +1]."""
        try:
            return self._det

        except AttributeError:
            self._det = _get_det(self.rot_r)
            return self._det

    @property
    def trace(self):
        """Trace of the rotation matrix."""
        try:
            return self._trace

        except AttributeError:
            self._trace = self.rot_r.trace()
            return self._trace

    @property
    def is_proper(self):
        """True if the rotational part has determinant == 1."""
        return self.det == +1

    @property
    def has_timerev(self):
        """True if symmetry contains the time-reversal operator."""
        return self.time_sign == -1

    @property
    def is_fm(self):
        """True if self if ferromagnetic symmetry."""
        return self.afm_sign == +1

    @property
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

    #def preserve_k(self, frac_coords):
    #    """
    #    Check if the operation preserves the k-point modulo a reciprocal lattice vector.

    #    Returns:
    #        bool, Sk - k
    #    """
    #    sk = self.rotate_k(frac_coords, wraps_tows=False)

    #    return issamek(sk, frac_coords), sk - frac_coords

    def rotate_r(self, frac_coords, in_ucell=False):
        """
        Apply the symmetry operation to point in real space given in reduced coordinates.

        .. note:
            We use the convention: symmop(r) = R^{-1] (r -tau)
        """
        rotm1_rmt = np.dot(self.rotm1_r, frac_coords - self.tau)

        return wrap_in_ucell(rotm1_rmt) if in_ucell else rotm1_rmt

    def rotate_gvecs(self, gvecs):
        """
        Apply the symmetry operation to the list of gvectors gvecs in reduced coordinates.
        """
        rot_gvecs = np.zeros_like(gvecs)
        for ig, gvec in enumerate(gvecs):
            rot_gvecs[ig] = np.dot(self.rot_g, gvec) * self.time_sign

        return rot_gvecs

#########################################################################################


class SymmOpList(collections.Sequence):

    def __init__(self, symrel, tnons, symafm, has_timerev, inord="C"):
        """
        Args:
            symrel:
                (nsym,3,3) array with the rotational part of the symmetries in real
                space (reduced coordinates are assumed, see also `inord` for the order.
            tnons:
                (nsym,3) array with fractional translation in reduced coordinates.
            symafm:
                (nsym) array with +1 for Ferromagnetic symmetry and -1 for AFM
            has_timerev:
                True if time-reversal symmetry is included.
            inord:
                storage order of mat in symrel[:]. If inord == "F", mat.T is stored
                as matrices are always stored in C-order. Use inord == "F" if you have 
                read symrel from an external file produced by abinit.

        .. note:
            All the arrays are store in C-order. Use as_fortran_arrays to extract data that
            can be passes to Fortran routines.
        """
        inord = inord.upper()
        assert inord in ["C", "F"]

        # Time reversal symmetry.
        self._has_timerev = has_timerev
        self._time_signs = [+1, -1] if self.has_timerev else [+1]

        self._symrel, self._tnons, self._symafm = map(np.asarray, (symrel, tnons, symafm))

        if len(self.symrel) != len(self.tnons) or len(self.symrel) != len(self.symafm):
            raise ValueError("symrel, tnons and symafm must have equal shape[0]")

        if inord == "F": # Fortran to C.
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
                                       rot_g=self.symrec[isym],
                                ))

        self._symmops = tuple(all_syms)

    def __len__(self):
        return len(self._symmops)

    def __iter__(self):
        return self._symmops.__iter__()
                                      
    def __getitem__(self, slice):
        return self._symmops[slice]

    def __contains__(self, symmop):
        return symmop in self._symmops

    def count(self, symmop):
        """Returns the number of occurences of symmop in self."""
        return self._symmops.count(symmop)

    def index(self, symmop):
        """Return the (first) index of symmop in self. Raises ValueError if not found.""" 
        return self._symmops.index(symmop)

    def find(self, symmop):
        """Return the (first) index of symmop in self. -1 if not found.""" 
        try:
            return self.index(symmop)
        except ValueError:
            return -1

    @property
    def has_timerev(self):
        """True if time-reversal symmetry is present."""
        return self._has_timerev

    @property
    def symrel(self):
        return self._symrel

    @property
    def tnons(self):
        return self._tnons

    @property
    def symrec(self):
        return self._symrec

    @property
    def symafm(self):
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
            time_sign:
                If specified, only symmetryes with time-reversal sign time_sign are returned.
            afm_sign:
                If specified, only symmetryes with anti-ferromagnetic part afm_sign are returned.

        returns:
            tuple of `SymmOp` instances.
        """
        symmops = []
        for sym in self._symmops:
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

    def is_group(self):
        """Returns True if self is a group."""
        check = 0
    
        # Identity must be present.
        if [op.is_identity for op in self].count(True) != 1:
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

    def asdict(self):
        """
        Returns a dictions where the keys are the symmetry operations and 
        the values are the indices of the operations in the iterable.
        """
        d = {op: idx for (idx, op) in enumerate(self)}
        assert len(d) == len(self)
        return d

    def to_fortran_arrays(self):
        fort_arrays = collections.namedtuple("FortranSpaceGroupArrays", "symrel symrec tnons symafm timrev")

        symrel = np.asfortranarray(self.symrel.T)
        symrec = np.asfortranarray(self.symrec.T)

        for isym in range(self.num_spatial_symmetries):
            symrel[:,:,isym] = symrel[:,:,isym].T
            symrec[:,:,isym] = symrec[:,:,isym].T

        return fort_arrays(
            symrel=symrel,
            symrec=symrec,
            tnons =np.asfortranarray(self.tnons.T),
            symafm=self.symafm,
            timrev = 2 if self.has_timerev else 1
        )

    #def mult_table(self):
    #    """
    #    Given a set of nsym 3x3 operations which are supposed to form a group, 
    #    this routine constructs the multiplication table of the group.
    #    mtab(nsym,nsym)=The index of the product S_i * S_j in the input set sym.
    #    """
    #    d = self.asdict()
    #    for op1 in self:
    #        for op2 in self:
    #            op12 = op1 * op2 
    #            table = d[op12]
    #    return table

    #@property
    #def classes(self):
    #    """
    #    A class is defined as the set of distinct elements obtained by 
    #    considering for each element, S, of the group all its conjugate
    #    elements X^-1 S X where X range over all the elements of the group.
    #    """

    #    try:
    #        return self._classes

    #    except AttributeError:
    #        
    #        num_classes, found, classes = -1, len(self) * [False], len(self) * [None]

    #        for (ii, op1) in enumerate(self):
    #            if found[ii]: continue
    #            num_classes += 1

    #            for (jj, op2) in enumerate(self):
    #                # Form conjugate
    #                op1_conj = op1.conjugate(op2)
    #                for (kk, op3) in enumerate(self):
    #                    if not found[kk] and op1_conj == op3:
    #                        found[kk] = True
    #                        classes[num_classes]


class SpaceGroup(SymmOpList):
    """Container storing the space group symmetries."""

    def __init__(self, spgid, symrel, tnons, symafm, has_timerev, inord="C"):
        """
        Args:
            spgid:
                space group number (from 1 to 232, 0 if cannot be specified).
            symrel:
                (nsym,3,3) array with the rotational part of the symmetries in real
                space (reduced coordinates are assumed, see also `inord` for the order.
            tnons:
                (nsym,3) array with fractional translation in reduced coordinates.
            symafm:
                (nsym) array with +1 for Ferromagnetic symmetry and -1 for AFM
            has_timerev:
                True if time-reversal symmetry is included.
            inord:
                storage order of mat in symrel[:]. If inord == "F", mat.T is stored
                as matrices are always stored in C-order. Use inord == "F" if you have 
                read symrel from an external file produced by abinit.

        .. note:
            All the arrays are store in C-order. Use as_fortran_arrays to extract data that
            can be passes to Fortran routines.
        """
        SymmOpList.__init__(self, symrel, tnons, symafm, has_timerev, inord=inord)

        self.spgid = spgid
        assert self.spgid in range(0, 233)


    @classmethod
    def from_file(cls, file, inord="F"):
        """Initialize the object from a Netcdf file."""
        file, closeit = as_etsfreader(file)

        new = cls(spgid=file.read_value("space_group"),
                  symrel=file.read_value("reduced_symmetry_matrices"),
                  tnons=file.read_value("reduced_symmetry_translations"),
                  symafm=file.read_value("symafm"),
                  has_timerev=True,  # FIXME not treated by ETSF-IO.
                  inord=inord,
                  )
        if closeit:
            file.close()

        return new

    def __str__(self):
        """String representation."""
        lines = ["spgid %d, num_spatial_symmetries %d, has_timerev %s" % (
            self.spgid, self.num_spatial_symmetries, self.has_timerev)]
        app = lines.append

        #app(["rot_r  rot_g"]
        for op in self.symmops(time_sign=+1):
            app(str(op))

        return "\n".join(lines)


#class LittleGroup(SymmOpList):

class Irrep(object):
    """
    .. attributes:

        name: 
            Name of the irreducible representation.

        dim:
            Dimension of the irreducible representation.

        nsym:
            Number of symmetries.

        mats:
            array of shape [nsym,dim,dim] with
            the irreducible representations of the group.

        traces:
            traces[nsym]. The trace of each matrix.
    """

    def __init__(self, name, mats):
        self.name = name
        assert len(mats.shape) == 3
        self.mats = mats
        self.nsym = len(mats)
        self.dim = mats.shape[1]
        assert self.dim == mats.shape[2]

        self.trace = tuple([m.trace() for m in mats])


#class IrrepsDatabase(dict)
#    _PTGROUP_NAMES = [
#      "1",   
#      "-1",
#      "2",
#      "m",
#      "2/m",
#      "222",
#      "mm2",
#      "mmm",
#      "4",
#      "-4",
#      "4/m",
#      "422",
#      "4mm",
#      "-42m",
#      "4/mmm",
#      "3",
#      "-3",
#      "32",
#      "3m",
#      "-3m",
#      "6",
#      "-6",
#      "6/m", 
#      "622",
#      "6mm",
#      "-62m",
#      "6/mmm",
#      "23",
#      "m-3",
#      "432",
#      "-43m ",
#      "m-3m",
#    ]
