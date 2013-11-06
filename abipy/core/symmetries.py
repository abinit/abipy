"""Objects used to deal with symmetry operations in crystals."""
from __future__ import division, print_function

import numpy as np
import warnings
import collections

from abipy.core.kpoints import wrap_to_ws, issamek
from abipy.iotools import as_etsfreader


__all__ = [
    "LatticeRotation",
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


def mati3inv(mat3, trans=True):
    """
    Invert and transpose orthogonal 3x3 matrix of INTEGER elements.

    Args:
        mat3:
            3x3 matrix-like object with integer elements

    Returns:
        ndarray with the TRANSPOSE of the inverse of mat3 if trans==True.
        If trans==False, the inverse of mat3 is returned.

    .. note::

       Used for symmetry operations. This function applies to *ORTHOGONAL* matrices only.
       Since these form a group, inverses are also integer arrays.
    """
    mat3 = np.array(mat3)
    assert mat3.dtype in [np.int, np.int8, np.int16, np.int32, np.int64]

    mit = np.empty((3,3), dtype=np.int)
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
        raise ValueError("Attempting to invert integer array: %s\n ==> determinant is zero." % mat3)

    mit = mit // dd
    if trans:
        return mit
    else:
        return mit.T


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
        raise ValueError("determinant must be \pm 1 while it is %s" % det)

    return det


class SymmOp(object):
    """Crystalline symmetry."""
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

    def __repr__(self):
        return str(self)

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
        return not (self == other)

    def __mul__(self, other):
        """
        Returns a new SymmOp which is equivalent to apply  the "other" `SymmOp` followed by this one.
        {R,t} {S,u} = {RS, Ru + t}
        """
        return SymmOp(rot_r=np.dot(self.rot_r, other.rot_r),
                      tau=self.tau + np.dot(self.rot_r, other.tau),
                      time_sign=self.time_sign * other.time_sign,
                      afm_sign=self.afm_sign * other.afm_sign)

    def __hash__(self):
        """
        `Symmop` can be used as keys in dictionaries.
        Note that the hash is computed from integer values. 
        """
        return 8 * self.trace + 4 * self.det + 2 * self.time_sign

    @property
    def is_symmorphic(self):
        """True if the fractional translation is non-zero."""
        return np.any(np.abs(self.tau) > 0.0)

    @property
    def is_identity(self):
        """True is self is the identity operator."""
        return (np.all(self.rot_r == np.eye(3, dtype=np.int)) and
                isinteger(self.tau, atol=self._ATOL_TAU) and
                self.time_sign == 1 and
                self.afm_sign == 1
               )

    def inverse(self):
        """Returns inverse of transformation i.e. {R^{-1}, -R^{-1} tau}."""
        return SymmOp(rot_r=self.rotm1_r,
                      tau=-np.dot(self.rotm1_r, self.tau),
                      time_sign=self.time_sign,
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

    def preserve_k(self, frac_coords, ret_g0=True):
        """
        Check if the operation preserves the k-point modulo a reciprocal lattice vector.

        Args:
            frac_coords:
                Fractional coordinates of the k-point
            ret_g0:
                False if only the boolean result is wanted.

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

        .. note:

            We use the convention: symmop(r) = R^{-1] (r - tau)
        """
        rotm1_rmt = np.dot(self.rotm1_r, frac_coords - self.tau)

        return wrap_in_ucell(rotm1_rmt) if in_ucell else rotm1_rmt

    def rotate_gvecs(self, gvecs):
        """
        Apply the symmetry operation to the list of gvectors gvecs in reduced coordinates.

        Args:
            gvecs:
                ndarray with shape [ng, 3] containing the reduced coordinates of the G-vectors.
        Returns:
            rot_gvecs:
                ndarray with shape [ng, 3] containing the result of self(G).
        """
        rot_gvecs = np.empty_like(gvecs)

        for ig, gvec in enumerate(gvecs):
            rot_gvecs[ig] = np.dot(self.rot_g, gvec) * self.time_sign

        return rot_gvecs


class OpList(collections.Sequence):

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

        .. warning: 

                The order of the operations in self and in other 
                is not relevant.
        """
        if other is None: return False
        if len(self) != len(other): return False

        # Check if each operation in self is also present 
        # in other. The order is irrelevant.
        founds = []
        for i, op1 in enumerate(self):
            if op1 not in other: return False
            founds.append(i)

        if len(set(founds)) == len(founds):
            return True
        
        warnings.warn("self contains duplicated symops! Likely a bug!")
        return False

    def __ne__(self, other):
        return not (self == other)

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
        Returns a dictionary where the keys are the symmetry operations and 
        the values are the indices of the operations in the iterable.
        """
        d = {op: idx for (idx, op) in enumerate(self)}
        assert len(d) == len(self)
        return d

    @property
    def mult_table(self):
        """
        Given a set of nsym 3x3 operations which are supposed to form a group, 
        this routine constructs the multiplication table of the group.
        mtable[i,j] gives the index of the product S_i * S_j.
        """
        try:
            return self._mult_table

        except AttributeError:
            mtable = np.empty((len(self), len(self)), dtype=np.int)
    
            d = self.asdict()
            for (i, op1) in enumerate(self):
                for (j, op2) in enumerate(self):
                    op12 = op1 * op2 
                    # Save the index of op12 in self
                    mtable[i,j] = d[op12]

            self._mult_table = mtable
            return self._mult_table

    @property
    def classes(self):
        """
        A class is defined as the set of distinct elements obtained by 
        considering for each element, S, of the group all its conjugate
        elements X^-1 S X where X ranges over all the elements of the group.

        Returns:
            Nested list l = [cls0_indices, cls1_indeces, ...] wheree each sublist 
            contains the indices of the class. len(l) equals the number of classes.
        """
        try:
            return self._classes

        except AttributeError:
            
            found, classes = len(self) * [False], [[] for i in range(len(self))]

            num_classes = -1
            for (ii, op1) in enumerate(self):
                if found[ii]: continue
                num_classes += 1

                for (jj, op2) in enumerate(self):
                    # Form conjugate and search it among the operations 
                    # that have not been found yet.
                    op1_conj = op1.conjugate(op2)

                    for (kk, op3) in enumerate(self):
                        if not found[kk] and op1_conj == op3:
                            found[kk] = True
                            classes[num_classes].append(kk)

            self._classes = classes[:num_classes+1]

            #assert sum(len(c) for c in self._classes) == len(self)
            return self._classes


class SpaceGroup(OpList):
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
            All the arrays are store in C-order. Use as_fortran_arrays to extract data 
            that can be passes to Fortran routines.
        """
        self.spgid = spgid
        assert self.spgid in range(0, 233)

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

        self._ops = tuple(all_syms)

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

    def __repr__(self):
        return str(self)

    def __str__(self):
        """String representation."""
        lines = ["spgid %d, num_spatial_symmetries %d, has_timerev %s" % (
            self.spgid, self.num_spatial_symmetries, self.has_timerev)]
        app = lines.append

        for op in self.symmops(time_sign=+1):
            app(str(op))

        return "\n".join(lines)

    @property
    def is_symmorphic(self):
        """True if there's at least one operation with non-zero fractional translation."""
        for op in self:
            if op.is_symmorphic:
                return True

        return False

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
                If specified, only symmetries with time-reversal sign time_sign are returned.
            afm_sign:
                If specified, only symmetries with anti-ferromagnetic part afm_sign are returned.

        returns:
            tuple of `SymmOp` instances.
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

    #def to_array(self):
    #    fort_arrays = collections.namedtuple("FortranSpaceGroupArrays", "symrel symrec tnons symafm timrev")
                                                                                                              
    #    symrel = np.asfortranarray(self.symrel.T)
    #    symrec = np.asfortranarray(self.symrec.T)
                                                                                                              
    #    for isym in range(self.num_spatial_symmetries):
    #        symrel[:,:,isym] = symrel[:,:,isym].T
    #        symrec[:,:,isym] = symrec[:,:,isym].T
                                                                                                              
    #    return fort_arrays(
    #        symrel=symrel,
    #        symrec=symrec,
    #        tnons =np.asfortranarray(self.tnons.T),
    #        symafm=self.symafm,
    #        timrev = 2 if self.has_timerev else 1
    #    )



    def get_little_group(self, kpoint):
        """
        Find the little group of the kpoint

        Args:
            kpoint:
                Accept vector with the reduced coordinates or `Kpoint` object.
    
        Returns:
            ltg_symmops, g0vecs, indices

            ltg_symmops is a tuple with the symmetry operations that preserve the k-point i.e. Sk = k + G0
            g0vecs is the tuple for G0 vectors for each operation in ltg_symmops
            indices gives the index of the little group operation in the initial spacegroup.
        """
        frac_coords = getattr(kpoint, "frac_coords", kpoint)

        isyms, g0vecs = [], []

        # Exclude AFM operations.
        for isym, symmop in enumerate(self.fm_symmops):
            is_same, g0 = symmop.preserve_k(frac_coords)
            if is_same:
                isyms.append(isym)
                g0vecs.append(g0)

        ltg_symmops = [self[isym] for isym in isyms]

        return ltg_symmops, g0vecs, isyms


_E3D = np.identity(3,  np.int)

#class LatticeRotation(np.ndarray):
class LatticeRotation(object):
    """
    This object defines a pure rotation of the lattice (proper, improper, mirror symmetry)
    that is a rotation which is compatible with a lattice. The rotation matrix is
    expressed in reduced coordinates, therefore its elements are integers.
    """
    #def __new__(cls, input_array, unit, unit_type=None):
    #    # Input array is an already formed ndarray instance
    #    # We first cast to be our class type
    #    obj = np.asarray(input_array).view(cls)
    #    # add the new attributes to the created instance
    #    obj._unit = Unit(unit)
    #    obj._unit_type = unit_type
    #    return obj

    #def __array_finalize__(self, obj):
    #    """
    #    See http://docs.scipy.org/doc/numpy/user/basics.subclassing.html for
    #    comments.
    #    """
    #    if obj is None:
    #        return
    #    self._unit = getattr(obj, "_unit", None)
    #    self._unit_type = getattr(obj, "_unit_type", None)

    #def __reduce__(self):
    #    #print("in reduce")
    #    reduce = list(super(ArrayWithUnit, self).__reduce__())
    #    #print("unit",self._unit)
    #    #print(reduce[2])
    #    reduce[2] = {"np_state": reduce[2], "_unit": self._unit}
    #    return tuple(reduce)

    #def __setstate__(self, state):
    #    #print("in setstate %s" % str(state))
    #    super(ArrayWithUnit, self).__setstate__(state["np_state"])
    #    self._unit = state["_unit"]

    def __init__(self, rotation):
        self.rotation = np.matrix(rotation, np.int)
        self.rotation.shape = (3,3)

    #def __repr__(self):
    #    return str(self)

    #def __str__(self):
    #    lines =  "Rotation: " + str(self.order) + ", versor: " + str(self.versor) + ", " + str(self.trcoords) + "\n"
    #    lines.append(str(self.rotation))
    #    return "\n".join(lines)

    # Might subclass np.matrix though.
    def __eq__(self, other):
        return np.allclose(self.rotation, other.rotation)

    def __ne__(self, other):
        return not (self == other)

    # Implement the unary arithmetic operations (+, -)
    def __pos__(self): 
        return self

    def __neg__(self): 
        return self.__class__(-self.rotation) 

    def __mul__(self, other):
        return self.__class__(self.rotation * other.rotation)

    def __pow__(self, intexp, modulo=1):
       if intexp ==  0: return self.__class__(_E3D)
       if intexp  >  0: return self.__class__(self.rotation**intexp)
       if intexp == -1: return self.inverse()
       if intexp  <  0: return self.__pow__(-intexp).inverse()

    @property
    def det(self):
        """Return the determinant of a symmetry matrix mat[3,3]. It must be +-1"""
        try:
            return self._det
        except AttributeError:
            self._det = _get_det(self.rotation)
            return self._det

    @property
    def trace(self):
        """The trace of the rotation"""
        return self.rotation.trace()[0,0]

    @property
    def is_proper(self):
        """True if proper rotation"""
        return self.det == 1

    def inverse(self): 
        """
        Invert an orthogonal 3x3 matrix of INTEGER elements.
        Note use of integer arithmetic. Raise ValueError if not invertible.
        """
        inv = mati3inv(self.rotation, trans=False)
        return self.__class__(inv)

    @property
    def rottype(self):
        """
        Receive a 3x3 orthogonal matrix and reports its type:
            1 Identity
            2 Inversion
            3 Proper rotation of an angle <> 180 degrees
            4 Proper rotation of 180 degrees
            5 Mirror symmetry
            6 Improper rotation
        """
        #rot = self.rotation # Just an alias. 
        # Treat identity and inversion first
        #identity  = Rotation(_E3D)

        if self.isE: return 1
        if self.isI: return 2
    
        if self.isproper: # Proper rotation
            t = 3 # try angle != 180
            #det180 = get_sym_det(rot+_E3D)
            if (self + identity).det == 0: t = 4 # 180 rotation
        else: 
            # Mirror symmetry or Improper rotation
            t = 6
            #detmirror = get_sym_det(rot-_E3D)
            if (self - identity).det == 0: 
                t = 5 # Mirror symmetry if an eigenvalue is 1

        return t

    @property
    def isE(self):
        """True if it is the identity"""
        return np.allclose(self.rotation, _E3D)

    @property
    def isI(self):
        """True if it is the inversion"""
        return np.allclose(self.rotation, -_E3D)

    @property
    def order(self):
        """Order and root of unit"""
        order, root_invers = None, 0
        for ior in range(1,7):
            rn = self ** ior

            if rn.isE:
                order = ior
                break

            if rn.isI: 
                root_invers = ior

        if order is None: 
            raise ValueError("symmetry is not a root of unit!")

        return order, root_invers

    #@property
    #def name(self):
    #    order, root_invers = self.info
    #    name = ""
    #    if self.det == -1: name = "-"
    #    name += str(self.order) # FIXME this one doesn't work yet.
    #    if root_invers != 0: 
    #        name += "-"
    #    else:
    #        name += "+"

    #    return name


#class Irrep(object):
#    """
#    This object represents an irreducible representation.
#
#    .. attributes:
#        nsym:
#            Number of symmetries.
#        dim:
#            Dimension of the irreducible representation.
#        all_traces:
#            all_traces[nsym]. The trace of each irrep.
#        character[num_classes]
#    """
#
#    def __init__(self, name, ops, mats):
#        """
#        Args:
#            name: 
#                Name of the irreducible representation.
#            ops:
#                List of symmetry operations packed in classes
#            mats:
#                Array of shape [nsym,dim,dim] with the irreducible 
#                representations of the group. mats are packed in classes.
#        """
#        self.name = name
#        self._ops = ops
#
#        assert len(mats.shape) == 3
#        self._mats = mats
#
#        self.dim = mats.shape[1]
#        assert self.dim == mats.shape[2]
#
#        self.all_traces = [m.trace() for m in mats]
#        # List of tuples, each tuple gives the start and stop index for the class.
#        [(0, 2), (2,4), (4,n)]
#        self.class_ranges = 

        # Compute character table.
        #character = []
        #for icls, (start, stop) in enumerate(self.class_ranges):
        #    t0 = self.all_traces[start]
        #    isok = all(t0 == self.all_traces[i] for i in range(start, stop)]
        #    character[icls] = t0

        #self._character = character

#    @property
#    def ops(self):
#        return self._ops

#    @property
#    def mats(self):
#        return self._mats
#
#    @property
#    def nsym(self)
#        return len(self.mats)
#
#    @property
#    def character(self)
#        return self._character


#ptg_filepath = os.path.join(..., "ptg_irreps.json")
#
#with open(ptg_filepath, "r") as fh:
#    d = json.load(fh)
#
#_PTG_IRREPS = {}
#for ptg_name, v in d.items():
#    _PTF_IRREPS[ptg_name] = Irrep(name, ops, mats)
#
#del d, ptg_filepath, ptg_name
#
#def ptgroup_irreps(ptg_name):
#    return _PTG_IRREPS[ptg_name]
#

class PointGroup(list):
    """
    A PointGroup is a list of Rotations and has irreducible representations
    """ 
    def __init__(self, rotations, name=None, irreprs=None):

        class_ids = mk_classes(rotations)

        # Always reorder rotations and irreprs according to class indeces.
        ord_rotations = [ None for ii in range(len(rotations)) ]
        idx = -1
        for ord_idx in rflat(class_ids): 
            idx += 1
            ord_rotations[idx] = rotations[ord_idx]

        ord_irreprs  = list() 
        for irr in irreprs:
            ord_matrices = [ None for ii in range(len(irr.matrices)) ]
            idx = -1
            for ord_idx in rflat(class_ids):
               idx += 1
               ord_matrices[idx] = irr.matrices[ord_idx]

            ord_irreprs.append( IrreducibleRepr(irr.name, irr.dim, ord_matrices) )

        list.__init__(self)
        for orot in ord_rotations: self.append(orot)

        self.class_ids = mk_classes(ord_rotations)
        self.nclass = len(self.class_ids)

        # Create name of each class.
        #self.class_names = [ "None" for ii in range(self.nclass) ]
        first_rot_ids = [ self.class_ids[ii][0] for ii in range(self.nclass) ]
        self.class_names = [ self[ii].name for ii in first_rot_ids ]

        self.nsym = len(self)
        self.name = str(name)

        self.irreprs = ord_irreprs
        #for ii in self.irreprs: print ii

        self.nirrepr = len(self.irreprs)

    #def find(self, rot):
    #    """Return the index of rot."""
    #    try:
    #        return self.index(rot)
    #    except ValueError:
    #        raise RotationNotFound(rot)

    #def findE(self):
    #    """Return the index of the identity."""
    #    try:
    #        return self.index(Rotation(_E3D))
    #    except RotationNotFound:
    #        raise

    #def show_character_table(self):
    #    vlen = 10
    #                                            
    #    print 100*"*"
    #    print ("Point Group" + self.name)
    #    cln = ""
    #    for clname in self.class_names:
    #        cln += str(clname).center(vlen)
    #    print "Class" + cln 
    #                                            
    #    mult = "Mult" 
    #    for cls in self.class_ids:
    #        mult += str(len(cls)).center(vlen)
    #    print mult
    #                                            
    #    for irrepr in self.irreprs:
    #        #print "irrepr ", irrepr
    #        row = irrepr.name.ljust(5)
    #        for icls in range(self.nclass): 
    #           sym_id = self.class_ids[icls][0] 
    #           mat = irrepr.matrices[sym_id]
    #           char = mat.trace()[0,0]
    #           row += str(char).center(vlen)
    #        print row
    #                                            
    #    print 100*"*"
    #    print 100*"*"

    def check(self):
        if not self.isgroup(): raise NotAGroup

        class_ids = mk_classes(self)
        #print class_ids
        check = -1
        for idx in rflat(class_ids):
            check = check +1
            if check!= idx: raise PointGroupException("Symmetries are not ordered by classes")

        mtable = self.mk_mtable()

        err = 0.0
        for idx1 in range(len(self)):
            for idx2 in range(len(self)):
                ij = (idx1, idx2)
                idx_prod = mtable[ij] 

                for irr in self.irreprs:
                    mat_prod = irr.matrices[idx1] * irr.matrices[idx2]
                    my_err = (mat_prod - irr.matrices[idx_prod]).max()
                    err = max(err, abs(my_err))

        print("Error in Group Representation", err)

        character_of = dict()
        for irr in self.irreprs:
            traces = irr.traces()
            #character = [ traces[ii] 
            chr = list()
            for clids in self.class_ids:
                idx = clids[0]
                chr.append(traces[idx])
            #character_of[irr.name] = N.array(chr)
            character_of[irr.name] = traces
            #irr.name 

        err_otrace = 0.0
        for k1, v1 in character_of.iteritems():
            for k2, v2 in character_of.iteritems():
                my_err = dotc(v1, v2) / self.nsym
                if k2 == k1: my_err -= 1.0 
                err_otrace = max(err_otrace, abs(my_err))
        print("Error in orthogonality relation of traces ", err)
