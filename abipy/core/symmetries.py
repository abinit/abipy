"""Objects used to deal with symmetry operations in crystals."""
from __future__ import division, print_function

import numpy as np

from abipy.kpoints.utils import wrap_to_ws
from abipy.iotools import as_etsfreader

__all__ = [
    "SymmOp",
    "SpaceGroup",
]


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

#########################################################################################


class SymmOp(object):

    def __init__(self, rot_r, tau, time_sgn, afm_sgn):
        """
        One crystalline symmetry.

        Args:
            rot_r:
                3x3 integer matrix with the rotational part in reduced coordinates (C order).
            tau:
                fractional translation in reduced coordinates.
            time_sgn:
                -1 if time reversal can be used, otherwise +1.
            afm_sgn:
                anti-ferromagnetic part [+1,-1].
        """
        rot_r = np.asarray(rot_r)
        self.rot_r = rot_r
        self.rotm1_r = mati3inv(rot_r, trans=False) # R^{-1}
        self.tau = np.asarray(tau)
        self.afm_sgn  = afm_sgn;  assert afm_sgn in  [-1, 1]
        self.time_sgn = time_sgn; assert time_sgn in [-1, 1]

        # Compute symmetry in reciprocal space: S = R^{-1t}
        self.rot_g = mati3inv(rot_r, trans=True)

    def __str__(self):
        return self.tostring()

    def tostring(self, prtvol=0):
        s = str(self.rot_r) + "\n"
        s += "tau = %s, time_sgn = %s, afm_sgn = %s\n" % (self.tau, self.time_sgn, self.afm_sgn)
        return s

    def __eq__(self, other):
        return (np.allclose(self.rot_r, other.rot_r) and
                np.allclose(self.tau, other.tau) and  # FIXME Should we allow for a Bravais lattice?
                self.afm_sgn  == other.afm_sgn and
                self.time_sgn == other.time_sgn
                )

    def __ne__(self, other):
        return not self == other

    @property
    def det(self):
        """Determinant of the rotation matrix [-1, +1]."""
        return _get_det(self.rot_r)

    @property
    def trace(self):
        """Trace of the rotation matrix."""
        return self.rot_r.trace()

    @property
    def isproper(self):
        """True if the rotational part has determinant == 1."""
        return self.det == +1

    @property
    def has_timerev(self):
        """True if symmetry contains the time-reversal operator."""
        return self.time_sgn == -1

    @property
    def isafm(self):
        """True if anti-ferromagnetic symmetry."""
        return self.afm_sgn == -1

    def rotate_k(self, kpoint, wrap_tows=True):
        """
        Apply the symmetry operation to the k-point kpoint given in reduced coordinates.

        Sk is wrapped to the first Brillouin zone if wrap is True.
        """
        sk = np.dot(self.rot_g, kpoint) * self.time_sgn
        if not wrap_tows:
            return sk
        else:
            return wrap_to_ws(sk)

    def rotate_g(self, gvecs):
        """
        Apply the symmetry operation to the list of gvectors gvecs in reduced coordinates.
        """
        rot_gvecs = np.zeros_like(gvecs)
        if self.time_sgn == 1:
            for ig, gvec in enumerate(gvecs):
                rot_gvecs[ig] = np.dot(self.rot_g, gvec)
        else:
            for ig, gvec in enumerate(gvecs):
                rot_gvecs[ig] = np.dot(self.rot_g, gvec) * self.time_sgn
        return rot_gvecs

#########################################################################################


class SpaceGroup(object):
    """Container storing the space group symmetries."""

    def __init__(self, spgid, symrel, tnons, symafm, has_timerev, inord="c"):
        """
        Creation method.

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
                storage order of mat in symrel[:]. If inord == "f", mat.T is stored
                as matrices are always stored in C-order. Use inord == "f" if you have read
                symrel from an external file produced by abinit.
        """
        self.spgid = spgid
        assert self.spgid in range(0, 233)

        time_signs = [+1, -1]
        if not has_timerev:
            time_signs = [+1]

        self._time_signs = time_signs
        self.has_timerev = has_timerev

        nsym = symrel.shape[0]
        if tnons.shape[0] != nsym or symafm.shape[0] != nsym:
            raise ValueError(" symrel, tnons and symafm must have equal shape[0]")

        sym_list = []
        for time_sign in time_signs:
            for idx in range(nsym):
                rot_r = np.asarray(symrel[idx])
                if inord.lower() == "f":
                    rot_r = rot_r.T # Fortran to C.
                tau = tnons[idx]
                afm_sgn = symafm[idx]

                sym = SymmOp(rot_r, tau, time_sign, afm_sgn)
                sym_list.append(sym)

        self._sym_tuple = tuple(sym_list)
        self.nsym = len(self._sym_tuple)

    @classmethod
    def from_file(cls, file):
        """Initialize the object from a Netcdf file."""
        file, closeit = as_etsfreader(file)

        new = cls(spgid=file.read_value("space_group"),
                  symrel=file.read_value("reduced_symmetry_matrices"),
                  tnons=file.read_value("reduced_symmetry_translations"),
                  symafm=file.read_value("symafm"),
                  has_timerev=True,  # FIXME not treated by ETSF-IO.
                  inord="f",
                  )
        if closeit:
            file.close()
        return new

    def __len__(self):
        return self.nsym

    def __getitem__(self, slice):
        return self._sym_tuple[slice]

    def __str__(self):
        return self.tostring()

    @property
    def afmsyms(self):
        """Tuple with antiferromagnetic symmetries."""
        return self.symmops(time_sgn=None, afm_sgn=-1)
                                                       
    @property
    def fmsyms(self):
        """Tuple of ferromagnetic symmetries."""
        return self.symmops(time_sgn=None, afm_sgn=+1)

    def tostring(self, prtvol=0):
        """String representation."""
        lines = []
        app = lines.append
        app("nsym = %d" % self.nsym)
        app(" has_timerev = %s" % self.has_timerev)
        for sym in self._sym_tuple:
            app(str(sym))
        return "\n".join(lines)

    def symmops(self, time_sgn=None, afm_sgn=None):
        """
        Args:
            time_sgn:
                If specified, only symmetryes with time-reversal sign time_sgn are returned.
            afm_sgn:
                If specified, only symmetryes with anti-ferromagnetic part afm_sgn are returned.
        returns:
            tuple of :class:`SymmOp` instances.
        """
        syms = []
        for s in self._sym_tuple:
            gotit = True
            if time_sgn:
                assert time_sgn in [-1,+1]
                gotit = gotit and s.time_sgn == time_sgn
            if afm_sgn:
                assert afm_sgn in [-1,+1]
                gotit = gotit and s.afm_sgn == afm_sgn
            if gotit: syms.append(s)

        return tuple(syms)

