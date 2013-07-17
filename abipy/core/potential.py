"""This module contains the class describing local potentials in real space on uniform 3D meshes."""
from __future__ import division, print_function

import numpy as np

from .structure import Structure
from .dftscalarfield import DFTScalarField

__all__ = [
    "Potential",
]

#########################################################################################
# Global variables

#: Flags denoting the potential type. Note bit order.
VNONE = 1
VLOC  = 2
VH    = 4
VX    = 8
VC    = 16

#: Mapping flag --> potential name.
_vnames =  {
VNONE: "None potential",
VLOC : "Local potential",
VH   : "Hartree potential",
VX   : "Exchange potential",
VC   : "Correlation potential",
}

#: List of allowed potential types.
VTYPES = _vnames.keys()

#########################################################################################


class PotentialInfo(object):

    def __init__(self, vtype, **kwargs):
        self.vtype = vtype
        self._dict = kwargs

    def __str__(self):
        s = self.vname() + "\n"
        for k, v in self._dict.items:
            s += str(k) + " = " + str(v) + "\n"
        return s

    def __add__(self, other):
        new_vtype = self.vtype + other.vtype

        new_dict = self._dict.copy()

        for k in other._dict:
            if k not in self._dict:
                new_dict[k] = other._dict[k]
            else:
                # Append values if key is already present.
                lst = list()
                lst.append(new_dict[k])
                lst.append(other._dict[k])
                new_dict[k] = lst

        return PotentialInfo(vtype, kwargs=new_dict)

    @property
    def vname(self):
        """Return the name of the potential."""
        s = bin(self.vtype)[2:][::-1]
        flags = [ 2**idx for idx, c in enumerate(s) if c == "1" ]

        plus = ""
        if len(flags) > 1: plus = "+"
        return plus.join( [_vnames[flag] for flag in flags] )

#########################################################################################


class Potential(DFTScalarField):
    #
    #: Attributes read from the netCDF file.
    _slots = [
      "cplex_pot",                # "real_or_complex_potential"
      "nsppol",                   # "number_of_spins"
      "nspinor",                  # "number_of_spinor_components"
      "nspden",                   # "number_of_components"
      "nfft1",                    # "number_of_grid_points_vector1"
      "nfft2",                    # "number_of_grid_points_vector2"
      "nfft3",                    # "number_of_grid_points_vector3"
      "nelect",                   # "number_of_electrons"
      "valence_charges",          # "valence_charges"
    ]
                                                                    
    _potopts = [
      "exchange_potential",              #"vx",
      "correlation_potential",           #"vc",
      "exchange_correlation_potential",  #"vxc",
      "exchange_functional",             #"exchange_functional"
      "correlation_functional",          #"correlation_functional"
    ]

    @classmethod
    def from_file(cls, path):
        """
        Read density from a the netCDF file.

        Args:
            path:
                string or file object.
        returns:
            :class:`Potential`
        """
        raise NotImplementedError("Potential must be rewritten from scrath")
        # Read all the keywords present in fname so that we know the potential type.
        #dim_names, var_names = ncread_keys(path)

        nfound = 0
        for pot in cls._potopts:
            if pot in var_names:
                nfound += 1
                v = Record()
                pot_vars = ["data"]
                map["data"] = pot
                pot_vars += _potopts[pot]

                #missing = ncread_varsdims(v, path, pot_vars, map_names=map)
                #
                # Handle the error
                if missing:
                    for miss in missing:
                        print("potential variables %s are missing!" % str(miss[0]))
                    raise ValueError("Fatal Eror")

        if nfound != 1:  # No potential or more than one.
            raise ValueError("Number of potential found in file is %s" % nfound )

        # Build potential descriptor.
        #vinfo = PotentialInfo(vtype, **kwargs)

        structure = Structure.from_file(path)

        rec = Record()
        #
        # Handle the error
        if missing:
            for miss in missing:
                print("internal name= " + miss[0] + ", ETSF name= " + miss[1] + " is missing!")

        # use iorder="f" to transpose the last 3 dimensions since ETSF
        # stores data in Fortran order while abipy uses C-ordering.

        if rec.cplex_pot == 1:
            # Get rid of fake last dimensions (cplex).
            vr = np.reshape(v.data, (rec.nspden, rec.nfft3, rec.nfft2, rec.nfft1))
            return Potential(rec.nspinor, rec.nsppol, rec.nspden, vtype, vr, structure, iorder="f")

        else:
            raise NotImplementedError("cplex_den = " + str(rec.cplex_den) + "not coded")


    def __init__(self, nspinor, nsppol, nspden, vtype, vr, structure, iorder="c"):
        """
        Args:
            nspinor:
                Number of spinorial components.
            nsppol:
                Number of spins.
            nspden:
                Number of spin density components.
            vtype:
                Flag defining the potential type.
            vr:
                numpy array with the potential in real space.
            structure:
                pymatgen structure
            iorder:
                Order of the array. "c" for C ordering, "f" for Fortran ordering.
        """
        super(Potential, self).__init__(nspinor, nsppol, nspden, vr, structure, iorder=iorder)

        if vtype not in VTYPES:
            raise ValueError("Unknow vtype: " + str(vtype))
        self.vtype = vtype

    @property
    def vname(self):
        """The name of the potential."""
        s = bin(self.vtype)[2:][::-1]
        flags = [2**idx for idx, c in enumerate(s) if c == "1"]

        plus = ""
        if len(flags) > 1: plus = "+"
        return plus.join([_vnames[flag] for flag in flags])

    def tostring(self, prtvol=0):
        s = self.vname + "\n"
        s += DFTScalarField.tostring(self, prtvol)
        return s

    #def make_vector_field(self):
    #  """Return vector field."""
    #  if self.iscollinear: return None

#########################################################################################
