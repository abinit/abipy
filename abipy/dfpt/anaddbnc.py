# coding: utf-8
"""
AnaddbNcFile provides a high-level interface to the data stored in the anaddb.nc file.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from monty.functools import lazy_property
from abipy.core.tensor import Tensor
from abipy.core.mixins import AbinitNcFile, Has_Structure
from abipy.iotools import ETSF_Reader
from abipy.dfpt.phonons import InteratomicForceConstants
from abipy.dfpt.ddb import Becs


class AnaddbNcFile(AbinitNcFile, Has_Structure):
    """
    AnaddbNcFile provides a high-level interface to the data stored in the anaddb.nc file.
    This object is usually instanciated with `abiopen("anaddb.nc")`.

    .. attribute:: structure

        Structure object.

    .. attribute:: emacro

        Macroscopic dielectric tensor. None if the file does not contain this information.

    .. attribute:: becs

        Born effective charges. None if the file does not contain this inf

    .. attribute:: ifc

        :class:`InteratomicForceConstants` object with the interatomic force constants calculated by anaddb.
        None, if the netcdf file does not contain the IFCs,
    """

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(AbinitNcFile, self).__init__(filepath)
        self.reader = ETSF_Reader(filepath)
        self._structure = self.reader.read_structure()

    def close(self):
        self.reader.close()

    @property
    def structure(self):
        return self._structure

    @lazy_property
    def emacro(self):
        """
        Macroscopic dielectric tensor. None if the file does not contain this information.
        """
        try:
            return Tensor.from_cartesian_tensor(self.reader.read_value("emacro_cart"),
                                                self.structure.lattice, space="r"),
        except Exception as exc:
            print(exc, "Returning None", sep="\n")
            return None

    @lazy_property
    def becs(self):
        """
        Born effective charges. None if the file does not contain this information.
        """
        try:
            chneut = -666 # TODO: anaddb.nc should contain the input file.
            return Becs(self.reader.read_value("becs_cart"), self.structure, chneut=chneut, order="f")
        except Exception as exc:
            print(exc, "Returning None", sep="\n")
            return None

    @lazy_property
    def ifc(self):
        """
        The interatomic force constants calculated by anaddb.
        The following anaddb variables should be used in the run: ifcflag, natifc, atifc, ifcout
        Return None, if the netcdf file does not contain the IFCs,
        """
        try:
            return InteratomicForceConstants.from_file(self.filepath)
        except Exception as exc:
            print(exc, "Returning None", sep="\n")
            return None
