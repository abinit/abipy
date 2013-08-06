"""GSR file."""
from __future__ import print_function, division

import numpy as np

from abipy.iotools import ETSF_Reader, AbinitNcFile, Has_Structure, Has_ElectronBands

__all__ = [
    "GSR_File",
]


class GSR_File(AbinitNcFile, Has_Structure, Has_ElectronBands):
    """
    File containing the results of a Ground-state calculation.
    """
    def __init__(self, filepath):
        super(GSR_File, self).__init__(filepath)

        with GSR_Reader(filepath) as r:
            # Initialize the  structure from file.
            self._structure = r.read_structure()

            # Initialize the band energies.
            self._bands = ElectronBands.from_file(filepath)

            #self.kpoints = kpoints_factory(filepath)
            #self.nkpt = len(self.kpoints)

            #self.nspinor = r.nspinor
            #self.nsppol = r.nsppol
            #self.nspden = r.nspden

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    @property
    def structure(self):
        """`Structure` object"""
        return self._structure
                                     
    @property
    def bands(self):
        """`ElectronBands` object"""
        return self._bands


class GSR_Reader(ETSF_Reader):
    """
    This object reads the results stored in the _GSR (Ground-State Results)
    file produced by ABINIT. It provides helper function to access the most
    important quantities.
    """
    def __init__(self, filepath):
        """Initialize the object from a filename."""
        super(GSR_Reader, self).__init__(filepath)

        #structure = self.read_structure()
        #self.kibz = kpoints_factory(path)

        #fermie = Ha2eV(self.read_value("fermie"))
        #np_eigvals = Ha2eV(self.read_value("eigenvalues"))
        # TODO
        #assert np_eigvals.units == "atomic units"
        #nsppol = np_eigvals.shape[0]

    #def read_forces(self):

    #def read_stress(self):

    #def read_energies(self):
    #    return AttrDict()
