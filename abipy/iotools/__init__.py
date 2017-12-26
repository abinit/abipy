# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

from .xsf import *
from .visualizer import *

import pymatgen.io.abinit.netcdf as ionc

as_etsfreader = ionc.as_etsfreader


class ETSF_Reader(ionc.ETSF_Reader):
    """
    Provides high-level API to read data from netcdf files written
    folloing the ETSF-IO specifications described in :cite:`Caliste2008`
    """

    def read_structure(self):
        """
        Overrides the ``read_structure`` method so that we always return
        an instance of AbiPy |Structure| object
        """
        from abipy.core.structure import Structure
        return Structure.from_file(self.path)
