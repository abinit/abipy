# coding: utf-8
from __future__ import print_function, division, unicode_literals

from .xsf import *
from .visualizer import *

import pymatgen.io.abinit.netcdf as ionc 

as_etsfreader = ionc.as_etsfreader

class ETSF_Reader(ionc.ETSF_Reader):
    """
    Overrides the read_structure method so that we always return
    an instance of our Structure object
    """
    def read_structure(self):
        from abipy.core.structure import Structure
        return Structure.from_file(self.path)
