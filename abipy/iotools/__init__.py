from .etsfio import *
from .xsf import *
from .visualizer import *
from .files import *

import pymatgen.io.abinitio.netcdf as ionc 

as_etsfreader = ionc.as_etsfreader

class ETSF_Reader(ionc.ETSF_Reader):
    """
    Overrides the read_structure method so that we always return
    an instance of our Structure object
    """
    def read_structure(self):
        from abipy.core.structure import Structure
        return Structure.from_file(self.path)

    def read_ebands(self):
        from abipy.electrons.electrons import ElectronBands
        return ElectronBands.from_file(self.path)

    def read_kpoints(self):
        from abipy.core.kpoints import kpoints_factory
        return kpoints_factory(self.path)
