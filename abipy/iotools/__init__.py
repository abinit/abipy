# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

from .xsf import *
from .visualizer import *

from monty.functools import lazy_property
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

    # Must overwrite implementation of pymatgen.io.abinit.netcdf
    # due to a possible bug introduced by initial whitespaces in symbol
    @lazy_property
    def chemical_symbols(self):
        """Chemical symbols char [number of atom species][symbol length]."""
        charr = self.read_value("chemical_symbols")
        symbols = []
        for v in charr:
            s = "".join(c.decode("utf-8") for c in v)
            # Strip to avoid possible whitespaces.
            symbols.append(s.strip())

        return symbols
