# coding: utf-8
"""Density/potential files in netcdf format."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from monty.string import marquee
from monty.termcolor import cprint
from monty.functools import lazy_property
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.fields import DensityReader
from abipy.electrons.ebands import ElectronsReader

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "DensityNcFile",
]

class DenNcReader(ElectronsReader, DensityReader):
    """Object used to read data from DEN.nc files."""


class DensityNcFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    Netcdf File containing the electronic density.

    Usage example:

    .. code-block:: python

        with DensityNcFile("foo_DEN.nc") as ncfile:
            ncfile.density
            ncfile.ebands.plot()
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def __init__(self, filepath):
        super(DensityNcFile, self).__init__(filepath)
        self.reader = DenNcReader(filepath)

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self):
        """Return string representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Structure", mark="="))
        app(str(self.structure))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))
        app("XC functional: %s" % str(self.xc))

        return "\n".join(lines)

    @lazy_property
    def density(self):
        return self.reader.read_density()

    @lazy_property
    def ebands(self):
        """:class:`ElectronBands` object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """:class:`Structure` object."""
        return self.ebands.structure

    @lazy_property
    def xc(self):
        """:class:`XcFunc` object with info on the exchange-correlation functional."""
        return self.reader.read_abinit_xcfunc()

    def close(self):
        self.reader.close()

    def write_chgcar(self, filename=None):
        """
        Write density in CHGCAR format. Return :class:`ChgCar` instance.
        """
        if filename is None:
            filename = self.basename.replace(".nc", "_CHGCAR")
            cprint("Writing density in CHGCAR format to file: %s" % filename, "yellow")
        return self.density.to_chgcar(filename=filename)

    def write_xsf(self, filename=None):
        """
        Write density in XSF format.
        """
        if filename is None:
            filename = self.basename.replace(".nc", ".xsf")
            cprint("Writing density in XSF format to file: %s" % filename, "yellow")
        return self.density.export(filename)

    def write_cube(self, filename=None, spin="total"):
        if filename is None:
            filename = self.basename.replace(".nc", ".cube")
            cprint("Writing density in CUBE format to file: %s" % filename, "yellow")
        return self.density.export_to_cube(filename, spin=spin)

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("denc = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(denc)"),
            nbv.new_code_cell("fig = denc.ebands.plot()"),
            nbv.new_code_cell("fig = denc.ebands.kpoints.plot()"),
            nbv.new_code_cell("fig = denc.ebands.get_edos().plot()"),
            nbv.new_code_cell("#cube = denc.write_cube(filename=None)"),
            nbv.new_code_cell("#xsf_path = denc.write_xsf(filename=None"),
            nbv.new_code_cell("#chgcar = denc.write_chgcar(filename=None"),
        ])

        return self._write_nb_nbpath(nb, nbpath)
