# coding: utf-8
"""
Interface to the ESKW.nc file storing the (star-function) interpolated band structure produced by Abinit.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from monty.functools import lazy_property
from monty.string import marquee
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.electrons.ebands import ElectronsReader


class EskwFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the (star-function) interpolated band structure produced by Abinit.
    It's similar to the GSR file but it does not contain the header and energies.

    Usage example:

    .. code-block:: python

        with EskwFile("foo_ESKW.nc") as eskw:
            eskw.ebands.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EskwFile
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(EskwFile, self).__init__(filepath)
        self.reader = ElectronsReader(filepath)

    @lazy_property
    def einterp(self):
        return self.reader.read_value("einterp")

    @lazy_property
    def band_block(self):
        # band_block(2)=Initial and final band index to be interpolated. [0, 0] if all bands are used.
        band_block = self.reader.read_value("band_block")
        if all(band_block != [0, 0]): band_block -= 1
        return band_block

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        app("band_block: %s" % str(self.band_block))
        app("einterp: %s" % str(self.einterp))

        return "\n".join(lines)

    def close(self):
        self.reader.close()

    @property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        od["einterp"] = self.interp
        od["einterp"] = self.einterp
        return od

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        #for fig in self.yield_structure_figs(**kwargs): yield fig
        for fig in self.yield_ebands_figs(**kwargs): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("eskw = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(eskw)"),
            nbv.new_code_cell("eskw.ebands.plot();"),
            nbv.new_code_cell("eskw.ebands.kpoints.plot();"),
            nbv.new_code_cell("# eskw.ebands.plot_transitions(omega_ev=3.0, qpt=(0, 0, 0), atol_ev=0.1);"),
            nbv.new_code_cell("""\
if eskw.ebands.kpoints.is_ibz:
    eskw.ebands.get_edos().plot();"""),
        ])

        return self._write_nb_nbpath(nb, nbpath)
