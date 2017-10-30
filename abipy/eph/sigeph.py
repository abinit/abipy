# coding: utf-8
"""
This module contains objects for the postprocessing of Sigma_eph calculations.

Warning:

    Work in progress, DO NOT USE THIS CODE
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

#from collections import OrderedDict
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
#from abipy.core.kpoints import Kpath
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims
from abipy.electrons.ebands import ElectronsReader, RobotWithEbands
#from abipy.dfpt.phonons import PhononBands, RobotWithPhbands, factor_ev2units, unit_tag, dos_label_from_units
from abipy.abio.robots import Robot


# TODO QPState and QPList from electrons.gw (Define base abstract class).


class SigEPhFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the Fan-Migdal self-energy, the ElectronBands on the k-mesh.
    Provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with SigEPhFile("out_SIGEPH.nc") as ncfile:
            print(ncfile)
            ncfile.ebands.plot()
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(SigEPhFile, self).__init__(filepath)
        self.reader = SigmaPhReader(filepath)

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(marquee("Structure", mark="="))
        app(str(self.structure))
        app("")
        app(self.ebands.to_string(with_structure=False, title="Electronic Bands"))
        app("")
        # SigmaPh section
        app(marquee("SigmaPh calculation", mark="="))

        return "\n".join(lines)

    @lazy_property
    def ebands(self):
        """:class:`ElectronBands` object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """:class:`Structure` object."""
        return self.ebands.structure

    #@property
    #def phbands(self):
    #    """:class:`PhononBands` object with frequencies along the q-path."""
    #    return self.reader.read_phbands_qpath()

    def close(self):
         """Close the file."""
         self.reader.close()

    #@lazy_property
    #def params(self):
    #    """AttrDict dictionary with the GW convergence parameters, e.g. ecuteps"""
    #    return self.reader.read_params()

    #@lazy_property
    #def qplist_spin(self):
    #    """Tuple of :class:`QPList` objects indexed by spin."""
    #    return self.reader.read_allqps()

    #def get_sigmaw(self, spin, kpoint, band):
    #    """"
    #    Read self-energy(w) for (spin, kpoint, band)
    #    Return :class:`Function1D` object
    #    """
    #    wmesh, sigxc_values = self.reader.read_sigmaw(spin, kpoint, band)
    #    wmesh, spf_values = self.reader.read_spfunc(spin, kpoint, band)

    #    return Sigmaw(spin, kpoint, band, wmesh, sigxc_values, spf_values)

    #def get_dataframe(self):
    #def get_dataframe_kpoint(self):

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.ebands.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class SigEphRobot(Robot, RobotWithEbands, NotebookWriter):
    """
    This robot analyzes the results contained in multiple SIGEPH.nc files.
    """
    EXT = "SIGEPH"

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.SigEPhRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class SigmaPhReader(ElectronsReader):
    """
    Reads data from file and constructs objects.
    """
    #def read_params(self):
    #    """Read sigeph input parameters. Return OrderedDict"""
    #    od = OrderedDict()
    #    return od

    #def read_phbands_qpath(self):
    #    """Read and return PhononBands."""
    #    structure = self.read_structure()

    #    # Build the list of q-points
    #    qpoints = Kpath(structure.reciprocal_lattice,
    #                    frac_coords=self.read_value("qpath"),
    #                    weights=None, names=None, ksampling=None)

    #    #nctkarr_t('phfreq_qpath', "dp", "natom3, nqpath, number_of_spins"),&
    #    phfreqs = self.read_value("phfreq_qpath")[0] * units.Ha_to_eV
    #    # TODO
    #    phdispl_cart = np.zeros((len(qpoints), 3*len(structure), 3*len(structure)))

    #    return PhononBands(structure=structure,
    #                       qpoints=qpoints,
    #                       phfreqs=phfreqs,
    #                       phdispl_cart=phdispl_cart,
    #                       non_anal_ph=None,
    #                       amu=self.read_value("atomic_mass_units"),
    #                       )
