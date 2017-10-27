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
from abipy.electrons.ebands import ElectronsReader
#from abipy.dfpt.phonons import PhononBands, factor_ev2units, unit_tag, dos_label_from_units
#from abipy.abio.robots import Robot, RobotWithEbands


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
         app(marquee("SigmaPh calculation", mark="="))

         #if verbose > 1:
         #    app(marquee("Abinit Header", mark="="))
         #    app(self.hdr.to_string(verbose=verbose))

         return "\n".join(lines)

     #@lazy_property
     #def hdr(self):
     #    """:class:`AttrDict` with the Abinit header e.g. hdr.ecut."""
     #    return self.reader.read_abinit_hdr()

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

     @add_fig_kwargs
     def plot(self, what="lambda", units="eV", ylims=None, ax=None, **kwargs):
         """
         Plot phonon bands with eph coupling strenght lambda(q, nu) or gamma(q, nu)

         Args:
             what: `lambda` for eph strength, gamma for ph linewidth.
             units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                 Case-insensitive.
             scale: float used to scale the marker size.
             ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                    or scalar e.g. `left`. If left (right) is None, default values are used
             ax: matplotlib :class:`Axes` or None if a new figure should be created.

         Returns:
             `matplotlib` figure
         """
         ax, fig, plt = get_ax_fig_plt(ax=ax)

         # Plot phonon bands
         #self.phbands.plot(ax=ax, units=units, show=False)

         # Add eph coupling.
         #xvals = np.arange(len(self.phbands.qpoints))
         #yvals = self.phbands.phfreqs * factor_ev2units(units)

         # [0] is for the number_of_spins
         # TODO units
         #if what == "lambda":
         #    s = self.reader.read_phlambda_qpath()[0]
         #    scale = 100
         #elif what == "gamma":
         #    s = self.reader.read_phgamma_qpath()[0]
         #    scale = 1
         #else:
         #    raise ValueError("Invalid value for what: `%s`" % what)

         #color = "blue"
         #for nu in self.phbands.branches:
         #    ax.scatter(xvals, yvals[:, nu], s=scale * np.abs(s[:, nu]),
         #            c=color, marker="o", #label=term if ib == 0 else None
         #    )

         #xvals = np.tile(xvals, 3 * len(self.structure)).T
         #ax.scatter(xvals.T, yvals.T, s=scale * np.abs(s).T)

         #set_axlims(ax, ylims, "y")
         #return fig

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
            #nbv.new_code_cell("robot = abilab.SigEPhRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
        ])

        return self._write_nb_nbpath(nb, nbpath)
