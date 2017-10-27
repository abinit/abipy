# coding: utf-8
"""
This module contains objects for the postprocessing of EPH calculations.

Warning:
    Work in progress, DO NOT USE THIS CODE.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pymatgen.core.units as units

from collections import OrderedDict
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import Kpath
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, set_axlims
from abipy.electrons.ebands import ElectronsReader
from abipy.dfpt.phonons import PhononBands, factor_ev2units, unit_tag, dos_label_from_units


class A2F(object):
    """
    Eliashberg function a2F(w).
    """

    def __init__(self, mesh, values_spin, values_spin_nu):
        """
        Eliashberg function
        vals(nomega,0:natom3,nsppol)
        vals(w,1:natom,1:nsppol): a2f(w) decomposed per phonon branch and spin
        vals(w,0,1:nsppol): a2f(w) summed over phonons modes, decomposed in spin
        """
        self.mesh = mesh
        # Spin dependent and total a2F(w)
        values_spin = np.atleast_2d(values_spin)
        values_spin_nu = np.atleast_3d(values_spin_nu)
        self.nsppol = len(values_spin)
        self.nmodes = self.values_spin_nu.shape[-1]

        if self.nsppol == 2:
            self.values = values_spin[0] + values_spin[1]
            self.values_nu = values_spin_nu[0] + values_spin_nu[1]
        elif self.nsppol == 1:
            self.values = 2 * values_spin[0]
            self.values_nu = 2 * values_spin_nu[0]
        else:
            raise ValueError("Invalid nsppol: %s" % self.nsppol)

        self.values_spin = values_spin
        self.values_spin_nu = values_spin_nu
        #self.lambdaw

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        lines = []
        return "\n".join(lines)

    #def get_momentum(self, order):
    #    """Computes the momenta of a2F(w)"""
    #    raise NotImplementedError()

    #def get_mcmillan_Tc(self, mustar):
    #    """
    #    Computes the critical temperature with the McMillan equation and the input mustar.
    #    """
    #    raise NotImplementedError()

    @add_fig_kwargs
    def plot(self, units="eV", exchange_xy=False, decompose=False, ax=None,
             xlims=None, ylims=None, label=None, **kwargs):
        """
        Plot a2F(w), it's primitive lambda(w) and optionally the individual contributions due
        to the phonon branches.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            exchange_xy: True to exchange x-y axes.
            decompose: True to plot the decomposition in terms of phonon branches.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            xlims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                    or scalar e.g. `left`. If left (right) is None, default values are used
            ylims: Limits for y-axis. See xlims for API.
            label: True to add legend label to each curve.

        Returns:
            `matplotlib` figure
        """""
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        def ax_plot(x, y, *args, **kwargs):
            if exchange_xy: x, y = y, x
            return ax.plot(x, y, *args, **kwargs)

        # TODO Handle styles
        # Plot a2f(w) and integral.
        ax.plot(self.mesh, self.values, label=lable)
        if self.nsppol == 2:
            for spin in range(self.nsppol):
                ax_plot(self.mesh, self.values_spin[spin])

        if decompose:
            for nu in range(self.nmodes):
                ax_plot(self.mesh, self.values_nu[nu])
            if self.nsppol == 2:
                for spin in range(self.nsppol):
                    for nu in range(self.nmodes):
                        ax_plot(self.mesh, self.values_spin_nu[spin, nu])

        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        return fig


#class A2Ftr(object):


# TODO Change name.
class EphFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the phonon linewidths, EliashbergFunction, the phonon bands,
    the ElectronBands on the k-mesh.
    Provides methods to plot results.

    Usage example:

    .. code-block:: python

        with EphFile("out_EPH.nc") as ncfile:
            print(ncfile)
            ncfile.ebands.plot()
            ncfile.phbands.plot()
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(EphFile, self).__init__(filepath)
        self.reader = EphReader(filepath)

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
        app(marquee("EPH calculation", mark="="))

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

    @property
    def phbands(self):
        """:class:`PhononBands` object with frequencies along the q-path."""
        return self.reader.read_phbands_qpath()

   #@lazy_property
   #def a2f(self):
   #    """:class:`A2f` with the Eliashberg function a2F(w)."""
   #    return self.reader.read_a2f()

   #@lazy_property
   #def a2ftr(self):
   #    """:class:`A2ftr` with the Eliashberg transport spectral function a2F_tr(w, x, x')."""
   #    return self.reader.read_a2ftr()

    def close(self):
        """Close the file."""
        self.reader.close()

    @add_fig_kwargs
    def plot_eph_strength(self, what="lambda", ylims=None, ax=None, **kwargs):
        """
        Plot phonon bands with eph coupling strenght lambda(q, nu)

        Args:
            what: `lambda` for eph strength, gamma for ph linewidth.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Plot phonon bands
        #self.phbands.plot(ax=ax, units=units, show=False)

        # Decorate the axis (e.g add ticks and labels).
        self.phbands.decorate_ax(ax, units="")

        # Add eph coupling.
        #natom3 = 3 * len(self.structure)
        #nqpt =  len(self.phbands.qpoints)
        #xvals = np.tile(np.arange(nqpt), natom3)
        #xvals = np.reshape(xvals, (natom3, nqpt))
        # This for the number_of_spins
        #yvals = self.reader.read_phlambda_qpath()[0]
        #s = self.reader.read_phgamma_qpath()[0]

        if what == "lambda":
            yvals = self.reader.read_phlambda_qpath()[0]
        elif what == "gamma":
            yvals = self.reader.read_phgamma_qpath()[0]
        else:
            raise ValueError("Invalid value for what: `%s`" % what)

        xvals = np.arange(len(self.phbands.qpoints))
        for nu in self.phbands.branches:
            ax.plot(xvals, yvals[:, nu])

        #if "color" not in kwargs and not match_bands:
        #    kwargs["color"] = "black"
        #if "linewidth" not in kwargs:
        #    kwargs["linewidth"] = 2.0
        #ax.plot(xvals, yvals.T)

        #set_axlims(ax, ylims, "y")
        return fig

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
        self.phbands.plot(ax=ax, units=units, show=False)

        # Add eph coupling.
        xvals = np.arange(len(self.phbands.qpoints))
        yvals = self.phbands.phfreqs * factor_ev2units(units)

        # [0] is for the number_of_spins
        # TODO units
        if what == "lambda":
            s = self.reader.read_phlambda_qpath()[0]
            scale = 100
        elif what == "gamma":
            s = self.reader.read_phgamma_qpath()[0]
            scale = 1
        else:
            raise ValueError("Invalid value fo what: `%s`" % what)

        color = "blue"
        for nu in self.phbands.branches:
            ax.scatter(xvals, yvals[:, nu], s=scale * np.abs(s[:, nu]),
                    c=color, marker="o", #label=term if ib == 0 else None
            )

        #xvals = np.tile(xvals, 3 * len(self.structure)).T
        #ax.scatter(xvals.T, yvals.T, s=scale * np.abs(s).T)

        set_axlims(ax, ylims, "y")
        return fig

    @add_fig_kwargs
    def plot_with_a2f(self, units="eV", ylims=None, **kwargs):
        """
        Plot phonon bands with lambda(q, nu) + a2F(w) + DOS.
        """
        # Build grid plot.
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        fig = plt.figure()
        ncols = 2
        #width_ratios = [2, 0.2, 0.2]
        width_ratios = [2, 0.2]
        gspec = GridSpec(1, ncols, width_ratios=width_ratios)
        gspec.update(wspace=0.05)

        ax_phbands = plt.subplot(gspec[0])
        ax_doses = []
        for i in range(ncols-1):
            ax = plt.subplot(gspec[i + 1], sharey=ax_phbands)
            ax.grid(True)
            set_axlims(ax, ylims, "y")
            ax_doses.append(ax)

        # Plot phonon bands with markers.
        self.plot(units=units, ylims=ylims, ax=ax_phbands, show=False)

        # Plot a2F(w)
        #self.a2f.plot(units=units, exchange_xy=True, ylims=ylims, ax=ax, show=False)

        # Plot a2Ftr(w)
        #self.a2ftr.plot(units=units, exchange_xy=True, ylims=ylims, ax=ax, show=False)

        # Plot DOS g(w)
        #phdos.plot(units=units, exchange_xy=True, ylims=ylims, ax=ax, show=False)

        return fig

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
            nbv.new_code_cell("ncfile.plot();"),
            nbv.new_code_cell("ncfile.plot_phlinewidths();"),
            nbv.new_code_cell("ncfile.plot_with_a2f();"),
            #nbv.new_code_cell("ncfile.a2f.plot();"),
            #nbv.new_code_cell("ncfile.plot_with_a2ftr();"),
            #nbv.new_code_cell("ncfile.a2ftr.plot();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class EphReader(ElectronsReader):
    """
    Reads data from file and constructs objects.
    """

    def read_phbands_qpath(self):
        """Read and return PhononBands."""
        structure = self.read_structure()

        # Build the list of q-points
        qpoints = Kpath(structure.reciprocal_lattice,
                        frac_coords=self.read_value("qpath"),
                        weights=None, names=None, ksampling=None)

        #nctkarr_t('phfreq_qpath', "dp", "natom3, nqpath, number_of_spins"),&
        phfreqs = self.read_value("phfreq_qpath")[0] * units.Ha_to_eV
        # TODO
        phdispl_cart = np.zeros((len(qpoints), 3*len(structure), 3*len(structure)))

        return PhononBands(structure=structure,
                           qpoints=qpoints,
                           phfreqs=phfreqs,
                           phdispl_cart=phdispl_cart,
                           non_anal_ph=None,
                           amu=self.read_value("atomic_mass_units"),
                           )

    def read_phlambda_qpath(self):
        """
        Reads the EPH coupling strength computed along the q-path.
        """
        #nctkarr_t('phlambda_qpath', "dp", "natom3, nqpath, number_of_spins")])
        return self.read_value("phlambda_qpath")

    def read_phgamma_qpath(self):
        """
        Reads the phonon linewidths computed along the q-path.
        """
        #nctkarr_t('phgamma_qpath', "dp", "natom3, nqpath, number_of_spins"),&
        return self.read_value("phgamma_qpath")

    #def read_a2f(self):
    #    """Read and return the Eliashberg function a2F(w)""".
    #    values_spin = self.read_value()
    #    values_spin_nu = self.read_value()
    #    return A2f(mesh=self.read_value(),
    #               values_spin,
    #               values_spin_nu
    #               )

    #def read_a2ftr(self):
    #    """Read and return the Eliashberg transport spectral function a2F_tr(w, x, x')."""
    #    if ...: return None
    #    return A2ftr()


from abipy.abio.robots import Robot, RobotWithEbands, RobotWithPhbands

class EphRobot(Robot, RobotWithEbands, RobotWithPhbands, NotebookWriter):
    """
    This robot analyzes the results contained in multiple EPH.nc files.
    """
    EXT = "EPH"

    @add_fig_kwargs
    def plot_a2f_convergence(self, sortby="nkpt", ax=None, xlims=None, **kwargs):
        """
        Plot the convergence of the Eliashberg function wrt to `sortby` parameter.

        Args:
            sortby: Define the convergence parameter, sort files and produce plot labels. Can be string or function.
                If string, it's assumed that the ncfile has an attribute with the same name and getattr is invoked.
                If callable, the output of callable(ncfile) is used.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used.

        Returns:
            `matplotlib` figure
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for label, ncfile, param in self.sortby(sortby):
            ncfile.a2f.plot(ax=ax, xlims=xlims,
                            label="%s %s" % (sortby, param) if not callable(sortby) else str(param),
                            show=False)
        return fig

    #@add_fig_kwargs
    #def plot_a2ftr_convergence(self, sortby="nkpt", ax=None, xlims=None, **kwargs):

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.EphRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("robot.plot_a2f_convergence();"),
        ])

        return self._write_nb_nbpath(nb, nbpath)
