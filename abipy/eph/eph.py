# coding: utf-8
"""
This module contains objects for the postprocessing of EPH calculations.

Warning:
    Work in progress, DO NOT USE THIS CODE.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pymatgen.core.units as units
import abipy.core.abinit_units as abu

from collections import OrderedDict
from scipy.integrate import cumtrapz, simps
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import Kpath
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims
from abipy.electrons.ebands import ElectronsReader, ElectronDos, RobotWithEbands
from abipy.dfpt.phonons import PhononBands, PhononDos, RobotWithPhbands
from abipy.abio.robots import Robot


class A2f(object):
    """
    Eliashberg function a2F(w). Energies are in eV.
    """
    # Markers used for up/down bands.
    marker_spin = {0: "^", 1: "v"}

    def __init__(self, mesh, values_spin, values_spin_nu):
        """

        Args:
            mesh: Energy mesh in eV

        values(nomega,0:natom3,nsppol)
        vals(w,1:natom,1:nsppol): a2f(w) decomposed per phonon branch and spin
        vals(w,0,1:nsppol): a2f(w) summed over phonons modes, decomposed in spin

        TODO: Add metadata (qsampling, broadening, possibility of computing a2f directly?
        """
        self.mesh = mesh

        # Spin dependent and total a2F(w)
        values_spin = np.atleast_2d(values_spin)
        values_spin_nu = np.atleast_3d(values_spin_nu)
        self.nsppol = len(values_spin)
        self.nmodes = values_spin_nu.shape[1]
        assert self.nmodes % 3 == 0
        self.natom = self.nmodes // 3

        if self.nsppol == 2:
            self.values = values_spin[0] + values_spin[1]
            self.values_nu = values_spin_nu[0] + values_spin_nu[1]
        elif self.nsppol == 1:
            self.values = values_spin[0]
            self.values_nu = values_spin_nu[0]
        else:
            raise ValueError("Invalid nsppol: %s" % self.nsppol)

        self.values_spin = values_spin
        self.values_spin_nu = values_spin_nu
        #self.lambdaw ?

    @lazy_property
    def iw0(self):
        """
        Index of the first point in the mesh whose value is >= 0
        Integrals are performed with wmesh[iw0 + 1, :] i.e. unstable modes are neglected.
        """
        for i, x in enumerate(self.mesh):
            if x >= 0.0: return i
        else:
            raise ValueError("Cannot find zero in energy mesh")

    def __str__(self):
        return self.to_string()

    def to_string(self, title=None, verbose=0):
        """String representation with verbosity level ``verbose`` and an optional ``title``."""
        lines = []; app = lines.append

        app("Eliashberg Function" if not title else str(title))
        app("Mesh from %.4f to %.4f [eV] with %d points" % (
            self.mesh[0], self.mesh[-1], len(self.mesh)))

        # TODO: Add ElectronDos
        app("Isotropic lambda: %.3f" % (self.lambda_iso))
        app("Omega_log: %s [eV], %s [K]" % (self.omega_log, self.omega_log * abu.eV_to_K))
        for mustar in (0.1, 0.12, 0.2):
            app("\tFor mustar %s: McMillan Tc: %s [K]" % (mustar, self.get_mcmillan_tc(mustar)))

        if verbose:
            # $\int dw [a2F(w)/w] w^n$
            for n in [0, 4]:
                app("Moment %s: %s" % (n, self.get_moment(n)))

        return "\n".join(lines)

    @lazy_property
    def lambda_iso(self):
        """Isotropic lambda."""
        return self.get_moment(n=0)

    @lazy_property
    def omega_log(self):
        r"""
        Logarithmic moment of alpha^2F: exp((2/\lambda) \int dw a2F(w) ln(w)/w)
        """
        #return 270 / abu.eV_to_K
        iw = self.iw0 + 1
        wmesh, a2fw = self.mesh[iw:], self.values[iw:]

        #ax, fig, plt = get_ax_fig_plt(ax=None)
        #ax.plot(wmesh, a2fw / wmesh * np.log(wmesh))
        #plt.show()

        fw = a2fw / wmesh * np.log(wmesh)
        integral = simps(fw, x=wmesh)
        return np.exp(1.0 / self.lambda_iso * integral)

    def get_moment(self, n, spin=None, cumulative=False):
        r"""
        Computes the moment of a2F(w) i.e. $\int dw [a2F(w)/w] w^n$
        From Allen PRL 59 1460 (See also Grimvall, Eq 6.72 page 175)
        """
        wmesh = self.mesh[self.iw0+1:]
        if spin is None:
            a2fw = self.values[self.iw0+1:]
        else:
            a2fw = self.values_spin[spin][self.iw0+1:]

        # Primitive is given on the same mesh as self.
        ff = a2fw * (wmesh ** (n - 1))
        vals = np.zeros(self.mesh.shape)
        vals[self.iw0+1:] = cumtrapz(ff, x=wmesh, initial=0.0)
        return vals if cumulative else vals[-1].copy()

    def get_moment_nu(self, n, nu, spin=None, cumulative=False):
        r"""
        Computes the moment of a2F(w) i.e. $\int dw [a2F(w)/w] w^n$
        From Allen PRL 59 1460 (See also Grimvall, Eq 6.72 page 175)
        """
        wmesh = self.mesh[self.iw0+1:]
        if spin is None:
            a2fw = self.values_nu[nu][self.iw0+1:]
        else:
            a2fw = self.values_spin_nu[spin][nu][self.iw0+1:]

        # Primitive is given on the same mesh as self.
        ff = a2fw * (wmesh ** (n - 1))
        vals = np.zeros(self.mesh.shape)
        vals[self.iw0+1:] = cumtrapz(ff, x=wmesh, initial=0.0)
        return vals if cumulative else vals[-1].copy()

    def get_mcmillan_tc(self, mustar):
        """
        Computes the critical temperature with the McMillan equation and the input mustar.

        Return: Tc in Kelvin.
        """
        tc = (self.omega_log / 1.2) * \
            np.exp(-1.04 * (1.0 + self.lambda_iso) / (self.lambda_iso - mustar * (1.0 + 0.62 * self.lambda_iso)))
        return tc * abu.eV_to_K

    def get_mustar_from_tc(self, tc):
        """
        Return the value of mustar that gives the critical temperature ``tc`` in the McMillan equation.

        Args:
            tc: Critical temperature in Kelvin.
        """
        l = self.lambda_iso
        num = l + (1.04 * (1 + l) / np.log(1.2 * abu.kb_eVK * tc / self.omega_log))
        return num / (1 + 0.62 * l)

    @add_fig_kwargs
    def plot(self, units="eV", with_lambda=True, exchange_xy=False, ax=None,
             xlims=None, ylims=None, label=None, fontsize=12, **kwargs):
        """
        Plot a2F(w), its primitive lambda(w) and optionally the individual
        contributions due to the phonon branches.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            with_lambda: True to display lambda(q, nu).
            exchange_xy: True to exchange x-y axes.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
		or scalar e.g. ``left``. If left (right) is None, default values are used
            ylims: Limits for y-axis. See xlims for API.
            label: True to add legend label to each curve.
            fontsize: Legend and title fontsize

        Returns: |matplotlib-Figure|
        """""
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        wfactor = abu.phfactor_ev2units(units)

        # TODO Better handling of styles
        style = dict(
            linestyle=kwargs.pop("linestyle", "-"),
            color=kwargs.pop("color", "k"),
            linewidth=kwargs.pop("linewidth", 1),
        )

        # Plot a2f(w)
        xx, yy = self.mesh * wfactor, self.values
        if exchange_xy: xx, yy = yy, xx
        ax.plot(xx, yy, label=label, **style)

        # Plot lambda(w)
        if with_lambda:
            lambda_w = self.get_moment(n=0, cumulative=True)
            l_ax = ax.twinx()
            xx, yy = self.mesh * wfactor, lambda_w
            if exchange_xy: xx, yy = yy, xx
            l_ax.plot(xx, yy, label=label, **style)
            l_ax.set_ylabel(r"$\lambda(\omega)$")

        if self.nsppol == 2:
            # Plot spin resolved a2f(w).
            for spin in range(self.nsppol):
                xx, yy = self.mesh * wfactor, self.values_spin[spin]
                if exchange_xy: xx, yy = yy, xx
                ax_plot(xx, yy, marker=self.marker_spin[spin], **style)

        xlabel = abu.wlabel_from_units(units)
        ylabel = r"$\alpha^2F(\omega)$"
        if exchange_xy: xlabel, ylabel = ylabel, xlabel
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        ax.grid(True)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        if label: ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_nuterms(self, units="eV", ax_mat=None, with_lambda=True, fontsize=12,
                     xlims=None, ylims=None, label=None, **kwargs):
        """
        Plot a2F(w), its primitive lambda(w) and optionally the individual
        contributions due to the phonon branches.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            ax_mat: Matrix of axis of shape [natom, 3]. None if a new figure should be created.
            fontsize: Legend and title fontsize.
            xlims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
		or scalar e.g. ``left``. If left (right) is None, default values are used
            ylims: Limits for y-axis. See xlims for API.
            label: True to add legend label to each curve.

        Returns: |matplotlib-Figure|
        """""
        # Get ax_mat and fig.
        nrows, ncols = self.natom, 3
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_mat = np.reshape(ax_mat, (self.natom, 3))

        wfactor = abu.phfactor_ev2units(units)
        wvals = self.mesh * wfactor

        if with_lambda:
            lax_nu = [ax.twinx() for ax in ax_mat.flat]
            # Share axis after creation. Based on
            # https://stackoverflow.com/questions/42973223/how-share-x-axis-of-two-subplots-after-they-are-created
            lax_nu[0].get_shared_x_axes().join(*lax_nu)
            lax_nu[0].get_shared_y_axes().join(*lax_nu)
            for i, ax in enumerate(lax_nu):
                if i == 2: continue
                ax.set_yticklabels([])
                #ax.set_xticklabels([])

        # TODO Better handling of styles
        a2f_style = dict(
            linestyle=kwargs.pop("linestyle", "-"),
            color=kwargs.pop("color", "k"),
            linewidth=kwargs.pop("linewidth", 1),
        )
        lambda_style = a2f_style.copy()
        lambda_style["color"] = "red"

        import itertools
        for idir, iatom in itertools.product(range(3), range(self.natom)):
            nu = idir + 3 * iatom
            ax = ax_mat[iatom, idir]
            ax.grid(True)
            ax.set_title(r"$\nu = %d$" % nu, fontsize=fontsize)
            if idir == 0:
                ax.set_ylabel(r"$\alpha^2F(\omega)$")
            else:
                pass
                # Turn off tick labels
                #ax.set_yticklabels([])
                #ax.set_yticks([])

            if iatom == self.natom -1:
                ax.set_xlabel(abu.wlabel_from_units(units))
            #set_axlims(ax, xlims, "x")
            #set_axlims(ax, ylims, "y")

            # Plot total a2f(w)
            ax.plot(wvals, self.values_nu[nu], **a2f_style)

            # Plot lambda(w)
            if with_lambda:
                lambdaw_nu = self.get_moment_nu(n=0, nu=nu, cumulative=True)
                lax = lax_nu[nu]
                lax.plot(wvals, lambdaw_nu, **lambda_style)
                if idir == 2:
                    lax.set_ylabel(r"$\lambda_{\nu}(\omega)$", color=lambda_style["color"])

            #if self.nsppol == 2:
            #   # Plot spin resolved a2f(w)
            #   for spin in range(self.nsppol):
            #       ax.plot(self.mesh, self.values_spin_nu[spin, nu],
            #               marker=self.marker_spin[spin], **a2f_style)

        return fig

    @add_fig_kwargs
    def plot_a2(self, phdos, atol=1e-12, **kwargs):
        """
        Grid with 3 plots (a2F, F, a2F).

        Args:
            phdos:
            atol:

        Returns: |matplotlib-Figure|
        """
        phdos = PhononDos.as_phdos(phdos)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=3, ncols=1,
                                                sharex=True, sharey=False, squeeze=True)
        ax_list = ax_list.ravel()

        # Spline phdos onto a2f mesh and compute a2F(w) / F(w)
        f = phdos.spline(self.mesh)
        f = self.values / np.where(np.abs(f) > atol, f, atol)
        ax = ax_list[0]
        ax.plot(self.mesh, f, color="k", linestyle="-")
        ax.grid(True)
        ax.set_ylabel(r"$\alpha^2(\omega)$ [1/eV]")

        # Plot F. TODO: This should not be called plot_dos_idos!
        ax = ax_list[1]
        phdos.plot_dos_idos(ax=ax, what="d", color="k", linestyle="-")
        ax.grid(True)
        ax.set_ylabel(r"$F(\omega)$ [states/eV]")

        # Plot a2f
        self.plot(ax=ax_list[2], show=False)

        return fig

    @add_fig_kwargs
    def plot_tc_vs_mustar(self, start=0.1, stop=0.5, num=50, ax=None, **kwargs):
        """
        Plot Tc(mustar)

        Args:
            start: The starting value of the sequence.
            stop: The end value of the sequence
            num (int): optional. Number of samples to generate. Default is 50. Must be non-negative.
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        # TODO start and stop to avoid singularity in Mc Tc
        mustar_values = np.linspace(start, stop, num=num)
        tc_vals = [self.get_mcmillan_tc(mustar) for mustar in mustar_values]

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.plot(mustar_values, tc_vals, **kwargs)
        ax.set_yscale("log")
        ax.grid(True)
        ax.set_xlabel(r"$\mu^*$")
        ax.set_ylabel(r"$T_c [K]$")

        return fig


class A2Ftr(object):

    # Markers used for up/down bands (collinear spin)
    marker_spin = {0: "^", 1: "v"}

    def __init__(self, mesh, vals_in, vals_out):
        """

        Args:
            mesh: Energy mesh in eV
	    vals_in(nomega,3,3,0:natom3,nsppol):
		Eliashberg transport functions for in and out scattering
	    vals_in(w,3,3,1:natom3,1:nsppol): a2f_tr(w) decomposed per phonon branch and spin
	    vals_in(w,3,3,0,1:nsppol): a2f_tr(w) summed over phonons modes, decomposed in spin
        """
        self.mesh = mesh

    @lazy_property
    def iw0(self):
        """
        Index of the first point in the mesh whose value is >= 0
        Integrals are performed with wmesh[iw0 + 1, :] i.e. unstable modes are neglected.
        """
        for i, x in enumerate(self.mesh):
            if x >= 0.0: return i
        else:
            raise ValueError("Cannot find zero in energy mesh")


# TODO Change name.
class EphFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the phonon linewidths, EliashbergFunction, the |PhononBands|,
    the |ElectronBands| and |ElectronDos| on the k-mesh.
    Provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with EphFile("out_EPH.nc") as ncfile:
            print(ncfile)
            ncfile.ebands.plot()
            ncfile.phbands.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EphFile
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
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
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        app("")
        app(self.phbands.to_string(with_structure=False, verbose=verbose, title="Phonon Bands"))
        app("")
        # E-PH section
        app(marquee("E-PH calculation", mark="="))
        app("Has transport a2Ftr(w): %s" % self.has_a2ftr)
        app(self.a2f_qcoarse.to_string(title="A2f(w) on the ab-initio q-mesh:", verbose=verbose))
        app("")
        app(self.a2f_qintp.to_string(title="A2f(w) interpolated on the dense q-mesh:", verbose=verbose))

        return "\n".join(lines)

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @lazy_property
    def edos(self):
        """|ElectronDos| object with e-DOS computed by Abinit."""
        return self.reader.read_edos()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    @property
    def phbands(self):
        """|PhononBands| object with frequencies along the q-path."""
        return self.reader.read_phbands_qpath()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        return od

    @lazy_property
    def a2f_qcoarse(self):
        """
        :class:`A2f` with the Eliashberg function a2F(w)
        computed on the (coarse) ab-initio q-mesh.
        """
        return self.reader.read_a2f(qsamp="qcoarse")

    @lazy_property
    def a2f_qintp(self):
        """
        :class:`A2f` with the Eliashberg function a2F(w)
        computed on the dense q-mesh by Fourier interpolation.
        """
        return self.reader.read_a2f(qsamp="qintp")

    def get_a2f_qsamp(self, qsamp):
        """Return the :class:`A2f` object associated to q-sampling ``qsamp``."""
        if qsamp == "qcoarse": return self.a2f_qcoarse
        if qsamp == "qintp": return self.a2f_qintp
        raise ValueError("Invalid value for qsamp `%s`" % str(qsamp))

    @lazy_property
    def has_a2ftr(self):
        """True if the netcdf file contains transport data."""
        return "a2ftr_qcoarse" in self.reader.rootgrp.variables

    @lazy_property
    def a2ftr_qcoarse(self):
        """
        :class:`A2ftr` with the Eliashberg transport spectral function a2F_tr(w, x, x')
        computed on the (coarse) ab-initio q-mesh
        """
        if not self.has_a2ftr: return None
        return self.reader.read_a2ftr(qsamp="qcoarse")

    @lazy_property
    def a2ftr_qinpt(self):
        """
        :class:`A2ftr` with the Eliashberg transport spectral function a2F_tr(w, x, x')
        computed on the dense q-mesh by Fourier interpolation.
        """
        if not self.has_a2ftr: return None
        return self.reader.read_a2ftr(qsamp="qintp")

    def get_a2ftr_qsamp(self, qsamp):
        """Return the :class:`A2ftr` object associated to q-sampling ``qsamp``."""
        if qsamp == "qcoarse": return self.a2ftr_qcoarse
        if qsamp == "qintp": return self.a2ftr_qintp
        raise ValueError("Invalid value for qsamp `%s`" % str(qsamp))

    def close(self):
        """Close the file."""
        self.reader.close()

    @add_fig_kwargs
    def plot_eph_strength(self, what="lambda", ylims=None, ax=None, label=None, fontsize=12, **kwargs):
        """
        Plot phonon bands with eph coupling strength lambda(q, nu)

        Args:
            what: ``lambda`` for the eph coupling strength, ``gamma`` for phonon linewidths.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            ax: |matplotlib-Axes| or None if a new figure should be created.
            label: String used to label the plot in the legend.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Plot phonon bands
        #self.phbands.plot(ax=ax, units=units, show=False)

        # Decorate the axis (e.g add ticks and labels).
        self.phbands.decorate_ax(ax, units="")

        # Add eph coupling.
        if what == "lambda":
            yvals = self.reader.read_phlambda_qpath()[0]
            ylabel = r"$\lambda(q,\nu)$"
        elif what == "gamma":
            yvals = self.reader.read_phgamma_qpath()[0]
            ylabel = r"$\gamma(q,\nu)$ [eV]"
        else:
            raise ValueError("Invalid value for what: `%s`" % str(what))

        style = dict(
            linestyle=kwargs.pop("linestyle", "-"),
            color=kwargs.pop("color", "k"),
            linewidth=kwargs.pop("linewidth", 1),
        )

        xvals = np.arange(len(self.phbands.qpoints))
        for nu in self.phbands.branches:
            ax.plot(xvals, yvals[:, nu],
                    label=label if (nu == 0 and label) else None,
                    **style)

        ax.set_ylabel(ylabel)
        set_axlims(ax, ylims, "y")
        if label: ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot(self, what="lambda", units="eV", alpha=0.8, ylims=None, ax=None, **kwargs):
        """
        Plot phonon bands with eph coupling strength lambda(q, nu) or gamma(q, nu)

        Args:
            what: ``lambda`` for eph coupling strength, ``gamma`` for phonon linewidths.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            scale: float used to scale the marker size.
            alpha: The alpha blending value for the markers between 0 (transparent) and 1 (opaque)
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Plot phonon bands.
        self.phbands.plot(ax=ax, units=units, show=False)

        # Add eph coupling.
        xvals = np.arange(len(self.phbands.qpoints))
        yvals = self.phbands.phfreqs * abu.phfactor_ev2units(units)

        # [0] is for the number_of_spins
        # TODO units
        #if what == "lambda":
        #    s = self.reader.read_phlambda_qpath()[0]
        #    scale = 100
        #elif what == "gamma":
        #    s = self.reader.read_phgamma_qpath()[0]
        #    scale = 1
        #else:
        #    raise ValueError("Invalid value fo what: `%s`" % what)

        scale = 1
        gammas = self.reader.read_phgamma_qpath()[0]
        gam_min, gam_max = gammas.min(), gammas.max()
        lambdas = self.reader.read_phlambda_qpath()[0]
        lamb_min, lamb_max = lambdas.min(), lambdas.max()
        cmap = "jet"

        for nu in self.phbands.branches:
            """
            scale = 100000
            ax.scatter(xvals, yvals[:, nu], s=(scale * np.abs(gammas[:, nu]))**2,
                       c=lambdas[:, nu],
                       vmin=lamb_min, vmax=lamb_max,
                       cmap=cmap,
                       marker="o",
                       #c=color,
                       #alpha=alpha
                       #label=term if ib == 0 else None
            )
            """

            scale = 500
            ax.scatter(xvals, yvals[:, nu], s=scale * np.abs(lambdas[:, nu]),
                       c=gammas[:, nu],
                       vmin=gam_min, vmax=gam_max,
                       cmap=cmap,
                       marker="o",
                       #c=color,
                       alpha=alpha,
                       #label=term if ib == 0 else None
            )

        # Make a color bar
        #plt.colorbar(ax, cmap=cmap)

        set_axlims(ax, ylims, "y")
        return fig

    @add_fig_kwargs
    def plot_a2f_interpol(self, units="eV", ax=None, ylims=None, **kwargs):
        """
        Plot

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        #linestyle_qsamp = dict(qcoarse="--", qintp="-")
        for qsamp in ["qcoarse", "qintp"]:
            a2f = self.get_a2f_qsamp(qsamp)
            a2f.plot(units=units, ylims=ylims, ax=ax, with_lambda=False, show=False)
            #ax.yaxis.set_ticks_position("right")
            #ax.yaxis.set_label_position("right")
            #ax.tick_params(labelbottom='off')
            #ax.set_ylabel("")

        return fig

    @add_fig_kwargs
    def plot_with_a2f(self, units="eV", qsamp="qintp", phdos=None, ylims=None, **kwargs):
        """
        Plot phonon bands with lambda(q, nu) + a2F(w) + phonon DOS.
        """
        # Max three additional axes with [a2F, a2F_tr, DOS]
        ncols = 2
        width_ratios = [1, 0.2]
        if self.has_a2ftr:
            ncols += 1
            width_ratios.append(0.2)

        if phdos is not None:
            phdos = PhononDos.as_phdos(phdos)
            ncols += 1
            width_ratios.append(0.2)

        # Build grid plot.
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        fig = plt.figure()
        gspec = GridSpec(1, ncols, width_ratios=width_ratios, wspace=0.05)
        ax_phbands = plt.subplot(gspec[0])

        ax_doses = []
        for i in range(ncols - 1):
            ax = plt.subplot(gspec[i + 1], sharey=ax_phbands)
            ax.grid(True)
            set_axlims(ax, ylims, "y")
            ax_doses.append(ax)

        # Plot phonon bands with markers.
        self.plot(units=units, ylims=ylims, ax=ax_phbands, show=False)

        # Plot a2F(w)
        a2f = self.get_a2f_qsamp(qsamp)
        ax = ax_doses[0]
        a2f.plot(units=units, exchange_xy=True, ylims=ylims, ax=ax, with_lambda=False, show=False)
        ax.yaxis.set_ticks_position("right")
        #ax.yaxis.set_label_position("right")
        #ax.tick_params(labelbottom='off')
        ax.set_ylabel("")

        # Plot a2Ftr(w)
        ix = 1
        if self.has_a2ftr:
            ax = ax_doses[ix]
            a2ftr = self.get_a2ftr_qsamp(qsamp)
            self.a2ftr.plot(units=units, exchange_xy=True, ylims=ylims, ax=ax, show=False)
            ax.yaxis.set_ticks_position("right")
            #ax.yaxis.set_label_position("right")
            #ax.tick_params(labelbottom='off')
            ax.set_ylabel("")
            ix += 1

        # Plot DOS g(w)
        if phdos is not None:
            ax = ax_doses[ix]
            phdos.plot_dos_idos(ax=ax, exchange_xy=True, what="d", color="k", linestyle="-")
            ax.yaxis.set_ticks_position("right")
            #ax.yaxis.set_label_position("right")
            #ax.tick_params(labelbottom='off')
            ax.set_xlabel(r"$F(\omega)$")
            #ax.set_ylabel("")

        return fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.ebands.plot();"),
            nbv.new_code_cell("ncfile.plot();"),
            #nbv.new_code_cell("ncfile.plot_phlinewidths();"),
            nbv.new_code_cell("ncfile.plot_with_a2f();"),
            nbv.new_code_cell("ncfile.a2f.plot();"),
        ])

        if self.has_a2ftr:
            nb.cells.extend([
                nbv.new_code_cell("ncfile.a2ftr.plot();"),
                #nbv.new_code_cell("ncfile.plot_with_a2ftr();"),
            ])

        return self._write_nb_nbpath(nb, nbpath)


class EphRobot(Robot, RobotWithEbands, RobotWithPhbands):
    """
    This robot analyzes the results contained in multiple EPH.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EphRobot
    """
    #TODO: Method to plot the convergence of DOS(e_F)
    EXT = "EPH"

    linestyle_qsamp = dict(qcoarse="--", qintp="-")

    all_qsamps = ["qcoarse", "qintp"]

    def get_dataframe(self, abspath=False, with_geo=False, funcs=None):
        """
        Build and return a |pandas-DataFrame| with the most important results.

        Args:
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().
            with_geo: True if structure info should be added to the dataframe
            funcs: Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a :class:`EphFile` object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.

        Return: |pandas-DataFrame|
        """
        rows, row_names = [], []
        for i, (label, ncfile) in enumerate(self.items()):
            row_names.append(label)
            d = OrderedDict()

            for qsamp in self.all_qsamps:
                a2f = ncfile.get_a2f_qsamp(qsamp)
                d["lambda_" + qsamp] = a2f.lambda_iso
                d["omegalog_" + qsamp] = a2f.omega_log

                # Add transport properties.
                if ncfile.has_a2ftr:
                    for qsamp in self.all_qsamps:
                        a2ftr = ncfile.get_a2ftr_qsamp(qsamp)
                        d["lambdatr_avg_" + qsamp] = a2f.lambda_tr

            # Add convergence parameters
            #d.update(ncfile.params)

            # Add info on structure.
            if with_geo:
                d.update(ncfile.structure.get_dict4pandas(with_spglib=True))

            # Execute functions.
            if funcs is not None: d.update(self._exec_funcs(funcs, ncfile))
            rows.append(d)

        import pandas as pd
        row_names = row_names if not abspath else self._to_relpaths(row_names)
        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    @add_fig_kwargs
    def plot_lambda_convergence(self, what="lambda", sortby=None, ylims=None, fontsize=8,
                                ax=None, colormap="jet", **kwargs):
        """
        Plot the convergence of the lambda(q, nu) parameters wrt to the ``sortby`` parameter.

        Args:
            what: "lambda" for eph strength, gamma for phonon linewidths.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Legend and title fontsize.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: matplotlib color map.

        Returns: |matplotlib-Figure|
        """
        # TODO Add hue
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)
        for i, (label, ncfile, param) in enumerate(self.sortby(sortby)):
            ncfile.plot_eph_strength(
                    ax=ax,
                    what=what, ylims=ylims,
                    label=self.sortby_label(sortby, param),
                    color=cmap(i / len(self)), fontsize=fontsize,
                    show=False,
                    )
        return fig

    @add_fig_kwargs
    def plot_a2f_convergence(self, sortby=None, qsamps="all", ax=None, xlims=None,
                            fontsize=8, colormap="jet", **kwargs):
        """
        Plot the convergence of the Eliashberg function wrt to the ``sortby`` parameter.

        Args:
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            qsamps:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            fontsize: Legend and title fontsize.
            colormap: matplotlib color map.

        Returns: |matplotlib-Figure|
        """
        # TODO Add hue
        qsamps = self.all_qsamps if qsamps == "all" else list_strings(qsamps)
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)

        for i, (label, ncfile, param) in enumerate(self.sortby(sortby)):
            for qsamp in qsamps:
                ncfile.get_a2f_qsamp(qsamp).plot(
                    ax=ax,
                    label=self.sortby_label(sortby, param) + " " + qsamp,
                    color=cmap(i / len(self)), fontsize=fontsize,
                    linestyle=self.linestyle_qsamp[qsamp],
                    show=False,
                )

        set_axlims(ax, xlims, "x")
        return fig

    #@add_fig_kwargs
    #def plot_a2ftr_convergence(self, sortby=None, qsamps="all", ax=None, xlims=None,
    #                           fontsize=8, colormap="jey", **kwargs):
    #    qsamps = self.all_qsamps if qsamps == "all" else list_strings(qsamps)
    #    ax, fig, plt = get_ax_fig_plt(ax=ax)
    #    cmap = plt.get_cmap(colormap)
    #    for i, (label, ncfile, param) in enumerate(self.sortby(sortby)):
    #        for qsamp in qsamps:
    #            ncfile.get_a2ftr_qsamp(qsamp).plot(
    #                    ax=ax=ax,
    #                    label=self.sortby_label(sortby, param),
    #                    color=cmap(i / len(self)),
    #                    show=False,
    #                    )
    #    set_axlims(ax, xlims, "x")
    #    return fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.EphRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("data = robot.get_dataframe()\ndata"),
            nbv.new_code_cell("robot.plot_lambda_convergence();"),
            nbv.new_code_cell("robot.plot_a2f_convergence();"),
        ])

        if all(ncf.has_a2ftr for ncf in self.abifiles):
            nb.cells.extend([
                nbv.new_code_cell("robot.plot_a2ftr_convergence();"),
            ])

        # Mixins.
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())
        nb.cells.extend(self.get_phbands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)


class EphReader(ElectronsReader):
    """
    Reads data from the EPH.nc file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EphReader
    """
    def read_edos(self):
        """
        Read the |ElectronDos| used to compute EPH quantities.
        """
        mesh = self.read_value("edos_mesh") * units.Ha_to_eV
        # [nsppol+1, nw] arrays with TOT_DOS, Spin_up, Spin_down in a.u.
        var = self.read_variable("edos_dos")
        if var.shape[0] == 3:
            # Spin polarized. Extract up-down components.
            spin_dos = var[1:, :] / units.Ha_to_eV
        else:
            # Spin unpolarized. Extract Tot DOS
            spin_dos = var[0, :] / units.Ha_to_eV

        #spin_idos = self.read_variable("edos_idos")[1:, :] / units.Ha_to_eV
        nelect = self.read_value("number_of_electrons")
        fermie = self.read_value("fermi_energy") * units.Ha_to_eV
        return ElectronDos(mesh, spin_dos, nelect, fermie=fermie)

    def read_phbands_qpath(self):
        """
        Read and return a |PhononBands| object with frequencies computed along the q-path.
        """
        structure = self.read_structure()

        # Build the list of q-points
        qpoints = Kpath(structure.reciprocal_lattice,
                        frac_coords=self.read_value("qpath"),
                        weights=None, names=None, ksampling=None)

        #nctkarr_t('phfreq_qpath', "dp", "natom3, nqpath"),&
        phfreqs = self.read_value("phfreq_qpath") * units.Ha_to_eV
        phdispl_cart = self.read_value("phdispl_cart_qpath", cmode="c") * units.bohr_to_ang

        amu_list = self.read_value("atomic_mass_units", default=None)
        if amu_list is not None:
            atom_species = self.read_value("atomic_numbers")
            amu = {at: a for at, a in zip(atom_species, amu_list)}
        else:
            raise ValueError("atomic_mass_units is not present!")
            amu = None

        return PhononBands(structure=structure,
                           qpoints=qpoints,
                           phfreqs=phfreqs,
                           phdispl_cart=phdispl_cart,
                           non_anal_ph=None,
                           amu=amu,
                           )

    def read_phlambda_qpath(self):
        """
        Reads the EPH coupling strength *interpolated* along the q-path.
        Return |numpy-array| with shape [nsppol, nqpath, natom3]
        """
        return self.read_value("phlambda_qpath")

    def read_phgamma_qpath(self):
        """
        Reads the phonon linewidths *interpolated* along the q-path.
        Return results in eV
        Return |numpy-array| with shape [nsppol, nqpath, natom3]
        """
        return self.read_value("phgamma_qpath") * units.Ha_to_eV

    def read_a2f(self, qsamp):
        """
        Read and return the Eliashberg function :class:`A2F`.
        """
        assert qsamp in ("qcoarse", "qintp")
        mesh = self.read_value("a2f_mesh_" + qsamp)  * units.Ha_to_eV
        # C shape [nsppol, natom + 1, nomega]
        data = self.read_value("a2f_values_" + qsamp) # * 0.25
        values_spin = data[:, 0, :].copy()
        values_spin_nu = data[:, 1:, :].copy()
        return A2f(mesh, values_spin, values_spin_nu)

    #def read_a2ftr(self, qsamp):
    #    """Read and return the Eliashberg transport spectral function a2F_tr(w, x, x')."""
    #    assert qsamp in ("qcoarse", "qintp")
    #    mesh = self.read_value("a2ftr_mesh_" + qsamp) * units.Ha_to_eV
    #    # Transpose tensor components F --> C
    #    vals_in = self.read_value("a2ftr_in_" + qsamp)
    #    vals_out = self.read_value("a2ftr_out_" + qsamp)
    #    return A2ftr(mesh=mesh, vals_in, vals_out)
