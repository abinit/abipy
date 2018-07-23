# coding: utf-8
"""
This module contains objects for postprocessing A2F calculations (phonon lifetimes in metals
and Eliashberg function).

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
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels)
from abipy.tools import duck
from abipy.electrons.ebands import ElectronDos, RobotWithEbands
from abipy.dfpt.phonons import PhononBands, PhononDos, RobotWithPhbands
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader


_LATEX_LABELS = {
    "lambda_iso": r"$\lambda_{iso}$",
    "omega_log": r"$\omega_{log}$",
    "a2f": r"$\alpha^2F(\omega)$",
    "lambda": r"$\lambda(\omega)$",
}


class A2f(object):
    """
    Eliashberg function a2F(w). Energies are in eV.
    """
    # Markers used for up/down bands.
    marker_spin = {0: "^", 1: "v"}

    def __init__(self, mesh, values_spin, values_spin_nu, ngqpt, meta):
        """
        Args:
            mesh: Energy mesh in eV
            values(nomega,0:natom3,nsppol)
            vals(w,1:natom,1:nsppol): a2f(w) decomposed per phonon branch and spin
            vals(w,0,1:nsppol): a2f(w) summed over phonons modes, decomposed in spin
            ngqpt: Q-mesh used to compute A2f.
            meta: Dictionary with metavariables.

        TODO:
            1. possibility of computing a2f directly from data on file?
        """
        self.mesh = mesh
        self.ngqpt = ngqpt
        self.meta = meta

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
        """
        String representation with verbosity level ``verbose`` and an optional ``title``.
        """
        lines = []; app = lines.append

        app("Eliashberg Function" if not title else str(title))
        # TODO: Add ElectronDos
        #app("Isotropic lambda: %.3f" % (self.lambda_iso))
        app("Isotropic lambda: %.2f, omega_log: %.3f (eV), %.3f (K)" % (self.lambda_iso, self.omega_log, self.omega_log * abu.eV_to_K))
        app("Q-mesh: %s" % str(self.ngqpt))
        app("Mesh from %.4f to %.4f (eV) with %d points" % (
            self.mesh[0], self.mesh[-1], len(self.mesh)))

        if verbose:
            for mustar in (0.1, 0.12, 0.2):
                app("\tFor mustar %s: McMillan Tc: %s [K]" % (mustar, self.get_mcmillan_tc(mustar)))
        if verbose > 1:
            # $\int dw [a2F(w)/w] w^n$
            for n in [0, 4]:
                app("Moment %s: %s" % (n, self.get_moment(n)))

            app("Meta: %s" % str(self.meta))

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
        iw = self.iw0 + 1
        wmesh, a2fw = self.mesh[iw:], self.values[iw:]

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
    def plot(self, what="a2f", units="eV", exchange_xy=False, ax=None,
             xlims=None, ylims=None, label=None, fontsize=12, **kwargs):
        """
        Plot a2F(w) or lambda(w) depending on the value of `what`.

        Args:
            what: a2f for a2F(w), lambda for lambda(w)
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            exchange_xy: True to exchange x-y axes.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
		or scalar e.g. ``left``. If left (right) is None, default values are used
            ylims: Limits for y-axis. See xlims for API.
            label: True to add legend label to each curve.
            fontsize: Legend and title fontsize
            kwargs: linestyle, color, linewidth passed to ax.plot.

        Returns: |matplotlib-Figure|
        """""
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        wfactor = abu.phfactor_ev2units(units)
        ylabel = _LATEX_LABELS[what]

        style = dict(
            linestyle=kwargs.pop("linestyle", "-"),
            color=kwargs.pop("color", "k"),
            linewidth=kwargs.pop("linewidth", 1),
        )

        # Plot a2f(w)
        if what == "a2f":
            xx, yy = self.mesh * wfactor, self.values
            if exchange_xy: xx, yy = yy, xx
            ax.plot(xx, yy, label=label, **style)

            if self.nsppol == 2:
                # Plot spin resolved a2f(w).
                for spin in range(self.nsppol):
                    xx, yy = self.mesh * wfactor, self.values_spin[spin]
                    if exchange_xy: xx, yy = yy, xx
                    ax.plot(xx, yy, marker=self.marker_spin[spin], **style)

        # Plot lambda(w)
        elif what == "lambda":
            lambda_w = self.get_moment(n=0, cumulative=True)
            xx, yy = self.mesh * wfactor, lambda_w
            if exchange_xy: xx, yy = yy, xx
            ax.plot(xx, yy, label=label, **style)

        else:
            raise ValueError("Invalid value for what: `%s`" % str(what))

        xlabel = abu.wlabel_from_units(units)
        if exchange_xy: xlabel, ylabel = ylabel, xlabel

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True)
        set_axlims(ax, xlims, "x")
        set_axlims(ax, ylims, "y")
        if label: ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_with_lambda(self, units="eV", ax=None, xlims=None, fontsize=12, **kwargs):
        """
        Plot a2F(w) and lambda(w) on the same figure.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            xlims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
		or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Legend and title fontsize

        Returns: |matplotlib-Figure|
        """""
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for i, what in enumerate(["a2f", "lambda"]):
            this_ax = ax if i == 0 else ax.twinx()
            self.plot(what=what, ax=this_ax, units=units, fontsize=fontsize, xlims=xlims, show=False, **kwargs)
            if i:
                this_ax.yaxis.set_label_position("right")
                this_ax.grid(True)

        return fig

    @add_fig_kwargs
    def plot_nuterms(self, units="eV", ax_mat=None, with_lambda=True, fontsize=12,
                     xlims=None, ylims=None, label=None, **kwargs):
        """
        Plot a2F(w), lambda(w) and optionally the individual contributions due to the phonon branches.

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
        Grid with 3 plots showing: a2F(w), F(w), a2F(w). Requires phonon DOS.

        Args:
            phdos: |PhononDos|
            atol: F(w) is replaced by atol in a2F(w) / F(w) ratio where :math:`|F(w)|` < atol

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

        # Plot F(w). TODO: This should not be called plot_dos_idos!
        ax = ax_list[1]
        phdos.plot_dos_idos(ax=ax, what="d", color="k", linestyle="-")
        ax.grid(True)
        ax.set_ylabel(r"$F(\omega)$ [states/eV]")

        # Plot a2f
        self.plot(ax=ax_list[2], color="k", linestyle="-", linewidths=2, show=False)

        return fig

    @add_fig_kwargs
    def plot_tc_vs_mustar(self, start=0.1, stop=0.3, num=50, ax=None, **kwargs):
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
        ax.set_ylabel(r"$T_c$ [K]")

        return fig


class A2Ftr(object):
    """
    Transport Eliashberg function a2F(w). Energies are in eV.
    """
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


class A2fFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the phonon linewidths, EliashbergFunction, the |PhononBands|,
    the |ElectronBands| and |ElectronDos| on the k-mesh.
    Provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with A2fFile("out_A2F.nc") as ncfile:
            print(ncfile)
            ncfile.ebands.plot()
            ncfile.phbands.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: A2fFile
    """
    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf_ file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(A2fFile, self).__init__(filepath)
        self.reader = A2fReader(filepath)

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
        app("K-mesh for electrons:")
        app(self.ebands.kpoints.ksampling.to_string(verbose=verbose))
        if verbose:
            app("Has transport a2Ftr(w): %s" % self.has_a2ftr)
        app("")
        a2f = self.a2f_qcoarse
        app("a2f(w) on the %s q-mesh (ddb_ngqpt|eph_ngqpt)" % str(a2f.ngqpt))
        app("Isotropic lambda: %.2f, omega_log: %.3f (eV), %.3f (K)" % (
            a2f.lambda_iso, a2f.omega_log, a2f.omega_log * abu.eV_to_K))
        #app(self.a2f_qcoarse.to_string(title=title, verbose=verbose))
        app("")
        a2f = self.a2f_qintp
        app("a2f(w) Fourier interpolated on the %s q-mesh (ph_ngqpt)" % str(a2f.ngqpt))
        app("Isotropic lambda: %.2f, omega_log: %.3f (eV), %.3f (K)" % (
            a2f.lambda_iso, a2f.omega_log, a2f.omega_log * abu.eV_to_K))
        #app(self.a2f_qintp.to_string(title=title, verbose=verbose))

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
        """
        |PhononBands| object with frequencies along the q-path.
        Contains (interpolated) linewidths.
        """
        return self.reader.read_phbands_qpath()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        od = self.get_ebands_params()
        # Add EPH parameters.
        od.update(self.reader.common_eph_params)

        return od

    @lazy_property
    def a2f_qcoarse(self):
        """
        :class:`A2f` with the Eliashberg function a2F(w) computed on the (coarse) ab-initio q-mesh.
        """
        return self.reader.read_a2f(qsamp="qcoarse")

    @lazy_property
    def a2f_qintp(self):
        """
        :class:`A2f` with the Eliashberg function a2F(w) computed on the dense q-mesh by Fourier interpolation.
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
    def a2ftr_qintp(self):
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

    #def interpolate(self, ddb, lpratio=5, vertices_names=None, line_density=20, filter_params=None, verbose=0):
    #    """
    #    Interpolate the phonon linewidths on a k-path and, optionally, on a k-mesh.

    #    Args:
    #        lpratio: Ratio between the number of star functions and the number of ab-initio k-points.
    #            The default should be OK in many systems, larger values may be required for accurate derivatives.
    #        vertices_names: Used to specify the k-path for the interpolated QP band structure
    #            when ``ks_ebands_kpath`` is None.
    #            It's a list of tuple, each tuple is of the form (kfrac_coords, kname) where
    #            kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
    #            the k-point. Each point represents a vertex of the k-path. ``line_density`` defines
    #            the density of the sampling. If None, the k-path is automatically generated according
    #            to the point group of the system.
    #        line_density: Number of points in the smallest segment of the k-path. Used with ``vertices_names``.
    #        filter_params: TO BE DESCRIBED
    #        verbose: Verbosity level

    #    Returns:
    #    """
    #    # Get symmetries from abinit spacegroup (read from file).
    #    abispg = self.structure.abi_spacegroup
    #    fm_symrel = [s for (s, afm) in zip(abispg.symrel, abispg.symafm) if afm == 1]

    #    phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(nqsmall=0, ndivsm=10, asr=2, chneut=1, dipdip=1,
    #        dos_method="tetra", lo_to_splitting="automatic", ngqpt=None, qptbounds=None, anaddb_kwargs=None, verbose=0,
    #        mpi_procs=1, workdir=None, manager=None)

    #    phbands = phbst_file.phbands
    #    phbst_file.close()

    #    # Read qibz and ab-initio linewidths from file.
    #    qcoords_ibz = self.reader.read_value("qibz")
    #    data_ibz = self.reader.read_value("phgamma_qibz") * units.Ha_to_eV
    #    import matplotlib.pyplot as plt
    #    plt.plot(data_ibz[0])
    #    plt.show()

    #    # Build interpolator.
    #    from abipy.core.skw import SkwInterpolator
    #    cell = (self.structure.lattice.matrix, self.structure.frac_coords, self.structure.atomic_numbers)

    #    has_timrev = True
    #    fermie, nelect = 0.0, 3 * len(self.structure)
    #    skw = SkwInterpolator(lpratio, qcoords_ibz, data_ibz, fermie, nelect,
    #                          cell, fm_symrel, has_timrev,
    #                          filter_params=filter_params, verbose=verbose)

    #    # Interpolate and set linewidths.
    #    qfrac_coords = [q.frac_coords for q in phbands.qpoints]
    #    phbands.linewidths = skw.interp_kpts(qfrac_coords).eigens

    #    return phbands

    @add_fig_kwargs
    def plot_eph_strength(self, what_list=("phbands", "gamma", "lambda"), ax_list=None,
                          ylims=None, label=None, fontsize=12, **kwargs):
        """
        Plot phonon bands with EPH coupling strength lambda(q, nu) and lambda(q, nu)
        These values have been Fourier interpolated by Abinit.

        Args:
            what_list: ``phfreqs`` for phonons, `lambda`` for the eph coupling strength,
                ``gamma`` for phonon linewidths.
            ax_list: List of |matplotlib-Axes| (same length as what_list)
                or None if a new figure should be created.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            label: String used to label the plot in the legend.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        what_list = list_strings(what_list)
        nrows, ncols = len(what_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()
        units = "eV"

        for i, (ax, what) in enumerate(zip(ax_list, what_list)):
            # Decorate the axis (e.g add ticks and labels).
            self.phbands.decorate_ax(ax, units="")

            if what == "phbands":
                # Plot phonon bands
                self.phbands.plot(ax=ax, units=units, show=False)
            else:
                # Add eph coupling.
                if what == "lambda":
                    yvals = self.reader.read_phlambda_qpath()
                    ylabel = r"$\lambda(q,\nu)$"
                elif what == "gamma":
                    yvals = self.reader.read_phgamma_qpath()
                    ylabel = r"$\gamma(q,\nu)$ (eV)"
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
    def plot(self, what="gamma", units="eV", scale=None, alpha=0.6, ylims=None, ax=None, colormap="jet", **kwargs):
        """
        Plot phonon bands with gamma(q, nu) or lambda(q, nu) depending on the vaue of `what`.

        Args:
            what: ``lambda`` for eph coupling strength, ``gamma`` for phonon linewidths.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            scale: float used to scale the marker size.
            alpha: The alpha blending value for the markers between 0 (transparent) and 1 (opaque)
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            ax: |matplotlib-Axes| or None if a new figure should be created.
            colormap: matplotlib color map.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)

        # Plot phonon bands.
        self.phbands.plot(ax=ax, units=units, show=False)

        # Add eph coupling.
        xvals = np.arange(len(self.phbands.qpoints))
        wvals = self.phbands.phfreqs * abu.phfactor_ev2units(units)

        # Sum contributions over nsppol (if spin-polarized)
        # TODO units
        gammas = self.reader.read_phgamma_qpath()
        lambdas = self.reader.read_phlambda_qpath()
        if what == "lambda":
            scale = 500 if scale is None else float(scale)
            sqn = scale * np.abs(lambdas)
            cqn = gammas
        elif what == "gamma":
            scale = 10 ** 6 if scale is None else float(scale)
            sqn = scale * np.abs(gammas)
            cqn = lambdas
        else:
            raise ValueError("Invalid what: `%s`" % str(what))
        vmin, vmax = cqn.min(), cqn.max()

        sc = ax.scatter(np.tile(xvals, len(self.phbands.branches)),
                        wvals.T, # [q, nu] --> [nu, q]
                        s=sqn.T,
                        c=cqn.T,
                        vmin=vmin, vmax=vmax,
                        cmap=cmap,
                        marker="o",
                        alpha=alpha,
                        #label=term if ib == 0 else None
        )

        # Make a color bar
        #plt.colorbar(sc, ax=ax, orientation="horizontal", pad=0.2)
        set_axlims(ax, ylims, "y")

        return fig

    @add_fig_kwargs
    def plot_a2f_interpol(self, units="eV", ylims=None, fontsize=8, **kwargs):
        """
        Compare ab-initio a2F(w) with interpolated values.

        Args:
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz").
                Case-insensitive.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Legend and title fontsize

        Returns: |matplotlib-Figure|
        """
        what_list = ["a2f", "lambda"]
        nrows, ncols = len(what_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        styles = dict(
            qcoarse={"linestyle": "--", "color": "b"},
            qintp={"linestyle": "-", "color": "r"},
        )

        for ix, (ax, what) in enumerate(zip(ax_list, what_list)):
            for qsamp in ["qcoarse", "qintp"]:
                a2f = self.get_a2f_qsamp(qsamp)
                a2f.plot(what=what, ax=ax, units=units, ylims=ylims, fontsize=fontsize,
                         label=qsamp if ix == 0 else None,
                         show=False, **styles[qsamp])

        return fig

    @add_fig_kwargs
    def plot_with_a2f(self, what="gamma", units="eV", qsamp="qintp", phdos=None, ylims=None, **kwargs):
        """
        Plot phonon bands with lambda(q, nu) + a2F(w) + phonon DOS.

        Args:
            what: ``lambda`` for eph coupling strength, ``gamma`` for phonon linewidths.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            qsamp:
            phdos: |PhononDos| object. Used to plot the PhononDos on the right.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
		or scalar e.g. ``left``. If left (right) is None, default values are used

        Returns: |matplotlib-Figure|
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
        self.plot(what=what, units=units, ylims=ylims, ax=ax_phbands, show=False)

        # Plot a2F(w)
        a2f = self.get_a2f_qsamp(qsamp)
        ax = ax_doses[0]
        a2f.plot(units=units, exchange_xy=True, ylims=ylims, ax=ax, show=False)
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

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        yield self.plot(show=False)
        #yield self.plot_eph_strength(show=False)
        yield self.plot_with_a2f(show=False)
        for qsamp in ["qcoarse", "qintp"]:
            a2f = self.get_a2f_qsamp(qsamp)
            yield a2f.plot_with_lambda(show=False)
        #yield self.plot_nuterms(show=False)
        #yield self.plot_a2(show=False)
        #yield self.plot_tc_vs_mustar(show=False)

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


class A2fRobot(Robot, RobotWithEbands, RobotWithPhbands):
    """
    This robot analyzes the results contained in multiple A2F.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: A2fRobot
    """
    #TODO: Method to plot the convergence of DOS(e_F)
    EXT = "A2F"

    linestyle_qsamp = dict(qcoarse="--", qintp="-")
    marker_qsamp = dict(qcoarse="^", qintp="o")

    all_qsamps = ["qcoarse", "qintp"]

    def get_dataframe(self, abspath=False, with_geo=False, with_params=True, funcs=None):
        """
        Build and return a |pandas-DataFrame| with the most important results.

        Args:
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().
            with_geo: True if structure info should be added to the dataframe
            funcs: Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a :class:`A2fFile` object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.
            with_params: False to exclude calculation parameters from the dataframe.

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

            # Add info on structure.
            if with_geo:
                d.update(ncfile.structure.get_dict4pandas(with_spglib=True))

            # Add convergence parameters
            if with_params:
                d.update(ncfile.params)

            # Execute functions.
            if funcs is not None: d.update(self._exec_funcs(funcs, ncfile))
            rows.append(d)

        import pandas as pd
        row_names = row_names if not abspath else self._to_relpaths(row_names)
        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    @add_fig_kwargs
    def plot_lambda_convergence(self, what="lambda", sortby=None, hue=None, ylims=None, fontsize=8,
                                colormap="jet", **kwargs):
        """
        Plot the convergence of the lambda(q, nu) parameters wrt to the ``sortby`` parameter.

        Args:
            what: "lambda" for eph strength, gamma for phonon linewidths.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Legend and title fontsize.
            colormap: matplotlib color map.

        Returns: |matplotlib-Figure|
        """
        # Build (1, ngroups) grid plot.
        if hue is None:
            labels_ncfiles_params = self.sortby(sortby, unpack=False)
            nrows, ncols = 1, 1
        else:
            groups = self.group_and_sortby(hue, sortby)
            nrows, ncols = 1, len(groups)

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        cmap = plt.get_cmap(colormap)

        if hue is None:
            # Plot all results on the same figure with different color.
            for i, (label, ncfile, param) in enumerate(labels_ncfiles_params):
                ncfile.plot_eph_strength(what_list=what,
                        ax_list=[ax_mat[0, 0]],
                        ylims=ylims,
                        label=self.sortby_label(sortby, param),
                        color=cmap(i / len(self)), fontsize=fontsize,
                        show=False,
                        )
        else:
            # ngroup figures
            for ig, g in enumerate(groups):
                ax = ax_mat[0, ig]
                label = "%s: %s" % (self._get_label(hue), g.hvalue)
                for ifile, ncfile in enumerate(g.abifiles):
                    ncfile.plot_eph_strength(what_list=what,
                        ax_list=[ax],
                        ylims=ylims,
                        label=label,
                        color=cmap(ifile / len(g)), fontsize=fontsize,
                        show=False,
                        )
                if ig != 0:
                    set_visible(ax, False, "ylabel")

        return fig

    @add_fig_kwargs
    def plot_a2f_convergence(self, sortby=None, hue=None, qsamps="all", xlims=None,
                            fontsize=8, colormap="jet", **kwargs):
        """
        Plot the convergence of the Eliashberg function wrt to the ``sortby`` parameter.

        Args:
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            qsamps:
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            fontsize: Legend and title fontsize.
            colormap: matplotlib color map.

        Returns: |matplotlib-Figure|
        """
        qsamps = self.all_qsamps if qsamps == "all" else list_strings(qsamps)
        #qsamps = ["qcoarse"]

        # Build (2, ngroups) grid plot.
        if hue is None:
            labels_ncfiles_params = self.sortby(sortby, unpack=False)
            nrows, ncols = len(qsamps), 1
        else:
            groups = self.group_and_sortby(hue, sortby)
            nrows, ncols = len(qsamps), len(groups)

        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)
        cmap = plt.get_cmap(colormap)

        for i, qsamp in enumerate(qsamps):
            if hue is None:
                ax = ax_mat[i, 0]
                for j, (label, ncfile, param) in enumerate(labels_ncfiles_params):
                    ncfile.get_a2f_qsamp(qsamp).plot(what="a2f", ax=ax,
                       label=self.sortby_label(sortby, param) + " " + qsamp,
                       color=cmap(j / len(self)), fontsize=fontsize,
                       linestyle=self.linestyle_qsamp[qsamp],
                       show=False,
                    )
                set_axlims(ax, xlims, "x")
            else:
                for ig, g in enumerate(groups):
                    ax = ax_mat[i, ig]
                    label = "%s: %s" % (self._get_label(hue), g.hvalue) + " " + qsamp
                    for ncfile in g.abifiles:
                        ncfile.get_a2f_qsamp(qsamp).plot(what="a2f", ax=ax,
                            label=label,
                            color=cmap(ig / len(g)), fontsize=fontsize,
                            linestyle=self.linestyle_qsamp[qsamp],
                            show=False,
                        )
                    set_axlims(ax, xlims, "x")
                    if ig != 0:
                        set_visible(ax, False, "ylabel")

                if i != len(qsamps) - 1:
                    set_visible(ax, False, "xlabel")

        return fig

    @add_fig_kwargs
    def plot_a2fdata_convergence(self, sortby=None, hue=None, qsamps="all", what_list=("lambda_iso", "omega_log"),
                                 fontsize=8, **kwargs):
        """
        Plot the convergence of the isotropic lambda and omega_log wrt the ``sortby`` parameter.

        Args:
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            qsamps:
            what_list:
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        what_list = list_strings(what_list)

        # Build grid with (n, 1) plots.
        nrows, ncols = len(what_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        if hue is None:
            labels, ncfiles, params = self.sortby(sortby, unpack=True)
        else:
            groups = self.group_and_sortby(hue, sortby)

        qsamps = self.all_qsamps if qsamps == "all" else list_strings(qsamps)
        marker = kwargs.pop("marker", "o")

        for ix, (ax, what) in enumerate(zip(ax_list, what_list)):
            #ax.set_title(what, fontsize=fontsize)
            if hue is None:
                params_are_string = duck.is_string(params[0])
                xvals = params if not params_are_string else range(len(params))
                for iq, qsamp in enumerate(qsamps):
                    a2f_list = [ncfile.get_a2f_qsamp(qsamp) for ncfile in ncfiles]
                    yvals = [getattr(a2f, what) for a2f in a2f_list]
                    l = ax.plot(xvals, yvals,
                                marker=self.marker_qsamp[qsamp],
                                linestyle=self.linestyle_qsamp[qsamp],
                                color=None if iq == 0 else l[0].get_color(),
                                )
                    if params_are_string:
                        ax.set_xticks(xvals)
                        ax.set_xticklabels(params, fontsize=fontsize)
            else:
                for g in groups:
                    for iq, qsamp in enumerate(qsamps):
                        a2f_list = [ncfile.get_a2f_qsamp(qsamp) for ncfile in g.abifiles]
                        yvals = [getattr(a2f, what) for a2f in a2f_list]
                        label = "%s: %s" % (self._get_label(hue), g.hvalue) if iq == 0 else None
                        l = ax.plot(g.xvalues, yvals, label=label,
                                    marker=self.marker_qsamp[qsamp],
                                    linestyle=self.linestyle_qsamp[qsamp],
                                    color=None if iq == 0 else l[0].get_color(),
                                    )

            ax.grid(True)
            ax.set_ylabel(_LATEX_LABELS[what])
            if ix == len(what_list) - 1:
                ax.set_xlabel("%s" % self._get_label(sortby))
                if sortby is None: rotate_ticklabels(ax, 15)
            if hue is not None:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def gridplot_a2f(self, xlims=None, fontsize=8, sharex=True, sharey=True, **kwargs):
        """
        Plot grid with a2F(w) and lambda(w) for all files treated by the robot.

        Args:
            xlims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
		or scalar e.g. ``left``. If left (right) is None, default values are used
            sharex, sharey: True to share x- and y-axis.
            fontsize: Legend and title fontsize
        """
        return self._gridplot_a2f_what("a2f", xlims=xlims, fontsize=fontsize, sharex=sharex, sharey=sharey, **kwargs)

    #@add_fig_kwargs
    #def gridplot_a2ftr(self, xlims=None, fontsize=8, sharex=True, sharey=True, **kwargs):
    #    return self._gridplot_a2f_what("a2ftr", xlims=xlims, fontsize=fontsize, sharex=sharex, sharey=sharey, **kwargs)

    def _gridplot_a2f_what(self, what, qsamps="all", xlims=None, fontsize=8, sharex=True, sharey=True, **kwargs):
        """Internal method to plot a2F or a2f_tr"""
        nrows, ncols, nplots = 1, 1, len(self)
        if nplots > 1:
            ncols = 2
            nrows = nplots // ncols + nplots % ncols

        # Build grid plot
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=sharex, sharey=sharey, squeeze=False)
        ax_list = ax_list.ravel()
        # don't show the last ax if nplots is odd.
        if nplots % ncols != 0: ax_list[-1].axis("off")

        qsamps = self.all_qsamps if qsamps == "all" else list_strings(qsamps)
        #qsamps = ["qcoarse"]
        for qsamp in qsamps:
            if what == "a2f":
                a2f_list = [ncfile.get_a2f_qsamp(qsamp) for ncfile in self.abifiles]
            elif what == "a2ftr":
                a2f_list = self.get_a2ftr_qsamp(qsamp)
            else:
                raise ValueError("Invalid value for what: `%s`" % what)

            a2f_list = [ncfile.get_a2f_qsamp(qsamp) for ncfile in self.abifiles]

            for i, (a2f, ax, title) in enumerate(zip(a2f_list, ax_list, self.keys())):
                irow, icol = divmod(i, ncols)
                # FIXME: Twinx is problematic
                a2f.plot_with_lambda(ax=ax, show=False,
                                     linestyle=self.linestyle_qsamp[qsamp],
                                     )

                set_axlims(ax, xlims, "x")
                ax.set_title(title, fontsize=fontsize)
                if (irow, icol) != (0, 0):
                    set_visible(ax, False, "ylabel")
                if irow != nrows - 1:
                    set_visible(ax, False, "xlabel")

        return fig

    #@add_fig_kwargs
    #def plot_a2ftr_convergence(self, sortby=None, qsamps="all", ax=None, xlims=None,
    #                           fontsize=8, colormap="jet", **kwargs):
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

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot_lambda_convergence(show=False)
        yield self.plot_a2f_convergence(show=False)
        yield self.plot_a2fdata_convergence(show=False)
        yield self.gridplot_a2f(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.A2fRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
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


class A2fReader(BaseEphReader):
    """
    Reads data from the EPH.nc file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: A2fReader
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

        linewidths = self.read_phgamma_qpath()
        if self.read_nsppol() == 2:
            # We have spin-resolved linewidths, sum over spins here.
            linewidths = linewidths.sum(axis=0)

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
                           linewidths=linewidths,
                           )

    def read_phlambda_qpath(self, sum_spin=True):
        """
        Reads the EPH coupling strength *interpolated* along the q-path.

        Return:
            |numpy-array| with shape [nqpath, natom3] if not sum_spin else [nsppol, nqpath, natom3]
        """
        vals = self.read_value("phlambda_qpath")
        return vals if not sum_spin else vals.sum(axis=0)

    def read_phgamma_qpath(self, sum_spin=True):
        """
        Reads the phonon linewidths (eV) *interpolated* along the q-path.

        Return:
            |numpy-array| with shape [nqpath, natom3] if not sum_spin else [nsppol, nqpath, natom3]
        """
        vals = self.read_value("phgamma_qpath") * units.Ha_to_eV
        return vals if not sum_spin else vals.sum(axis=0)

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

        # Extract q-mesh and meta variables.
        ngqpt = self.ngqpt if qsamp == "qcoarse" else self.ph_ngqpt
        meta = {k: self.common_eph_params[k] for k in
                ["eph_intmeth", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie"]}

        return A2f(mesh, values_spin, values_spin_nu, ngqpt, meta)

    #def read_a2ftr(self, qsamp):
    #    """Read and return the Eliashberg transport spectral function a2F_tr(w, x, x')."""
    #    assert qsamp in ("qcoarse", "qintp")
    #    mesh = self.read_value("a2ftr_mesh_" + qsamp) * units.Ha_to_eV
    #    # Transpose tensor components F --> C
    #    vals_in = self.read_value("a2ftr_in_" + qsamp)
    #    vals_out = self.read_value("a2ftr_out_" + qsamp)
    #    return A2ftr(mesh=mesh, vals_in, vals_out)

    #def read_phgamma_ibz_data(self):
    #     ! linewidths in IBZ
    #     nctkarr_t('qibz', "dp", "number_of_reduced_dimensions, nqibz"), &
    #     nctkarr_t('wtq', "dp", "nqibz"), &
    #     nctkarr_t('phfreq_qibz', "dp", "natom3, nqibz"), &
    #     nctkarr_t('phdispl_cart_qibz', "dp", "two, natom3, natom3, nqibz"), &
    #     nctkarr_t('phgamma_qibz', "dp", "natom3, nqibz, number_of_spins"), &
    #     nctkarr_t('phlambda_qibz', "dp", "natom3, nqibz, number_of_spins") &
