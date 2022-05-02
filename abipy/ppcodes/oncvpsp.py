# coding: utf-8
"""
Classes and functions for parsing the ONCVPSP output file and plotting the results.
"""
from __future__ import annotations

import io
import os
import abc
import re
import tempfile
import numpy as np

from collections import namedtuple
from typing import Any, Union, List, Optional
from monty.functools import lazy_property
from monty.collections import AttrDict, dict2namedtuple
from monty.os.path import which
from monty.termcolor import cprint
from abipy.core.atom import NlkState, RadialFunction, RadialWaveFunction, l2char
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.tools.derivatives import finite_diff


def is_integer(s: Any) -> bool:
    """True if object `s` in an integer."""
    try:
        c = float(s)
        return int(c) == c
    except (ValueError, TypeError):
        return False


class OncvPlotter:
    """
    Plots the results produced by a pseudopotential generator.
    """
    # TODO: Improve support for fully-relativistic case.

    # matplotlib options.
    linewidth, markersize = 2, 2

    linestyle_aeps = dict(ae="solid", ps="dashed")

    markers_aeps = dict(ae=".", ps="o")

    color_l = {0: "black",
               1: "red",
               -1: "magenta",
               2: "blue",
               -2: "cyan",
               3: "orange",
               -3: "yellow",
    }

    @classmethod
    def from_file(cls, filepath: str) -> OncvPlotter:
        """Build the plotter from an output file."""
        parser = OncvOutputParser(filepath)
        parser.scan()
        if not parser.run_completed:
            raise RuntimeError("oncvpsp output is not completed")

        return cls(parser)

    def __init__(self, parser):
        self.parser = parser

    @staticmethod
    def decorate_ax(ax, xlabel=None, ylabel=None, title=None, fontsize=8):
        """
        Decorate a `matplotlib` Axis adding xlabel, ylabel, title, grid and legend
        """
        if title: ax.set_title(title, fontsize=fontsize)
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
        ax.grid(True)
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

    def _mpl_opts_laeps(self, l: int, aeps: str) -> dict:
        """
        Return dict with matplotlib ax.plot options to plot AE/PS quantities depending on l.
        """
        return dict(
            color=self.color_l[l],
            linestyle=self.linestyle_aeps[aeps],
            #marker=self.markers_aeps[aeps],
            linewidth=self.linewidth,
            markersize=self.markersize
        )

    #def _mpl_opts_nlk(self, nlk) -> dict:

    #@property
    #def has_scattering_states(self):

    def _add_rc_vlines(self, ax):
        """
        Add vertical lines to axis `ax` showing the core radii.
        """
        for l, rc in self.parser.rc_l.items():
            ax.axvline(rc, lw=1, color=self.color_l[l], ls="--")

    @add_fig_kwargs
    def plot_atan_logders(self, ax=None, with_xlabel=True, fontsize: int = 8, **kwargs):
        """
        Plot arctan of logder on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ae, ps = self.parser.atan_logders.ae, self.parser.atan_logders.ps
        ax, fig, plt = get_ax_fig_plt(ax)

        # Note that l can be negative if FR pseudo.
        # This corresponds to ikap 2 in Fortran.

        for l, ae_alog in ae.items():
            ps_alog = ps[l]

            if not self.parser.relativistic:
                lch = f"${l2char[abs(l)]}$"
            else:
                lch = l2char[abs(l)]
                if l >= 0: lch = f"${l2char[abs(l)]}^+$"
                if l < 0: lch  = f"${l2char[abs(l)]}^-$"

            # Add pad to avoid overlapping curves.
            pad = (abs(l) + 1) * 1.0

            ae_line, = ax.plot(ae_alog.energies, ae_alog.values + pad,
                               label=f"AE {lch}",
                               **self._mpl_opts_laeps(l, "ae"))

            ps_line, = ax.plot(ps_alog.energies, ps_alog.values + pad,
                               label=f"PS {lch}",
                               **self._mpl_opts_laeps(l, "ps"))


        xlabel = "Energy (Ha)" if with_xlabel else ""
        self.decorate_ax(ax, xlabel=xlabel, ylabel="ATAN(LogDer)", title="",
                         fontsize=fontsize,
                         )

        return fig

    @add_fig_kwargs
    def plot_radial_wfs(self, ax=None, what="bound_states", fontsize: int = 8, **kwargs):
        """
        Plot AE and PS radial wavefunctions on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            what: "bound_states" or "scattering_states".
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        if what == "bound_states":
            ae_wfs, ps_wfs = self.parser.radial_wfs.ae, self.parser.radial_wfs.ps
        elif what == "scattering_states":
            ae_wfs, ps_wfs = self.parser.scattering_wfs.ae, self.parser.scattering_wfs.ps
        else:
            raise ValueError(f"Invalid value for what: {what}")

        for nlk, ae_wf in ae_wfs.items():
            ps_wf, l, k = ps_wfs[nlk], nlk.l, nlk.k

            ax.plot(ae_wf.rmesh, ae_wf.values, label=f"AE {nlk.latex}",
                    **self._mpl_opts_laeps(l, "ae"))
            ax.plot(ps_wf.rmesh, ps_wf.values, label=f"PS {nlk.latex}",
                    **self._mpl_opts_laeps(l, "ps"))

        self.decorate_ax(ax, xlabel="r (Bohr)", ylabel=r"$\phi(r)$",
                         title="Wave Functions" if what == "bound_states" else "Scattering States",
                         fontsize=fontsize,
                         )

        self._add_rc_vlines(ax)

        return fig

    @add_fig_kwargs
    def plot_projectors(self, ax=None, fontsize: int = 8, **kwargs):
        """
        Plot projectors on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        linestyle_n = {1: "solid", 2: "dashed", 3: "dotted", 4: "dashdot"}

        for nlk, proj in self.parser.projectors.items():
            ax.plot(proj.rmesh, proj.values,
                    color=self.color_l.get(nlk.l, 'black'),
                    linestyle=linestyle_n[nlk.n],
                    linewidth=self.linewidth,
                    markersize=self.markersize,
                    label=f"Proj {nlk.n}, l={nlk.latex_l}",
                    )

        self.decorate_ax(ax, xlabel="r (Bohr)", ylabel="$p(r)$", title="Projectors",
                         fontsize=fontsize,
                         )

        self._add_rc_vlines(ax)

        return fig

    @add_fig_kwargs
    def plot_densities(self, ax=None, timesr2=False, fontsize: int = 8, **kwargs):
        """
        Plot AE, PS and model densities on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        for name, rho in self.parser.densities.items():
            d = rho.values if not timesr2 else rho.values * rho.rmesh ** 2
            line, = ax.plot(rho.rmesh, d, label=name,
                            linewidth=self.linewidth, markersize=self.markersize)

        ylabel = "$n(r)$" if not timesr2 else "$r^2 n(r)$"
        self.decorate_ax(ax, xlabel="r (Bohr)", ylabel=ylabel,
                         title="Charge densities",
                         fontsize=fontsize,
                         )

        return fig

    @add_fig_kwargs
    def plot_der_densities(self, ax=None, order=1, fontsize=8, **kwargs):
        """
        Plot the radial derivatives of the densities on axis ax.
        Used to analyze possible discontinuities or strong oscillations in r-space.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        from scipy.interpolate import UnivariateSpline

        for name, rho in self.parser.densities.items():
            # Only model core charge is shown.
            if name != "rhoM": continue

            # Need linear mesh for finite_difference --> Spline input densities on lin_rmesh
            lin_rmesh, h = np.linspace(rho.rmesh[0], rho.rmesh[-1], num=len(rho.rmesh) * 4, retstep=True)
            spline = UnivariateSpline(rho.rmesh, rho.values, s=0)
            lin_values = spline(lin_rmesh)
            vder = finite_diff(lin_values, h, order=order, acc=4)
            ax.plot(lin_rmesh, vder, label="%s-order derivative of %s" % (order, name))

        self.decorate_ax(ax, xlabel="r (Bohr)", ylabel="$D^%s \n(r)$" % order,
                         title="Derivative of the charge densities",
                         fontsize=fontsize,
                         )

        return fig

    @add_fig_kwargs
    def plot_potentials(self, ax=None, fontsize: int = 8, **kwargs):
        """
        Plot v_l and v_loc potentials on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        for l, pot in self.parser.potentials.items():

            ax.plot(pot.rmesh, pot.values,
                    label="$V_{loc}$" if l == -1 else "PS $V_{%s}$" % l2char[l],
                    **self._mpl_opts_laeps(l, "ae"))

        self.decorate_ax(ax, xlabel="r (Bohr)", ylabel="$v_l(r)$", title="Ion Pseudopotentials",
                         fontsize=fontsize,
                         )

        color = "k"
        if self.parser.lloc == 4: color = "magenta"
        ax.axvline(self.parser.rc5, lw=1, color=color, ls="--")

        self._add_rc_vlines(ax)

        return fig

    @add_fig_kwargs
    def plot_der_potentials(self, ax=None, order=1, fontsize: int = 8, **kwargs):
        """
        Plot the derivatives of vl and vloc potentials on axis ax.
        Used to analyze the derivative discontinuity introduced by the RRKJ method at rc.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        from abipy.tools.derivatives import finite_diff
        from scipy.interpolate import UnivariateSpline

        ax, fig, plt = get_ax_fig_plt(ax)

        for l, pot in self.parser.potentials.items():
            # Need linear mesh for finite_difference hence spline input potentials on lin_rmesh.
            lin_rmesh, h = np.linspace(pot.rmesh[0], pot.rmesh[-1], num=len(pot.rmesh) * 4, retstep=True)
            spline = UnivariateSpline(pot.rmesh, pot.values, s=0)
            lin_values = spline(lin_rmesh)
            vder = finite_diff(lin_values, h, order=order, acc=4)

            if l == -1:
                label = "%s-order derivative Vloc" % order
            else:
                label = "$s-order derivative PS l=%s" % str(l)

            line, = ax.plot(lin_rmesh, vder, label=label,
                             **self._mpl_opts_laeps(l, "ae"))

        self.decorate_ax(ax, xlabel="r (Bohr)", ylabel=r"$D^%s \phi(r)$" % order,
                         title="Derivative of the ion Pseudopotentials",
                         fontsize=fontsize,
                         )
        return fig

    @add_fig_kwargs
    def plot_kene_vs_ecut(self, ax=None, fontsize: int = 8, **kwargs):
        """
        Plot the convergence of the kinetic energy wrt ecut on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        for l, data in self.parser.kene_vs_ecut.items():
            ax.plot(data.energies, data.values, label="Conv l=%s" % l2char[l],
                    **self._mpl_opts_laeps(l, "ae"))

        for nlk, data in self.parser.kinerr_nlk.items():
            ax.plot(data.ecuts, data.values_ha,
                    **self._mpl_opts_laeps(nlk.l, "ps"))

        self.decorate_ax(ax, xlabel="Ecut (Ha)", ylabel=r"$\Delta E_{kin}$ (Ha)",
                         title="",
                         fontsize=fontsize,
                         )

        ax.set_yscale("log")

        return fig

    @add_fig_kwargs
    def plot_atanlogder_econv(self, ax_list=None, fontsize: int = 6, **kwargs):
        """
        Plot atan(logder) and converge of kinetic energy on the same figure.

        Return: matplotlib Figure
        """
        # Build grid of plots.
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=2, ncols=1,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        self.plot_atan_logders(ax=ax_list[0], fontsize=fontsize, show=False)
        ax_list[0].xaxis.set_label_position('top')

        self.plot_kene_vs_ecut(ax=ax_list[1], fontsize=fontsize, show=False)

        return fig

    @add_fig_kwargs
    def plot_den_formfact(self, ecut=120, ax=None, fontsize: int = 8, **kwargs):
        """
        Plot the density form factor as a function of ecut in Ha.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Return: matplotlib Figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        for name, rho in self.parser.densities.items():
            if name == "rhoC": continue
            form = rho.get_intr2j0(ecut=ecut) / (4 * np.pi)
            ax.plot(form.mesh, form.values, label=name,
                    linewidth=self.linewidth, markersize=self.markersize)

            intg = rho.r2f_integral()[-1]
            #print("r2 f integral: ", intg)
            #print("form_factor(0): ", name, form.values[0])

        # Plot vloc(q)
        #for l, pot in self.potentials.items():
        #    if l != -1: continue
        #    form = pot.get_intr2j0(ecut=ecut)
        #    mask = np.where(np.abs(form.values) > 20); form.values[mask] = 20
        #    line, = ax.plot(form.mesh, form.values, linewidth=self.linewidth, markersize=self.markersize)
        #    lines.append(line); legends.append("Vloc(q)")

        self.decorate_ax(ax, xlabel="Ecut (Ha)", ylabel="$n(q)$",
                         fontsize=fontsize,
                         title="Form factor, l=0 ",
                         )

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        Generate a predefined list of matplotlib figures with minimal input from the user.
        """
        verbose = kwargs.get("verbose", 0)

        yield self.plot_atanlogder_econv(show=False)
        yield self.plot_potentials(show=False)
        yield self.plot_radial_wfs(show=False)
        yield self.plot_radial_wfs(what="scattering_states", show=False)
        yield self.plot_projectors(show=False)
        yield self.plot_densities(show=False)
        #yield self.plot_densities(timesr2=True, show=False)
        yield self.plot_den_formfact(show=False)
        if verbose:
            yield self.plot_der_potentials(show=False)
            for order in [1, 2, 3, 4]:
                yield self.plot_der_densities(order=order, show=False)


class PseudoGenOutputParserError(Exception):
    """Exceptions raised by OutputParser."""


class PseudoGenOutputParser(metaclass=abc.ABCMeta):
    """
    Abstract class defining the interface that must be provided
    by the parsers used to extract results from the output file of
    a pseudopotential generator a.k.a. ppgen

    Attributes:

        errors: List of strings with errors reported by the pp generator
        warnings: List of strings with the warnings reported by the pp generator.
    """

    Error = PseudoGenOutputParserError

    def __init__(self, filepath: str) -> None:
        self.filepath = os.path.abspath(filepath)
        self.run_completed = False
        self._errors = []
        self._warnings = []

    @property
    def errors(self) -> List[str]:
        """
        List of strings with possible errors reported by the generator at run-time.
        """
        return self._errors

    @property
    def warnings(self) -> List[str]:
        """
        List of strings with possible errors reported by the generator at run-time.
        """
        return self._warnings

    @abc.abstractmethod
    def get_results(self):
        """
        Return the most important results in a dictionary.
        """

    @abc.abstractmethod
    def get_input_str(self) -> str:
        """Returns a string with the input file."""


# Object returned by self._grep
GrepResults = namedtuple("GrepResults", "data, start, stop")

# Used to store ae and pp quantities (e.g wavefunctions) in a single object.
AePsNamedTuple = namedtuple("AePsNamedTuple", "ae, ps")

ConvData = namedtuple("ConvData", "l energies values")

AtanLogDer = namedtuple("AtanLogDer", "l, energies, values")


class OncvOutputParser(PseudoGenOutputParser):
    """
    Object to read and extract data from the output file of oncvpsp.

    Attributes:
        atsym
        Z
        nc
        nv
        iexc
        psfile

    Example:

        parser = OncvOutputParser(filename)
        parser.scan()

        # To access data:
        parser.radial_wavefunctions

        # To plot data with matplotlib.
        p = parser.get_plotter()
        p.plot_atanlogder_econv()

    """
    # TODO Improve fully-relativistic case.

    def scan(self, verbose: int = 0) -> None:
        """
        Scan the output file, set `run_completed` attribute.

        Raises: self.Error if invalid file.
        """
        try:
            return self._scan(verbose=verbose)
        except Exception as exc:
            raise self.Error(f"Exception while parsing: {self.filepath}") from exc

    def _scan(self, verbose: int = 0) -> None:
        if not os.path.exists(self.filepath):
            raise self.Error(f"File {self.filepath} does not exist")

        # Read data and store it in lines
        self.lines = []
        import io
        with io.open(self.filepath, "rt") as fh:
            for i, line in enumerate(fh):
                #print(line)
                line = line.strip()
                self.lines.append(line)

                if verbose and line.startswith("fcfact*="):
                    print(line)

                if line.startswith("DATA FOR PLOTTING"):
                    self.run_completed = True

                # lines that contain the word ERROR but do not seem to indicate an actual teminating error
                acceptable_error_markers = [
                  'run_config: ERROR for fully non-local  PS atom,'
                ]

                if "ERROR" in line:
                    # Example:
                    # test_data: must have fcfact>0.0 for icmod= 1
                    # ERROR: test_data found   1 errors; stopping
                    if line in acceptable_error_markers:
                        self._warnings.append("\n".join(self.lines[i-1:i+1]))
                    else:
                        self._errors.append("\n".join(self.lines[i-1:i+1]))

                if "WARNING" in line:
                    self._warnings.append("\n".join(self.lines[i:i+2]))

                if "GHOST(+)" in line:
                    self._warnings.append(line)
                if "GHOST(-)" in line:
                    self._errors.append(line)

        # Get gendate, calc_type and version
        # scalar-relativistic version 2.1.1, 03/26/2014
        # scalar-relativistic version 3.0.0 10/10/2014
        toks = self.lines[1].replace(",", " ").split()
        self.gendate = toks.pop(-1)
        self.calc_type, self.version = toks[0], toks[2]

        self.major_version, self.minor_version, self.patch_level = tuple(map(int, self.version.split(".")))[:3]

        # Read configuration (not very robust because we assume the user didn't change the template but oh well)
        header = "# atsym  z    nc    nv    iexc   psfile"
        for i, line in enumerate(self.lines):
            if line.startswith("# atsym"):
                values = self.lines[i+1].split()
                keys = header[1:].split()
                # assert len(keys) == len(values)
                # Store values in self.
                for k, v in zip(keys, values):
                    # Convert nc and nv to int.
                    if k in ("nc", "nv", "iexc"): v = int(v)
                    if k in ("z", ): v = float(v)
                    setattr(self, k, v)
                break

        # Parse pseudization options for the local part.
        header = "# lloc, lpopt,  rc(5),   dvloc0"
        self.rc5 = None
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                tokens = self.lines[i + 1].split()
                #print("tokens", tokens)
                self.lloc = int(tokens[0])
                self.lptopt = int(tokens[1])
                self.rc5 = float(tokens[2])
                self.dvloc0 = float(tokens[3])
                break

        if self.rc5 is None:
            raise self.Error(f"Cannot find magic line starting with `{header}` in: {self.filepath}")

        # Parse ATOM and Reference configuration
        # Example:
        """
        #
        #   n    l    f        energy (Ha)
            1    0    2.00    -6.5631993D+01
            2    0    2.00    -5.1265474D+00
            2    1    6.00    -3.5117357D+00
            3    0    2.00    -3.9736459D-01
            3    1    2.00    -1.4998149D-01

        full rel
        in version 4, there no difference between FR and SR
        in version 3, the FR version has:

        #   n    l    f              energy (Ha)
        #   n    l    f        l+1/2             l-1/2
            1    0    2.00    -2.4703720D+03
            2    0    2.00    -4.2419865D+02

        """
        if self.relativistic and self.major_version <= 3:
            header = "#   n    l    f        l+1/2             l-1/2"
        else:
            header = "#   n    l    f        energy (Ha)"

        nc, nv = self.nc, self.nv

        for i, line in enumerate(self.lines):
            if line.startswith(header):
                beg, core = i + 1, [],
                for c in range(nc):
                    n, l, f = self.lines[beg+c].split()[:3]
                    if is_integer(f):
                        f = str(int(float(f)))
                    else:
                        f = "%.1f" % f
                    core.append(n + l2char[l] + "^%s" % f)
                #self.core = "$" + " ".join(core) + "$"

                beg, valence = i + nc + 1, []
                for v in range(nv):
                    #print("lines[beg+v]", self.lines[beg+v])
                    n, l, f = self.lines[beg+v].split()[:3]
                    if is_integer(f):
                        f = str(int(float(f)))
                    else:
                        #print(f)
                        f = "%.1f" % float(f)

                    valence.append(n + l2char[l] + "^{%s}" % f)

                #self.valence = "$" + " ".join(valence) + "$"
                #print("core", self.core, "valence",self.valence)
                break
        else:
            raise self.Error(f"Cannot find header:\n`{header}`\nin output file {self.filepath}")

    @lazy_property
    def lmax(self) -> int:
        # Read lmax (not very robust because we assume the user didn't change the template but oh well)
        header = "# lmax"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                return int(self.lines[i+1])
                break
        else:
            raise self.Error(f"Cannot find line with `#lmax` in: {self.filepath}")

    def to_string(self, verbose: int = 0 ) -> str:
        """
        String representation.
        """
        lines = []
        app = lines.append

        if not hasattr(self, "calc_type"):
            app("Object is empty. Call scan method to analyze output file")
            return "\n".join(lines)

        if not self.run_completed:
            app("completed: %s" % self.run_completed)
            return "\n".join(lines)

        app("%s, oncvpsp version: %s, date: %s" % (self.calc_type, self.version, self.gendate))

        from pprint import pformat
        app(pformat(self.get_results()))

        if self.warnings:
            lines.extend(self.warnings)

        if self.errors:
            lines.extend(self.errors)

        return "\n".join(lines)

    def __str__(self) -> str:
        return self.to_string()

    @property
    def relativistic(self) -> bool:
        """True if fully-relativistic calculation."""
        return self.calc_type in ("fully-relativistic", "relativistic")

    @lazy_property
    def rc_l(self) -> dict[int, float]:
        """
        Core radii as a function of l extracted from the output file.
        """
        rc_l = {}
        header = "#   l,   rc,"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                beg = i + 1
                nxt = 0
                while True:
                    ln = self.lines[beg + nxt]
                    if ln.startswith("#"): break
                    tokens = ln.split()
                    #print("line:", ln, "\ntokens", tokens)
                    l, rc = int(tokens[0]), float(tokens[1])
                    rc_l[l] = rc
                    nxt += 1

        if not rc_l:
            raise self.Error(f"Cannot find magic line starting with `{header}` in: {self.filepath}")

        return rc_l

    @lazy_property
    def kinerr_nlk(self) -> dict[NlkState, namedtuple]:
        """
        Dictionary with the error on the kinetic energy indexed by nlk.
        """

        # In relativist mode we write data inside the following loops:

        #do l1=1,lmax+1
        #   ll=l1-1
        #   if(ll==0) then
        #    mkap=1
        #   else
        #    mkap=2
        #   end if
        #   do ikap=1,mkap
        #       if(ikap==1) then
        #         kap=-(ll+1)
        #       else
        #         kap=  ll
        #       end if

        kinerr_nlk = {}

        if self.major_version > 3:
            # Calculating optimized projector #   1
            #
            #  for l=   0

            re_start = re.compile(r"^Calculating optimized projector #\s+(?P<iproj>\d+)")

        else:
            # Calculating first optimized projector for l=   0
            re_start = re.compile(r"^Calculating (?P<iproj>(first|second)) optimized projector for l=\s+(?P<l>\d+)")
            # TODO: In FR mode, we have
            #Calculating first optimized projector for l=   0
            #Calculating second optimized projector for l=   0

        nlk = None
        iproj_l_seen = set()

        for i, line in enumerate(self.lines):

            m = re_start.match(line)

            if m:
                # Extract iproj and l.
                if self.major_version > 3:
                    # for l=   0
                    iproj = int(m.group("iproj"))
                    l = int(self.lines[i+2].split("=")[-1].strip())
                else:
                    iproj = m.group("iproj")
                    iproj = {"first": 0, "second": 1}[iproj]
                    l = int(m.group("l"))

                k = None
                if self.relativistic:
                    k = 1
                    if (iproj, l) in iproj_l_seen: k= 2
                    iproj_l_seen.add((iproj, l))

                # Use n index to store iprj index.
                nlk = NlkState(n=iproj, l=l, k=k)
                #print("nlk:", nlk)
                continue

            # Now parse the following section associated to nlk

            #Energy error per electron        Cutoff
            #     Ha          eV             Ha
            #     0.01000     0.27211       27.01
            #     0.00100     0.02721       52.82
            #     0.00010     0.00272       66.22
            #     0.00001     0.00027       75.37

            if line.startswith("Energy error per electron        Cutoff"):
                values_ha, ecuts = [], []
                for j in range(4):
                    tokens = self.lines[i+2+j].split()
                    #print("tokens:", tokens)
                    if not tokens: break
                    err_ha, err_ev, ecut = map(float, tokens)
                    values_ha.append(err_ha)
                    ecuts.append(ecut)

                if nlk is None:
                    raise self.Error("Cannot find nlk quantum numbers")

                self._check_nlk_key(nlk, kinerr_nlk, "kinerr_nlk")

                kinerr_nlk[nlk] = dict2namedtuple(ecuts=ecuts, values_ha=values_ha)

        if not kinerr_nlk:
            raise self.Error(f"Cannot parse convergence profile in: {self.filepath}")

        return kinerr_nlk

    @staticmethod
    def _check_nlk_key(nlk, d, dict_name):

        if nlk in d:
            ks = "\n\t".join(str(k) for k in d)
            raise RuntimeError(f"nlk state `{nlk}` is already in {dict_name}:\nKeys:\n\t{ks}")

    @lazy_property
    def potentials(self) -> dict[int, RadialFunction]:
        """
        Dict with radial functions with the non-local and local potentials indexed by l.
        l = -1 corresponds to the local part (if present).
        """
        #radii, charge, pseudopotentials (ll=0, 1, lmax)
        #!p   0.0099448   4.7237412  -7.4449470 -14.6551019
        vl_data = self._grep("!p").data
        lmax = len(vl_data[0]) - 3
        assert lmax == self.lmax

        # From 0 up to lmax
        ionpots_l = {}
        for l in range(lmax + 1):
            ionpots_l[l] = RadialFunction("Ion Pseudopotential, l=%d" % l, vl_data[:, 0], vl_data[:, 2+l])

        # Local part is stored with l == -1 if lloc=4, not present if lloc = l
        vloc = self._grep("!L").data
        if vloc is not None:
            ionpots_l[-1] = RadialFunction("Local part, l=%d" % -1, vloc[:, 0], vloc[:, 1])

        return ionpots_l

    @lazy_property
    def densities(self) -> dict[str, RadialFunction]:
        """
        Dictionary with charge densities on the radial mesh.
        """
        # radii, charge, core charge, model core charge
        # !r   0.0100642   4.7238866  53.4149287   0.0000000
        rho_data = self._grep("!r").data

        return dict(
            rhoV=RadialFunction("Valence charge", rho_data[:, 0], rho_data[:, 1]),
            rhoC=RadialFunction("Core charge", rho_data[:, 0], rho_data[:, 2]),
            rhoM=RadialFunction("Model charge", rho_data[:, 0], rho_data[:, 3])
        )

    @lazy_property
    def radial_wfs(self) -> AePsNamedTuple:
        """
        Read and set the radial wavefunctions for the bound states.

        Usage:

            ae_wfs, ps_wfs = self.radial_wfs.ae, self.radial_wfs.ps

            for nlk, ae_wf in ae_wfs.items():
                ps_wf, l, k = ps_wfs[nlk], nlk.l, nlk.k
        """
        return self._get_radial_wavefunctions(what="bound_states")

    @lazy_property
    def scattering_wfs(self) -> AePsNamedTuple:
        """
        Read and set the scattering wavefunctions.
        """
        return self._get_radial_wavefunctions(what="scattering_states")

    def _get_radial_wavefunctions(self, what: str) -> AePsNamedTuple:
        # For scalar-relativistic bound states, we have
        #
        #    n= 1,  l= 0, all-electron wave function, pseudo w-f
        #
        #    &     0    0.009945   -0.092997    0.015273

        # bound states, fully-relativistic case:
        #
        #    n= 1,  l= 0  kap=-1, all-electron wave function, pseudo w-f
        #
        #    &     0    0.009955    0.066338    0.000979

        # For the scattering states (scalar and relativistic case)
        #
        # scattering, iprj= 2,  l= 1, all-electron wave function, pseudo w-f
        #
        # scattering, iprj= 2,  l= 1, kap= 1, all-electron wave function, pseudo w-f

        ae_waves, ps_waves = {}, {}

        beg = 0
        while True:
            g = self._grep("&", beg=beg)
            if g.data is None: break
            beg = g.stop + 1

            # Get header two lines above.
            header = self.lines[g.start - 2]

            if what == "bound_states":
                if header.startswith("scattering,"):
                    continue
            elif what == "scattering_states":
                if not header.startswith("scattering,"):
                    continue
                header = header.replace("scattering,", "")
            else:
                raise ValueError(f"Invalid value of what: `{what}`")

            #print(header)

            if not self.relativistic:
                # n= 1,  l= 0, all-electron wave function, pseudo w-f
                n, l = header.split(",")[0:2]
                n = int(n.split("=")[1])
                l = int(l.split("=")[1])
                kap = None
            else:
                # n= 1,  l= 0,  kap=-1, all-electron wave function, pseudo w-f
                if self.major_version <= 2: header = header.replace("kap=", ", kap=")
                n, l, kap = header.split(",")[0:3]
                n = int(n.split("=")[1])
                l = int(l.split("=")[1])
                kap = int(kap.split("=")[1])

            nlk = NlkState.from_nlkap(n=n, l=l, kap=kap)
            #print("Got nlk state:", nlk)

            rmesh = g.data[:, 1]
            ae_wf = g.data[:, 2]
            ps_wf = g.data[:, 3]

            self._check_nlk_key(nlk, ae_waves, "ae_waves")

            ae_waves[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, ae_wf)
            ps_waves[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, ps_wf)

        return AePsNamedTuple(ae=ae_waves, ps=ps_waves)

    @lazy_property
    def projectors(self) -> dict[NlkState, RadialFunction]:
        """
        Dict with projector wave functions indexed by nlk.
        """
        #
        #@     0    0.009945    0.015274   -0.009284
        beg = 0
        magic = "@"
        if self.major_version > 3: magic = "!J"

        # if(ikap==1) then
        #   write(6,'(a,i6,6(f12.6,1x))') '!J',-ll,rr(ii), &
        #        (vkb(ii,jj,l1,ikap),jj=1,nproj(l1))
        # else
        #   write(6,'(a,i6,6(f12.6,1x))') '!J',ll,rr(ii), &
        #        (vkb(ii,jj,l1,ikap),jj=1,nproj(l1))

        projectors_nlk = {}
        while True:
            g = self._grep(magic, beg=beg)
            if g.data is None: break
            beg = g.stop + 1

            rmesh = g.data[:, 1]
            l = int(g.data[0, 0])

            k = None
            if self.relativistic:
                k = 2
                if l <= 0: k = 1

            for n in range(len(g.data[0]) - 2):
                nlk = NlkState(n=n + 1, l=abs(l), k=k)
                #print("Got projector with: %s" % str(nlk))

                if nlk in projectors_nlk:
                    raise self.Error("nlk state `{nlk}` is already in projectors_nlk")

                projectors_nlk[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, g.data[:, n + 2])

        return projectors_nlk

    @lazy_property
    def atan_logders(self) -> AePsNamedTuple:
        """
        Atan of the log derivatives for different l-values.
        """
        #log derivativve data for plotting, l= 0
        #atan(r * ((d psi(r)/dr)/psi(r))), r=  1.60
        #l, energy, all-electron, pseudopotential
        #
        #!      0    2.000000    0.706765    0.703758
        ae_atan_logder_l, ps_atan_logder_l = {}, {}

        lstop = self.lmax + 1
        if self.major_version > 3:
            lstop = min(self.lmax + 2, 4)

        if not self.relativistic:
            l_list = list(range(lstop))
        else:
            # Order with (l, -l) for plotting purposes.
            l_list = []
            for l in range(lstop):
                if l != 0:
                    l_list.extend((l, -l))
                else:
                    l_list.append(0)

        for l in l_list:
            tag = "!      %d" % l if l >= 0 else "!     %d" % l
            data = self._grep(tag=tag).data
            if data is None:
                raise self.Error(f"Cannot find logder for l: {l}")
            assert l == int(data[0, 0])

            ae_atan_logder_l[l] = AtanLogDer(l=l, energies=data[:, 1], values=data[:, 2])
            ps_atan_logder_l[l] = AtanLogDer(l=l, energies=data[:, 1], values=data[:, 3])

        return AePsNamedTuple(ae=ae_atan_logder_l, ps=ps_atan_logder_l)

    @lazy_property
    def kene_vs_ecut(self) -> dict[int, ConvData]:
        """
        Dict with the convergence of the kinetic energy versus ecut for different l-values.
        """
        #convergence profiles, (ll=0,lmax)
        #!C     0    5.019345    0.010000
        #...
        #!C     1   19.469226    0.010000
        # TODO: This does not take into account scattering states or n > 1
        conv_l = {}

        for l in range(self.lmax + 1):
            data = self._grep(tag="!C     %d" % l).data
            conv_l[l] = ConvData(l=l, energies=data[:, 1], values=data[:, 2])

        return conv_l

    @lazy_property
    def hints(self) -> dict:
        """
        Hints for the cutoff energy as provided by oncvpsp.
        """
        # Extract the hints
        hints = 3 * [-np.inf]
        for i in range(3):
            for l in range(self.lmax + 1):
                hints[i] = max(hints[i], self.kene_vs_ecut[l].energies[-i-1])
        hints.reverse()

        # Truncate to the nearest int
        hints = [np.rint(h) for h in hints]
        # print("hints:", hints)

        return dict(
            low={"ecut": hints[0], "pawecutdg": hints[0]},
            normal={"ecut": hints[1], "pawecutdg": hints[1]},
            high={"ecut": hints[2], "pawecutdg": hints[2]}
        )

    def get_results(self) -> AttrDict:
        """
        Return the most important results extracted from the output file.
        """
        # Init return values
        #d = AttrDict(
        #    max_ecut=None,
        #    max_atan_logder_l1err=None,
        #    max_psexc_abserr=None,
        #    herm_err=None,
        #    nwarns=len(self.warnings)
        #    nerrs=len(self.errors)
        #)

        # Get the max ecut estimated by oncvpsp.
        # TODO: Should take into account scattering states.
        max_ecut = max(self.kene_vs_ecut[l].energies[-1] for l in self.kene_vs_ecut)

        # Compute the l1 error in atag(logder) between AE and PS
        from scipy.integrate import cumtrapz
        max_l1err = 0.0
        for l in self.atan_logders.ae:
            f1, f2 = self.atan_logders.ae[l], self.atan_logders.ps[l]

            abs_diff = np.abs(f1.values - f2.values)
            integ = cumtrapz(abs_diff, x=f1.energies) / (f1.energies[-1] - f1.energies[0])
            max_l1err = max(max_l1err, integ[-1])

        # Read Hermiticity error and compute the max value of PSP excitation error=
        # Hermiticity error    4.8392D-05
        # PSP excitation error=  1.56D-10
        herm_tag, pspexc_tag = "Hermiticity error", "PSP excitation error="
        herm_err, max_psexc_abserr = None, -np.inf

        for line in self.lines:
            i = line.find(herm_tag)
            if i != -1:
                herm_err = float(line.split()[-1].replace("D", "E"))

            i = line.find(pspexc_tag)
            if i != -1:
                max_psexc_abserr = max(max_psexc_abserr, abs(float(line.split()[-1].replace("D", "E"))))

        return AttrDict(
            max_ecut=max_ecut,
            max_atan_logder_l1err=max_l1err,
            max_psexc_abserr=max_psexc_abserr,
            herm_err=herm_err,
            nwarns=len(self.warnings),
            nerrs=len(self.errors),
        )

    def find_string(self, s: str) -> int:
        """
        Returns the index of the first line containing string s.
        Raises self.Error if s cannot be found.
        """
        for i, line in enumerate(self.lines):
            if s in line:
                return i
        else:
            raise self.Error(f"Cannot find `{s}` in lines")

    def get_input_str(self) -> str:
        """String with the ONCVPSP input file."""
        try:
            # oncvpsp 3.2.3
            i = self.find_string("<INPUT>")
            j = self.find_string("</INPUT>")
            return "\n".join(self.lines[i+1:j]) + "\n"
        except self.Error:
            # oncvpsp => 4
            i = self.find_string("Reference configufation results")
            return "\n".join(self.lines[:i])

    def get_psp8_str(self) -> Union[str, None]:
        """
        Return string with the pseudopotential data in psp8 format.
        Return None if field is not present.
        """
        start, stop = None, None
        for i, line in enumerate(self.lines):
            if 'Begin PSPCODE8' in line: start = i
            if start is not None and 'END_PSP' in line:
                stop = i
                break

        if start is None and stop is None: return None
        ps_data = "\n".join(self.lines[start+1:stop])

        if "<INPUT>" not in ps_data:
            # oncvpsp <= 3.2.2 --> Append the input to ps_data (note XML markers)
            # oncvpsp >= 3.2.3 --> Input is already there
            ps_data += "\n\n<INPUT>\n" + self.get_input_str() + "</INPUT>\n"

        return ps_data

    def get_upf_str(self) -> Union[str, None]:
        """
        Return string with the pseudopotential data in upf format.
        Return None if field is not present.
        """
        start, stop = None, None
        for i, line in enumerate(self.lines):
            if "Begin PSP_UPF" in line: start = i
            if start is not None and 'END_PSP' in line:
                stop = i
                break

        if start is None and stop is None: return None
        return "\n".join(self.lines[start+1:stop])

    def get_plotter(self) -> Union[OncvPlotter, None]:
        """
        Return an instance of OncvPlotter or None
        """
        try:
            return OncvPlotter(self)
        except Exception as exc:
            #raise
            return None

    def _grep(self, tag: str, beg: int = 0) -> GrepResults:
        """
        Finds the first field in the file with the specified tag.
        `beg` gives the initial position in the file.
        """
        data, stop, intag = [], None, -1

        if beg >= len(self.lines):
            raise ValueError(f"beg {beg} > len(lines) ({len(self.lines)})")

        for i, l in enumerate(self.lines[beg:]):
            l = l.lstrip()
            if l.startswith(tag):
                if intag == -1:
                    intag = beg + i
                data.append([float(c) for c in l.split()[1:]])
            else:
                # Exit because we know there's only one section starting with 'tag'
                if intag != -1:
                    stop = beg + i
                    break

        if not data:
            return GrepResults(data=None, start=intag, stop=stop)
        else:
            return GrepResults(data=np.array(data), start=intag, stop=stop)

    def gnuplot(self) -> None:
        """
        Plot the results with gnuplot.
        Based on the `replot.sh` script provided by the oncvpsp code.
        """
        outfile = self.filepath
        base = os.path.basename(outfile)
        gnufile = base + ".scr"
        plotfile = base + ".plot"
        temp = base + ".tmp"

        from monty.os import cd
        from subprocess import check_call
        workdir = tempfile.mkdtemp()
        print(f"Working in: {workdir}")

        with cd(workdir):
            check_call("awk 'BEGIN{out=0};/GNUSCRIPT/{out=0}; {if(out == 1) {print}}; \
                                /DATA FOR PLOTTING/{out=1}' %s > %s" % (outfile, plotfile), shell=True)

            check_call("awk 'BEGIN{out=0};/END_GNU/{out=0}; {if(out == 1) {print}}; \
                                /GNUSCRIPT/{out=1}' %s > %s" % (outfile, temp), shell=True)

            check_call('sed -e 1,1000s/t1/"%s"/ %s > %s' % (plotfile, temp, gnufile), shell=True)

            try:
                check_call(["gnuplot", gnufile])
            except KeyboardInterrupt:
                print("Received KeyboardInterrupt")

        os.rmdir(workdir)


def oncv_make_open_notebook(outpath: str,
                            foreground: bool = False,
                            classic_notebook: bool = False,
                            no_browser: bool = False) -> int:  # pragma: no cover
    """
    Generate an ipython notebook and open it in the browser.

    Args:
        foreground: By default, jupyter is executed in background and stdout, stderr are redirected
            to devnull. Use foreground to run the process in foreground
        classic_notebook: True to use the classic notebook instead of jupyter-lab (default)
        no_browser: Start the jupyter server to serve the notebook but don't open the notebook in the browser.
            Use this option to connect remotely from localhost to the machine running the kernel

    Return: system exit code.

    Raise: `RuntimeError` if jupyter executable is not in $PATH
    """
    nbpath = oncv_write_notebook(outpath, nbpath=None)

    if not classic_notebook:
        # Use jupyter-lab.
        app_path = which("jupyter-lab")
        if app_path is None:
            raise RuntimeError("""
Cannot find jupyter-lab application in $PATH. Install it with:

    conda install -c conda-forge jupyterlab

or:

    pip install jupyterlab

See also https://jupyterlab.readthedocs.io/
""")

    else:
        # Use classic notebook
        app_path = which("jupyter")
        if app_path is None:
            raise RuntimeError("""
Cannot find jupyter application in $PATH. Install it with:

    conda install -c conda-forge jupyter

or:

    pip install jupyterlab

See also https://jupyter.readthedocs.io/en/latest/install.html
""")
        app_path = app_path + " notebook "

    if not no_browser:
        if foreground:
            return os.system("%s %s" % (app_path, nbpath))
        else:
            fd, tmpname = tempfile.mkstemp(text=True)
            print(tmpname)
            cmd = "%s %s" % (app_path, nbpath)
            print("Executing:", cmd, "\nstdout and stderr redirected to %s" % tmpname)
            import subprocess
            process = subprocess.Popen(cmd.split(), shell=False, stdout=fd, stderr=fd)
            cprint("pid: %s" % str(process.pid), "yellow")
            return 0

    else:
        # Based on https://github.com/arose/nglview/blob/master/nglview/scripts/nglview.py
        notebook_name = os.path.basename(nbpath)
        dirname = os.path.dirname(nbpath)
        print("nbpath:", nbpath)

        import socket

        def find_free_port():
            """https://stackoverflow.com/questions/1365265/on-localhost-how-do-i-pick-a-free-port-number"""
            from contextlib import closing
            with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
                s.bind(('', 0))
                s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                return s.getsockname()[1]

        username = os.getlogin()
        hostname = socket.gethostname()
        port = find_free_port()

        client_cmd = "ssh -NL localhost:{port}:localhost:{port} {username}@{hostname}".format(
            username=username, hostname=hostname, port=port)

        print(f"""
Using port: {port}

\033[32m In your local machine, run: \033[0m

                {client_cmd}

\033[32m NOTE: you might want to replace {hostname} by full hostname with domain name \033[0m
\033[32m Then open your web browser, copy and paste the URL: \033[0m

http://localhost:{port}/notebooks/{notebook_name}
""")
        if not classic_notebook:
            cmd = f'{app_path} {notebook_name} --no-browser --port {port} --notebook-dir {dirname}'
        else:
            cmd = f'{app_path} notebook {notebook_name} --no-browser --port {port} --notebook-dir {dirname}'

        print("Executing:", cmd)
        print('NOTE: make sure to open `{}` in your local machine\n'.format(notebook_name))

        return os.system(cmd)


def oncv_write_notebook(outpath: str, nbpath: Optional[str] = None) -> str:
    """
    Write an ipython notebook to nbpath
    If nbpath is None, a temporay file is created.
    Return path to the notebook.
    """
    outpath = os.path.abspath(outpath)

    import nbformat
    nbf = nbformat.v4
    nb = nbf.new_notebook()

    nb.cells.extend([
        nbf.new_markdown_cell("## This is an auto-generated notebook for %s" % os.path.basename(outpath)),
        nbf.new_code_cell("""\
%matplotlib notebook

# Use this magic for jupyterlab.
# For installation instructions, see https://github.com/matplotlib/jupyter-matplotlib
#%matplotlib widget

"""),

        nbf.new_code_cell("""\
# Parse output file
from abipy.ppcodes.oncvpsp import OncvOutputParser
onc_parser = OncvOutputParser('%s')""" % outpath),

        nbf.new_code_cell("""\
# Parse the file and build the plotter
onc_parser.scan()
if not onc_parser.run_completed:
    raise RuntimeError("Cannot parse output file")

plotter = onc_parser.get_plotter()"""),

        nbf.new_markdown_cell(r"# AE and PS radial wavefunctions $\phi(r)$:"),
        nbf.new_code_cell("fig = plotter.plot_radial_wfs(show=False)"),

        nbf.new_markdown_cell("# Arctan of the logarithmic derivatives:"),
        nbf.new_code_cell("fig = plotter.plot_atan_logders(show=False)"),

        nbf.new_markdown_cell("# Convergence in $G$-space estimated by ONCVPSP:"),
        nbf.new_code_cell("fig = plotter.plot_kene_vs_ecut(show=False)"),

        nbf.new_markdown_cell("# Projectors:"),
        nbf.new_code_cell("fig = plotter.plot_projectors(show=False)"),

        nbf.new_markdown_cell("# Core-Valence-Model charge densities:"),
        nbf.new_code_cell("fig = plotter.plot_densities(show=False)"),

        nbf.new_markdown_cell("# Local potential and $l$-dependent potentials:"),
        nbf.new_code_cell("fig = plotter.plot_potentials(show=False)"),

        #nbf.new_markdown_cell("# 1-st order derivative of $v_l$ and $v_{loc}$ computed via finite differences:"),
        #nbf.new_code_cell("""fig = plotter.plot_der_potentials(order=1, show=False)"""),
        #nbf.new_markdown_cell("# 2-nd order derivative of $v_l$ and $v_{loc}$ computed via finite differences:"),
        #nbf.new_code_cell("""fig = plotter.plot_der_potentials(order=2, show=False)"""),

        #nbf.new_markdown_cell("Model core charge and form factors computed by ABINIT"),
        #nbf.new_code_cell("""\
#with pseudo.open_pspsfile() as psps:
#psps.plot()"""),
    ])

    # Plot data
    #plotter.plot_der_potentials()
    #for order in [1,2,3,4]:
    #    plotter.plot_der_densities(order=order)
    #plotter.plot_densities(timesr2=True)
    #plotter.plot_den_formfact()

    if nbpath is None:
        _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)

    with io.open(nbpath, 'wt', encoding="utf8") as f:
        nbformat.write(nb, f)

    return nbpath


class MultiOncvPlotter:
    """
    Class for comparing multiple pseudos.

    Usage example:

    .. code-block:: python

        plotter = MultiOncvPlotter.from_files(filepaths)
    """

    @classmethod
    def from_files(cls, files: List[str]) -> MultiOncvPlotter:
        """
        Create an instance from a list of oncvpsp output files.
        """
        new = cls()
        for file in files:
            new.add_file(file, file)

        return new

    def __init__(self):
        self._plotters_dict = {}

    def add_file(self, label: str, filepath: str) -> None:
        """
        Add a oncvps output file to the plotter with label
        """
        if label in self._plotters_dict:
            raise ValueError(f"Cannot overwrite label: {label}")

        parser = OncvOutputParser(filepath)
        parser.scan()
        plotter = parser.get_plotter()
        if plotter is not None:
            self._plotters_dict[label] = plotter

    def __len__(self) -> int:
        return len(self._plotters_dict)

    @property
    def plotters(self):
        """"List of registered `Plotters`."""
        return list(self._plotters_dict.values())

    @property
    def labels(self) -> List[str]:
        """List of labels."""
        return list(self._plotters_dict.keys())

    def items(self):
        return self._plotters_dict.items()

    def _get_ax_list(self, ax_list, ravel=True):

        num_plots, ncols, nrows = len(self), 1, len(self)

        # Build grid of plots.
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        if ravel:
            ax_list = ax_list.ravel()

        return ax_list, fig, plt

    @add_fig_kwargs
    def plot_atan_logders(self, ax_list=None, with_xlabel=True, fontsize: int = 8, **kwargs):
        """
        Plot arctan of logder on axis ax.

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_atan_logders(ax=ax, with_xlabel=with_xlabel, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_radial_wfs(self, ax_list=None, what="bound_states", fontsize: int = 8, **kwargs):
        """
        Plot AE and PS radial wavefunctions of ax_list.

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            what: "bound_states" or "scattering_states".
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_radial_wfs(ax=ax, what=what, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_projectors(self, ax_list=None, fontsize: int = 8, **kwargs):
        """
        Plot projectors on ax_list.

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_projectors(ax=ax, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_densities(self, ax_list=None, timesr2=False, fontsize: int = 8, **kwargs):
        """
        Plot AE, PS and model densities on ax_list

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_densities(ax=ax, timesr2=timesr2, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_der_densities(self, ax_list=None, order=1, fontsize=8, **kwargs):
        """
        Plot the radial derivatives of the densities on ax_list
        Used to analyze possible discontinuities or strong oscillations in r-space.

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_der_densities(ax=ax, order=order, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_potentials(self, ax_list=None, fontsize: int = 8, **kwargs):
        """
        Plot v_l and v_loc potentials on ax_list

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_potentials(ax=ax, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_der_potentials(self, ax_list=None, order=1, fontsize: int = 8, **kwargs):
        """
        Plot the derivatives of vl and vloc potentials on ax_list.
        Used to analyze the derivative discontinuity introduced by the RRKJ method at rc.

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_der_potentials(ax=ax, order=order, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_kene_vs_ecut(self, ax_list=None, fontsize: int = 8, **kwargs):
        """
        Plot the convergence of the kinetic energy wrt ecut on ax_list

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_kene_vs_ecut(ax=ax, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_atanlogder_econv(self, ax_mat=None, fontsize: int = 6, **kwargs):
        """
        Plot atan(logder) and converge of kinetic energy on the same figure.
        Return: matplotlib Figure
        """
        num_plots, ncols, nrows = 2 * len(self), 2, len(self)

        # Build grid of plots.
        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        for i, (label, plotter) in enumerate(self.items()):
            ax_list = ax_mat[i]
            plotter.plot_atanlogder_econv(ax_list=ax_list, fontsize=fontsize, show=False)
            for ax in ax_list:
                ax.set_title(label, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_den_formfact(self, ecut=120, ax_list=None, fontsize: int = 8, **kwargs):
        """
        Plot the density form factor as a function of ecut in Ha on ax_list

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.

        Return: matplotlib Figure.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list)

        for ax, (label, plotter) in zip(ax_list, self.items()):
            plotter.plot_den_formfact(ax=ax, ecut=ecut, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        Generate a predefined list of matplotlib figures with minimal input from the user.
        """
        verbose = kwargs.get("verbose", 0)

        #yield self.plot_atanlogder_econv(show=False)
        yield self.plot_atan_logders(show=False)
        yield self.plot_kene_vs_ecut()
        yield self.plot_radial_wfs(show=False)
        yield self.plot_radial_wfs(what="scattering_states", show=False)
        yield self.plot_projectors(show=False)
        yield self.plot_potentials(show=False)
        yield self.plot_densities(show=False)
        ##yield self.plot_densities(timesr2=True, show=False)
        yield self.plot_den_formfact(show=False)
        #if verbose:
        #    yield self.plot_der_potentials(show=False)
        #    for order in [1, 2, 3, 4]:
        #        yield self.plot_der_densities(order=order, show=False)


def psp8_get_densities(path, fc_file=None, ae_file=None, plot=False):
    """
    Extract AE-core, AE-valence and PS-valence from a psp8 file.
    Optionally write `.fc` and `.AE` file in format suitable for AIM and cut3d Hirshfeld code.

    Args:
        path: path of the psp8 file
        fc_file, ae_file: File-like object to `.fc.` and `.AE` file
            Set to None if files are now wanted.
        plot: If true, call matplotlib to plot densities.

    Return:
        namedtuple with numpy arrays: (rmesh, psval, aeval, aecore)
            Densities are multiplied by 4 pi i.e. int aecore * r**2 dr ~= N_core

    .. warning::

        The densities do not integrate exactly to the number of core/valence electrons.
        in particular the valence charge in elements with large Z as valence electrons are
        more delocalized and the finite radial mesh does not capture the contribution due to the tail.
        Client code is responsible for extrapolating the valence charge if needed.
        The error in the integral of the core charge is usually much smaller and should be ok
        for qualitative analysis.
    """
    # From 64_psp/psp8in.F90
    # (1) title (character) line
    # (2) znucl,zion,pspdat
    # (3) pspcod,pspxc,lmax,lloc,mmax,r2well  (r2well not used)
    # (4) rchrg,fchrg,qchrg  (fchrg /=0 if core charge, qchrg not used)
    # (5) nproj(0:lmax)  (several projectors allowed for each l)
    # Then, for ll=0,lmax:
    # if nproj(ll)>0
    #   1/<u1|vbkb1>, 1/<u2|vbkb2>, ...
    #   for irad=1,mmax: irad, r(irad), vbkb1(irad,ll), vbkb2(irad,ll), ...
    # elif ll=lloc
    #   for  irad=1,mmax: irad, r(irad), vloc(irad)
    # end if
    #
    # If (lloc>lmax)
    # for  irad=1,mmax: irad, r(irad), vloc(irad)
    # end if
    #
    # vbkb are Bloechl-Kleinman-Bylander projectors,(vpsp(r,ll)-vloc(r))*u(r,ll),
    # unnormalized
    # Note that an arbitrary local potential is allowed.  Set lloc>lmax, and
    # provide projectors for all ll<=lmax
    #
    # Finally, if fchrg>0,
    # for  irad=1,mmax: irad, r(irad), xccc(irad),
    # xccc'(irac), xccc''(irad), xccc'''(irad), xccc''''(irad)
    #
    # Model core charge for nonlinear core xc correction, and 4 derivatives
    from abipy.flowtk.pseudos import Pseudo
    pseudo = Pseudo.from_file(path)

    from pymatgen.io.abinit.pseudos import _dict_from_lines
    with open(path, "rt") as fh:
        lines = [fh.readline() for _ in range(6)]
        header = _dict_from_lines(lines[1:3], [3, 6])

        # Number of points on the linear mesh.
        mmax = int(header["mmax"])

        # Non-linear core correction parameters.
        rchrg, fchrg, qchrg = [float(t) for t in lines[3].split()[:3]]

        # Read Number of projectors(l) and extension switch
        nproj = [int(t) for t in lines[4].split()[:5]]
        assert len(nproj) == 5
        tokens = lines[5].split()
        extension_switch = int(tokens[0])

        # Old format, Densities are not available.
        if len(tokens) == 1 or tokens[1] == "extension_switch":
            raise RuntimeError("psp8 file does not contain density records")

        den_flag = int(tokens[1])
        if den_flag != 1:
            raise ValueError("Expecting den_flag 1 but got %s" % den_flag)

        # Read SOC projectors
        has_soc = extension_switch in (2, 3)
        if has_soc:
            line = fh.readline()
            # Start at l=1
            nproj_soc = [int(t) for t in line.split()[:4]]
            nproj_soc.insert(0, 0)
            #print("nproj_soc", nproj_soc)
            raise NotImplementedError("SOC not tested")

        lmax = int(header["lmax"])
        lloc = int(header["lloc"])
        nso = 1 if not has_soc else 2

        # Will now proceed at the reading of pots and projectors
        # rad(:)=radial grid r(i)
        # vpspll(:,1),...,vpspll(:,lnmax)=nonlocal projectors
        # vloc(:)=local potential

        #for nn in range(nso):
        # Skip projectors (scalar relativistic, always present).
        for l, npl in enumerate(nproj):
            #if npl == 0 and l != lloc: continue
            if npl == 0: continue
            line = fh.readline() # l, ekb[:npl]
            l_file = int(line.split()[0])
            if l != l_file:
                #print("For l=%s, npl=%s" % (l, npl), "wrong line", line)
                raise RuntimeError("l != l_file (%s != %s)" % (l, l_file))

            for ir in range(mmax):
                line = fh.readline()
                assert int(line.split()[0]) == ir + 1

        # Skip local potential.
        if lloc == 4:
            lloc_file = int(fh.readline())
            assert lloc_file == lloc
            for ir in range(mmax):
                fh.readline()

        # Skip model core charge function and derivatives, if present.
        if fchrg > 1e-15:
            for ir in range(mmax):
                fh.readline()

        # Read pseudo valence charge in real space on the linear mesh.
        # [i, r, PS_val, AE_val, AE_core]
        rmesh, psval, aeval, aecore = [np.empty(mmax) for _ in range(4)]
        for ir in range(mmax):
            l = fh.readline()
            #print("denline", l)
            findx, rad, v1, v2, v3 = l.split()
            assert ir + 1 == int(findx)
            rmesh[ir] = float(rad.replace("D", "E"))
            psval[ir] = float(v1.replace("D", "E"))
            aeval[ir] = float(v2.replace("D", "E"))
            aecore[ir] = float(v3.replace("D", "E"))

        #fact = 1 / (4 * np.pi)
        #aeval *= fact
        #psval *= fact
        #aecore *= fact
        from scipy.integrate import simps
        r2 = rmesh ** 2

        meta = dict(
            aeval_integral=simps(aeval * r2, x=rmesh),
            psval_integral=simps(psval * r2, x=rmesh),
            aecore_integral=simps(aecore * r2, x=rmesh),
            aeden_integral=simps((aecore + aeval) * r2, x=rmesh),
            symbol=pseudo.symbol,
            Z=pseudo.Z,
            Z_val=pseudo.Z_val,
            l_max=pseudo.l_max,
            md5=pseudo.md5,
        )

        import json
        if fc_file is not None:
            # Write results to file in "fc" format. The files contains:
            # header with number of points and unknown parameter (not used)
            # for each point in the radial mesh:
            #    4 columns with the radial r coordinate, the core density at r,
            #    and the first and second derivatives of the core density.
            # See http://www.abinit.org/downloads/core_electron
            from scipy.interpolate import UnivariateSpline
            #aeval_spline = UnivariateSpline(rmesh, aeval)
            #psval_spline = UnivariateSpline(rmesh, psval)
            aecore_spline = UnivariateSpline(rmesh, aecore)
            f1 = aecore_spline.derivative(1)
            f2 = aecore_spline.derivative(2)
            header = "%d 0.0 %s %s %s # nr, dummy, symbol, Z, Z_val" % (
                mmax, pseudo.symbol, pseudo.Z, pseudo.Z_val)
            print(header, file=fc_file)
            for ir in range(mmax):
                r = rmesh[ir]
                print(4 * "%.14E  " % (r, aecore[ir], f1(r), f2(r)), file=fc_file)

            print("\n\n<JSON>", file=fc_file)
            print(json.dumps(meta, indent=4), file=fc_file)
            print("</JSON>", file=fc_file)

        if ae_file is not None:
            # Write results to file in "AE" format
            # The files contains:
            #   header with number of points and unknown parameter (not used)
            #   then 2 columns with the radial r coordinate and the AE density at r
            #   See http://www.abinit.org/downloads/all_core_electron
            header = "%d 0.0 %s %s %s # nr, dummy, symbol, Z, Z_val" % (
                mmax, pseudo.symbol, pseudo.Z, pseudo.Z_val)
            print(header, file=ae_file)
            for ir in range(mmax):
                print(2 * "%.14E  " % (rmesh[ir], aecore[ir] + aeval[ir]), file=ae_file)

            print("\n\n<JSON>", file=ae_file)
            print(json.dumps(meta, indent=4), file=ae_file)
            print("</JSON>", file=ae_file)

        if plot:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(rmesh, r2 * aecore, label="AE core * r**2")
            ax.plot(rmesh, r2 * psval, label="PS valence * r**2")
            ax.plot(rmesh, r2 * aeval, label="AE valence * r**2")
            ax.grid(True)
            ax.legend(loc="best")
            plt.show()

        return dict2namedtuple(rmesh=rmesh, psval=psval, aeval=aeval, aecore=aecore)
