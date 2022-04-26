# coding: utf-8
"""
Classes and functions for post-processing the results produced by ONCVPSP.
"""
from __future__ import annotations

import io
import os
import abc
import re
import tempfile
import numpy as np

from collections import namedtuple, OrderedDict
from typing import Any, Union, List, Optional
from monty.functools import lazy_property
from monty.collections import AttrDict, dict2namedtuple
from monty.os.path import which
from monty.termcolor import cprint
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from abipy.tools.derivatives import finite_diff
from abipy.core.atom import NlkState, RadialFunction, RadialWaveFunction


_l2char = {
    "0": "s",
    "1": "p",
    "2": "d",
    "3": "f",
    "4": "g",
    "5": "h",
    "6": "i",
}


def is_integer(s: Any) -> bool:
    """True if object `s` in an integer."""
    try:
        c = float(s)
        return int(c) == c
    except (ValueError, TypeError):
        return False


def decorate_ax(ax, xlabel, ylabel, title, lines, legends, fontsize=8):
    """
    Decorate a `matplotlib` Axis adding xlabel, ylabel, title, grid and legend
    """
    if title: ax.set_title(title, fontsize=fontsize)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.legend(lines, legends, loc="best", fontsize=fontsize, shadow=True)


class PseudoGenDataPlotter:
    """
    Plots the results produced by a pseudopotential generator.
    """
    # TODO: Improve support for fully-relativistic case.
    # List of results supported by the plotter (initialized via __init__)
    all_keys = [
        "radial_wfs",
        "scattering_wfs",
        "projectors",
        "densities",
        "potentials",
        "atan_logders",
        "ene_vs_ecut",
        "kinerr_nlk",
        "rc_l",
        "rc5",
        "lloc",
    ]

    # matplotlib options.
    linewidth, markersize = 2, 2

    linestyle_aeps = dict(ae="solid", ps="dashed")

    markers_aeps = dict(ae=".", ps="o")

    color_l = { 0: "black",
                1: "red",
               -1: "magenta",
                2: "blue",
               -2: "cyan",
                3: "orange",
               -3: "yellow",
              }

    def __init__(self, **kwargs):
        """Store kwargs in self if k is in self.all_keys."""
        import matplotlib.pyplot as _mplt
        self._mplt = _mplt

        for k in self.all_keys:
            setattr(self, k, kwargs.pop(k, {}))

        if kwargs:
            raise ValueError("Unknown keys: %s" % list(kwargs.keys()))

    def keys(self):
        """Iterator with the keys stored in self."""
        return (k for k in self.all_keys if getattr(self, k))

    def plot_key(self, key, ax=None, show=True, **kwargs):
        """
        Plot quantity specified by key.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        # key --> self.plot_key()
        getattr(self, "plot_" + key)(ax=ax, **kwargs)
        if show:
            self._mplt.show()

    def _wf_pltopts(self, l, aeps):
        """Plot options for wavefunctions."""
        return dict(
            color=self.color_l[l],
            linestyle=self.linestyle_aeps[aeps],
            # marker=self.markers_aeps[aeps],
            linewidth=self.linewidth, markersize=self.markersize)

    #@property
    #def has_scattering_states(self):

    def _add_rc_vlines(self, ax):
        """
        Add vertical lines to axis `ax` showing the core radii.
        """
        for l, rc in self.rc_l.items():
            ax.axvline(rc, lw=1, color=self.color_l[l], ls="--") # ls=(0, (1, 10)))

    @add_fig_kwargs
    def plot_atan_logders(self, ax=None, with_xlabel=True, **kwargs):
        """
        Plot arctan of logder on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ae, ps = self.atan_logders.ae, self.atan_logders.ps
        ax, fig, plt = get_ax_fig_plt(ax)

        lines, legends = [], []
        for l, ae_alog in ae.items():
            ps_alog = ps[l]

            # Add pad  to avoid overlapping curves.
            pad = (l + 1) * 1.0

            ae_line, = ax.plot(ae_alog.energies, ae_alog.values + pad, **self._wf_pltopts(l, "ae"))
            ps_line, = ax.plot(ps_alog.energies, ps_alog.values + pad, **self._wf_pltopts(l, "ps"))

            lines.extend([ae_line, ps_line])
            legends.extend(["AE l=%s" % str(l), "PS l=%s" % str(l)])

        xlabel = "Energy (Ha)" if with_xlabel else ""
        decorate_ax(ax, xlabel=xlabel, ylabel="ATAN(LogDer)", title="",
                    #title="ATAN(Log Derivative)",
                    lines=lines, legends=legends)

        return fig

    @add_fig_kwargs
    def plot_radial_wfs(self, ax=None, what="bound_states", **kwargs):
        """
        Plot AE and PS radial wavefunctions on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            what: "bound_states" or "scattering_states".
            lselect: List to select l channels.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        if what == "bound_states":
            ae_wfs, ps_wfs = self.radial_wfs.ae, self.radial_wfs.ps
        elif what == "scattering_states":
            ae_wfs, ps_wfs = self.scattering_wfs.ae, self.scattering_wfs.ps
        else:
            raise ValueError(f"Invalid value for what: {what}")

        lselect = kwargs.get("lselect", [])

        lines, legends = [], []
        for nlk, ae_wf in ae_wfs.items():
            ps_wf, l, k = ps_wfs[nlk], nlk.l, nlk.k
            if l in lselect: continue
            #print(nlk)

            ae_line, = ax.plot(ae_wf.rmesh, ae_wf.values, **self._wf_pltopts(l, "ae"))
            ps_line, = ax.plot(ps_wf.rmesh, ps_wf.values, **self._wf_pltopts(l, "ps"))

            lines.extend([ae_line, ps_line])
            if k is None:
                legends.extend(["AE l=%s" % l, "PS l=%s" % l])
            else:
                legends.extend(["AE l=%s, k=%s" % (l, k), "PS l=%s, k=%s" % (l, k)])

        decorate_ax(ax, xlabel="r (Bohr)", ylabel=r"$\phi(r)$",
                    title="Wave Functions" if what == "bound_states" else "Scattering States",
                    lines=lines, legends=legends)

        self._add_rc_vlines(ax)

        return fig

    @add_fig_kwargs
    def plot_projectors(self, ax=None, **kwargs):
        """
        Plot oncvpsp projectors on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            lselect: List to select l channels
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        lselect = kwargs.get("lselect", [])

        linestyle = {1: "solid", 2: "dashed", 3: "dotted", 4: "dashdot"}
        lines, legends = [], []
        for nlk, proj in self.projectors.items():
            #print(nlk)
            if nlk.l in lselect: continue
            line, = ax.plot(proj.rmesh, proj.values,
                            color=self.color_l.get(nlk.l, 'black'), linestyle=linestyle[nlk.n],
                            linewidth=self.linewidth, markersize=self.markersize)
            lines.append(line); legends.append("Proj %s" % str(nlk))

        decorate_ax(ax, xlabel="r (Bohr)", ylabel="$p(r)$", title="Projectors",
                    lines=lines, legends=legends)

        self._add_rc_vlines(ax)

        return fig

    @add_fig_kwargs
    def plot_densities(self, ax=None, timesr2=False, **kwargs):
        """
        Plot AE, PS and model densities on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        lines, legends = [], []
        for name, rho in self.densities.items():
            d = rho.values if not timesr2 else rho.values * rho.rmesh ** 2
            line, = ax.plot(rho.rmesh, d, linewidth=self.linewidth, markersize=self.markersize)
            lines.append(line); legends.append(name)

        ylabel = "$n(r)$" if not timesr2 else "$r^2 n(r)$"
        decorate_ax(ax, xlabel="r (Bohr)", ylabel=ylabel, title="Charge densities",
                    lines=lines, legends=legends)

        return fig

    @add_fig_kwargs
    def plot_der_densities(self, ax=None, order=1, **kwargs):
        """
        Plot the derivatives of the densitiers on axis ax.
        Used to analyze possible derivative discontinuities

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        from scipy.interpolate import UnivariateSpline

        lines, legends = [], []
        for name, rho in self.densities.items():
            if name != "rhoM": continue
            # Need linear mesh for finite_difference --> Spline input densities on lin_rmesh
            lin_rmesh, h = np.linspace(rho.rmesh[0], rho.rmesh[-1], num=len(rho.rmesh) * 4, retstep=True)
            spline = UnivariateSpline(rho.rmesh, rho.values, s=0)
            lin_values = spline(lin_rmesh)
            vder = finite_diff(lin_values, h, order=order, acc=4)
            line, = ax.plot(lin_rmesh, vder) #, **self._wf_pltopts(l, "ae"))
            lines.append(line)

            legends.append("%s-order derivative of %s" % (order, name))

        decorate_ax(ax, xlabel="r (Bohr)", ylabel="$D^%s \n(r)$" % order, title="Derivative of the charge densities",
                    lines=lines, legends=legends)
        return fig

    @add_fig_kwargs
    def plot_potentials(self, ax=None, **kwargs):
        """
        Plot v_l and v_loc potentials on axis ax

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        lines, legends = [], []
        for l, pot in self.potentials.items():
            line, = ax.plot(pot.rmesh, pot.values, **self._wf_pltopts(l, "ae"))
            lines.append(line)

            if l == -1:
                legends.append("Vloc")
            else:
                legends.append("PS l=%s" % str(l))

        decorate_ax(ax, xlabel="r (Bohr)", ylabel="$v_l(r)$", title="Ion Pseudopotentials",
                    lines=lines, legends=legends)

        color = "k"
        if self.lloc == 4: color = "magenta"
        ax.axvline(self.rc5, lw=1, color=color, ls="--")

        self._add_rc_vlines(ax)

        return fig

    @add_fig_kwargs
    def plot_der_potentials(self, ax=None, order=1, **kwargs):
        """
        Plot the derivatives of vl and vloc potentials on axis ax.
        Used to analyze the derivative discontinuity introduced by the RRKJ method at rc.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        from abipy.tools.derivatives import finite_diff
        from scipy.interpolate import UnivariateSpline
        lines, legends = [], []

        for l, pot in self.potentials.items():
            # Need linear mesh for finite_difference --> Spline input potentials on lin_rmesh
            lin_rmesh, h = np.linspace(pot.rmesh[0], pot.rmesh[-1], num=len(pot.rmesh) * 4, retstep=True)
            spline = UnivariateSpline(pot.rmesh, pot.values, s=0)
            lin_values = spline(lin_rmesh)
            vder = finite_diff(lin_values, h, order=order, acc=4)
            line, = ax.plot(lin_rmesh, vder, **self._wf_pltopts(l, "ae"))
            lines.append(line)

            if l == -1:
                legends.append("%s-order derivative Vloc" % order)
            else:
                legends.append("$s-order derivative PS l=%s" % str(l))

        decorate_ax(ax, xlabel="r (Bohr)", ylabel=r"$D^%s \phi(r)$" % order,
                    title="Derivative of the ion Pseudopotentials",
                    lines=lines, legends=legends)
        return fig

    @add_fig_kwargs
    def plot_ene_vs_ecut(self, ax=None, **kwargs):
        """
        Plot the converge of ene wrt ecut on axis ax.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        lines, legends = [], []
        for l, data in self.ene_vs_ecut.items():
            line, = ax.plot(data.energies, data.values, **self._wf_pltopts(l, "ae"))
            lines.append(line)
            legends.append("Conv l=%s" % str(l))

        for nlk, data in self.kinerr_nlk.items():
            line, = ax.plot(data.ecuts, data.values_ha, **self._wf_pltopts(nlk.l, "ps"))

        decorate_ax(ax, xlabel="Ecut (Ha)", ylabel=r"$\Delta E_{kin}$ (Ha)", title="",
                    #, title="Energy error per electron (Ha)",
                    lines=lines, legends=legends)

        ax.set_yscale("log")
        return fig

    @add_fig_kwargs
    def plot_atanlogder_econv(self, **kwargs):
        """
        Plot atan(logder) and ecut converge on the same figure.
        Returns matplotlib Figure
        """
        fig, ax_list = self._mplt.subplots(nrows=2, ncols=1, sharex=False, squeeze=False)
        ax_list = ax_list.ravel()

        self.plot_atan_logders(ax=ax_list[0], show=False)
        self.plot_ene_vs_ecut(ax=ax_list[1], show=False)

        return fig

    @add_fig_kwargs
    def plot_dens_and_pots(self, **kwargs):
        """Plot densities and potentials on the same figure. Returns matplotlib Figure"""
        fig, ax_list = self._mplt.subplots(nrows=2, ncols=1, sharex=False, squeeze=False)
        ax_list = ax_list.ravel()

        self.plot_densities(ax=ax_list[0], show=False)
        self.plot_potentials(ax=ax_list[1], show=False)

        return fig

    @add_fig_kwargs
    def plot_waves_and_projs(self, **kwargs):
        """Plot ae-ps wavefunctions and projectors on the same figure. Returns matplotlib Figure"""
        lmax = max(nlk.l for nlk in self.radial_wfs.ae.keys())
        fig, ax_list = self._mplt.subplots(nrows=lmax + 1, ncols=2, sharex=True, squeeze=False)

        for l in range(lmax + 1):
            ax_idx = lmax - l
            self.plot_radial_wfs(ax=ax_list[ax_idx][0], lselect=[l], show=False)
            self.plot_projectors(ax=ax_list[ax_idx][1], lselect=[l], show=False)

        return fig

    @add_fig_kwargs
    def plot_den_formfact(self, ecut=120, ax=None, **kwargs):
        """
        Plot the density form factor as a function of ecut in Ha.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Return: matplotlib Figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        lines, legends = [], []
        for name, rho in self.densities.items():
            if name == "rhoC": continue
            form = rho.get_intr2j0(ecut=ecut) / (4 * np.pi)
            line, = ax.plot(form.mesh, form.values, linewidth=self.linewidth, markersize=self.markersize)
            lines.append(line); legends.append(name)

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

        decorate_ax(ax, xlabel="Ecut (Ha)", ylabel="$n(q)$", title="Form factor, l=0 ", lines=lines, legends=legends)
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        verbose = kwargs.get("verbose", 0)

        yield self.plot_radial_wfs(show=False)
        yield self.plot_radial_wfs(what="scattering_states", show=False)
        yield self.plot_atanlogder_econv(show=False)
        yield self.plot_projectors(show=False)
        yield self.plot_potentials(show=False)
        yield #self.plot_der_potentials(show=False)
        yield #for order in [1,2,3,4]:
        yield #    e(self.plot_der_densities(order=order, show=False))
        yield self.plot_densities(show=False)
        yield #self.plot_densities(timesr2=True, show=False)
        yield self.plot_den_formfact(show=False)


class MultiPseudoGenDataPlotter:
    """
    Class for plotting data produced by multiple pp generators on separated plots.

    Usage example:

    .. code-block:: python

        plotter = MultiPseudoGenPlotter()
        plotter.add_psgen("bar.nc", label="bar bands")
        plotter.add_plotter("foo.nc", label="foo bands")
        plotter.plot()
    """
    _LINE_COLORS = ["b", "r"]
    _LINE_STYLES = ["-", ":", "--", "-."]
    _LINE_WIDTHS = [2]

    def __init__(self):
        self._plotters_odict = OrderedDict()

    def __len__(self) -> int:
        return len(self._plotters_odict)

    @property
    def plotters(self):
        """"List of registered `Plotters`."""
        return list(self._plotters_odict.values())

    @property
    def labels(self) -> list:
        """List of labels."""
        return list(self._plotters_odict.keys())

    def keys(self) -> list:
        """List of strings with the quantities that can be plotted."""
        keys_set = set()
        for plotter in self.plotters:
            keys_set.update(plotter.keys())
        keys_set = list(keys_set)

        return list(sorted(keys_set))

    def iter_lineopt(self):
        """Generates style options for lines."""
        import itertools
        for o in itertools.product(self._LINE_WIDTHS,  self._LINE_STYLES, self._LINE_COLORS):
            yield {"linewidth": o[0], "linestyle": o[1], "color": o[2]}

    def add_psgen(self, label: str, psgen) -> None:
        """Add a plotter of class plotter_class from a `PseudoGenerator` instance."""
        oparser = psgen.parse_output()
        self.add_plotter(label, oparser.make_plotter())

    def add_plotter(self, label: str, plotter) -> None:
        """
        Adds a plotter.

        Args:
            label: label for the plotter. Must be unique.
            plotter: :class:`PseudoGenDataPlotter` object.
        """
        if label in self.labels:
            raise ValueError("label %s is already in %s" % (label, self.labels))

        self._plotters_odict[label] = plotter

    @add_fig_kwargs
    def plot_key(self, key, **kwargs):
        r"""
        Plot the band structure and the DOS.

        Args:
            klabels: dictionary whose keys are tuple with the reduced
                coordinates of the k-points. The values are the labels.
                e.g. klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        sharex          True if subplots should share the x axis
        ==============  ==============================================================

        Returns:
            matplotlib figure.
        """
        import matplotlib.pyplot as plt

        # Build grid of plots.
        fig, ax_list = plt.subplots(nrows=len(self), ncols=1, sharex=kwargs.pop("sharex", True), squeeze=False)
        ax_list = ax_list.ravel()

        for ax in ax_list:
            ax.grid(True)

        # Plot key on the each axis in ax_list.
        lines, legends = [], []
        #my_kwargs, opts_label = kwargs.copy(), {}
        i = -1
        for (label, plotter), lineopt in zip(self._plotters_odict.items(), self.iter_lineopt()):
            i += 1
            plotter.plot_key(key, ax=ax_list[i], show=False)

        return fig


class PseudoGenResults(AttrDict):

    _KEYS = [
        "max_ecut",
        "max_atan_logder_l1err",
    ]

    def __init__(self, *args, **kwargs):
        super(PseudoGenResults, self).__init__(*args, **kwargs)
        for k in self._KEYS:
            if k not in self:
                self[k] = None


class AtanLogDer(namedtuple("AtanLogDer", "l, energies, values")):

    @property
    def to_dict(self) -> dict:
        return dict(
            l=self.l,
            energies=list(self.energies),
            values=list(self.values))


class PseudoGenOutputParserError(Exception):
    """Exceptions raised by OuptputParser."""


class PseudoGenOutputParser:
    """
    Abstract class defining the interface that must be provided
    by the parsers used to extract results from the output file of
    a pseudopotential generator a.k.a. ppgen

    Attributes:

        errors: List of strings with errors reported by the pp generator
        warnings: List of strings with the warnings reported by the pp generator.
        results:
    """
    __metaclass__ = abc.ABCMeta

    Error = PseudoGenOutputParserError

    def __init__(self, filepath: str) -> None:
        self.filepath = os.path.abspath(filepath)
        self.run_completed = False
        self._errors = []
        self._warnings = []
        self._results = None

    @property
    def errors(self) -> List[str]:
        """List of strings with possible errors reported by the generator at run-time."""
        return self._errors

    @property
    def warnings(self) -> List[str]:
        """List of strings with possible errors reported by the generator at run-time."""
        return self._warnings

    @property
    def results(self):
        return self._results

    @abc.abstractmethod
    def get_results(self):
        """
        Return the most important results of the run in a dictionary.
        Set self.results attribute
        """

    @abc.abstractmethod
    def get_input_str(self) -> str:
        """Returns a string with the input file."""

    @abc.abstractmethod
    def get_pseudo_str(self) -> str:
        """Returns a string with the pseudopotential file."""


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
        p.radial_wavefunctions

        # To plot data with matplotlib.
        p = parser.make_plotter()
        p.plot_atanlogder_econv()

    """
    # TODO Add fully-relativistic case.

    # Used to store ae and pp quantities (e.g wavefunctions) in a single object.
    AePsNamedTuple = namedtuple("AePsNamedTuple", "ae, ps")

    # Object returned by self._grep
    GrepResults = namedtuple("GrepResults", "data, start, stop")

    # Object storing the final results.
    Results = PseudoGenResults

    # Class used to instantiate the plotter.
    Plotter = PseudoGenDataPlotter

    def scan(self, verbose: int = 0) -> None:
        """
        Scan the output file, set `run_completed` attribute.

        Raises: self.Error if invalid file.
        """
        try:
            return self._scan(verbose=verbose)
        except Exception as exc:
            print(f"Exception while parsing output file: {self.filepath}")
            #print(self.lines)
            raise exc

    def _scan(self, verbose: int = 0) -> None:
        if not os.path.exists(self.filepath):
            raise self.Error("File %s does not exist" % self.filepath)

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

        #if self.errors:
        #    return 1

        # scalar-relativistic version 2.1.1, 03/26/2014
        # scalar-relativistic version 3.0.0 10/10/2014
        # toks, self.gendate = self.lines[1].split(",")
        # toks = toks.split()
        toks = self.lines[1].replace(",", " ").split()
        self.gendate = toks.pop(-1)
        self.calc_type, self.version = toks[0], toks[2]

        #print("version", self.version)
        self.major_version, self.minor_version, self.patch_level = tuple(map(int, self.version.split(".")))[:3]

        fullr = self.calc_type not in ["scalar-relativistic", "non-relativistic"]
        if fullr:
            print("Parsing of wavefunctions and projectors is not yet coded for fully-relativisic pseudos")

        # Read configuration (not very robust because we assume the user didn't change the template but oh well)
        header = "# atsym  z    nc    nv    iexc   psfile"
        for i, line in enumerate(self.lines):
            if line.startswith("# atsym"):
                values = self.lines[i+1].split()
                keys = header[1:].split()
                # assert len(keys) == len(values)
                # Store them in self.
                for k, v in zip(keys, values):
                    setattr(self, k, v)
                break

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
        if fullr and self.major_version <= 3:
            header = "#   n    l    f        l+1/2             l-1/2"
        else:
            header = "#   n    l    f        energy (Ha)"

        # Convert nc and nv to int.
        self.nc, self.nv = int(self.nc), int(self.nv)
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
                    core.append(n + _l2char[l] + "^%s" % f)
                self.core = "$" + " ".join(core) + "$"

                beg, valence = i + nc + 1, []
                for v in range(nv):
                    #print("lines[beg+v]", self.lines[beg+v])
                    n, l, f = self.lines[beg+v].split()[:3]
                    if is_integer(f):
                        f = str(int(float(f)))
                    else:
                        #print(f)
                        f = "%.1f" % float(f)

                    valence.append(n + _l2char[l] + "^{%s}" % f)

                self.valence = "$" + " ".join(valence) + "$"
                #print("core", self.core, "valence",self.valence)
                break
        else:
            raise self.Error(f"Cannot find header:\n`{header}`\nin output file {self.filepath}")

        # Read lmax (not very robust because we assume the user didn't change the template but oh well)
        header = "# lmax"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                self.lmax = int(self.lines[i+1])
                break
        else:
            raise self.Error(f"Cannot find line with `#lmax` in: {self.filepath}")

        # Get core radii as a function of l from the output file
        self.rc_l = {}
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
                    self.rc_l[l] = rc
                    nxt += 1

        #print("rc_l", self.rc_l)
        if not self.rc_l:
            raise self.Error(f"Cannot find magic line starting with `{header}` in: {self.filepath}")


        # Get pseudization radius for the local part
        header = "# lloc, lpopt,  rc(5),   dvloc0"
        self.rc5 = None
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                ln = self.lines[i + 1]
                #print("line: ", line, "\nrc line: ", ln)
                tokens = ln.split()
                self.lloc = int(tokens[0])
                self.rc5 = float(tokens[2])
                break

        if self.rc5 is None:
            raise self.Error(f"Cannot find magic line starting with `{header}` in: {self.filepath}")

        self.kinerr_nlk = {}
        if self.major_version > 3:
            # Calculating optimized projector #   1
            #
            #  for l=   0

            re_start = re.compile(r"^Calculating optimized projector #\s+(?P<iproj>\d+)")
            #header = "Calculating optimized projector #"
            # FIXME: In FR mode, we have
            #Calculating first optimized projector for l=   0
            #Calculating second optimized projector for l=   0
        else:
            # Calculating first optimized projector for l=   0
            re_start = re.compile(r"^Calculating (?P<iproj>(first|second)) optimized projector for l=\s+(?P<l>\d+)")

        for i, line in enumerate(self.lines):
            m = re_start.match(line)
            #if line.startswith(header):
            #iproj = int(line.replace(header, "").strip())
            if m:
                if self.major_version > 3:
                    # for l=   0
                    iproj = int(m.group("iproj"))
                    l = int(self.lines[i+2].split("=")[-1].strip())
                    k = None
                else:
                    iproj = m.group("iproj")
                    iproj = {"first": 0, "second": 1}[iproj]
                    l = int(m.group("l"))
                    k = None

                nlk = NlkState(n=iproj, l=l, k=k)
                continue

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
                    if not tokens: break
                    #print("tokens:", tokens)
                    err_ha, err_ev, ecut = map(float, tokens)
                    #print(err_ha, ecut)
                    values_ha.append(err_ha)
                    ecuts.append(ecut)

                self.kinerr_nlk[nlk] = dict2namedtuple(ecuts=ecuts, values_ha=values_ha)

        if not self.kinerr_nlk:
            raise self.Error(f"Cannot parse convergence profile from: {self.filepath}")


    def __str__(self) -> str:
        """String representation."""
        lines = []
        app = lines.append

        if hasattr(self, "calc_type"):
            app("%s, oncvpsp version: %s, date: %s" % (self.calc_type, self.version, self.gendate))
            app("oncvpsp calculation: %s: " % self.calc_type)
            app("completed: %s" % self.run_completed)
        else:
            app("Object is empty. Call scan method to analyze output file")

        return "\n".join(lines)

    @property
    def fully_relativistic(self) -> bool:
        """True if fully-relativistic calculation."""
        return self.calc_type == "fully-relativistic"

    @lazy_property
    def potentials(self) -> dict:
        """Radial functions with the non-local and local potentials."""
        #radii, charge, pseudopotentials (ll=0, 1, lmax)
        #!p   0.0099448   4.7237412  -7.4449470 -14.6551019
        vl_data = self._grep("!p").data
        lmax = len(vl_data[0]) - 3
        assert lmax == self.lmax

        # From 0 up to lmax
        ionpots_l = {}
        for l in range(lmax + 1):
            ionpots_l[l] = RadialFunction("Ion Pseudopotential, l=%d" % l, vl_data[:, 0], vl_data[:, 2+l])

        # Local part is stored with l == -1 if lloc=4, not present if lloc=l
        vloc = self._grep("!L").data
        if vloc is not None:
            ionpots_l[-1] = RadialFunction("Local part, l=%d" % -1, vloc[:, 0], vloc[:, 1])

        return ionpots_l

    @lazy_property
    def densities(self) -> dict:
        """
        Dictionary with charge densities on the radial mesh.
        """
        # radii, charge, core charge, model core charge
        # !r   0.0100642   4.7238866  53.4149287   0.0000000
        rho_data = self._grep("!r").data

        return dict(
            rhoV=RadialFunction("Valence charge", rho_data[:, 0], rho_data[:, 1]),
            rhoC=RadialFunction("Core charge", rho_data[:, 0], rho_data[:, 2]),
            rhoM=RadialFunction("Model charge", rho_data[:, 0], rho_data[:, 3]))

    @lazy_property
    def radial_wfs(self):
        """
        Read and set the radial wavefunctions.
        """
        return self._get_radial_wavefunctions(what="bound_states")

    @lazy_property
    def scattering_wfs(self):
        """
        Read and set the scattering wavefunctions.
        """
        return self._get_radial_wavefunctions(what="scattering_states")

    def _get_radial_wavefunctions(self, what: str):
        # scalar-relativistic
        #n= 1,  l= 0, all-electron wave function, pseudo w-f
        #
        #&     0    0.009945   -0.092997    0.015273

        # Fully-relativistic
        #n= 1,  l= 0  kap=-1, all-electron wave function, pseudo w-f
        #
        #&     0    0.009955    0.066338    0.000979

        # scalar-relativistic
        #scattering, iprj= 2,  l= 1, all-electron wave function, pseudo w-f

        ae_waves, ps_waves = OrderedDict(), OrderedDict()

        beg = 0
        while True:
            g = self._grep("&", beg=beg)
            if g.data is None: break
            beg = g.stop + 1

            header = self.lines[g.start-2]
            #print(header)

            if what == "bound_states":
                if header.startswith("scattering,"):
                    continue
                    print(f"ignoring header {header}")

            elif what == "scattering_states":
                if not header.startswith("scattering,"):
                    continue
                header = header.replace("scattering,", "")
                print(header)
            else:
                raise ValueError(f"Invalid value of what: `{what}`")

            if not self.fully_relativistic:
                # n= 1,  l= 0, all-electron wave function, pseudo w-f
                n, l = header.split(",")[0:2]
                n = int(n.split("=")[1])
                l = int(l.split("=")[1])
                k = None
            else:
                # n= 1,  l= 0  kap=-1, all-electron wave function, pseudo w-f
                n, lk = header.split(",")[0:2]
                n = int(n.split("=")[1])
                toks = lk.split("=")
                l = int(toks[1].split()[0])
                k = int(toks[-1])

            nlk = NlkState(n=n, l=l, k=k)
            #print("Got state: %s" % str(nlk))

            rmesh = g.data[:, 1]
            ae_wf = g.data[:, 2]
            ps_wf = g.data[:, 3]

            #assert nlk not in ae_waves
            ae_waves[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, ae_wf)
            ps_waves[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, ps_wf)

        return self.AePsNamedTuple(ae=ae_waves, ps=ps_waves)

    @lazy_property
    def projectors(self):
        """Read the projector wave functions."""
        #n= 1 2  l= 0, projecctor pseudo wave functions, well or 2nd valence
        #
        #@     0    0.009945    0.015274   -0.009284
        projectors_nlk = OrderedDict()
        beg = 0
        magic = "@"
        if self.major_version > 3: magic = "!J"
        while True:
            g = self._grep(magic, beg=beg)
            if g.data is None: break
            beg = g.stop + 1
            # TODO: Get n, l, k from header.
            header = self.lines[g.start-2]

            #n= 1 2  l= 0, projecctor pseudo wave functions, well or 2nd valence
            # n= 1 2  l= 0  kap=-1, projecctor pseudo wave functions, well or 2nd valence

            rmesh = g.data[:, 1]
            l = int(g.data[0, 0])
            #print("header", header)
            #print("g.data", g.data[0, 0])
            #print("l",l)

            for n in range(len(g.data[0]) - 2):
                nlk = NlkState(n=n+1, l=l, k=None)
                #print("Got projector with: %s" % str(nlk))
                assert nlk not in projectors_nlk
                projectors_nlk[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, g.data[:, n+2])

        return projectors_nlk

    @lazy_property
    def atan_logders(self):
        """Atan of the log derivatives for different l-values."""
        #log derivativve data for plotting, l= 0
        #atan(r * ((d psi(r)/dr)/psi(r))), r=  1.60
        #l, energy, all-electron, pseudopotential
        #
        #!      0    2.000000    0.706765    0.703758
        ae_atan_logder_l, ps_atan_logder_l = OrderedDict(), OrderedDict()

        lstop = self.lmax + 1
        if self.major_version > 3:
            lstop = min(self.lmax + 2, 4)

        for l in range(lstop):
            data = self._grep(tag="!      %d" % l).data
            assert l == int(data[0, 0])
            ae_atan_logder_l[l] = AtanLogDer(l=l, energies=data[:, 1], values=data[:, 2])
            ps_atan_logder_l[l] = AtanLogDer(l=l, energies=data[:, 1], values=data[:, 3])

        return self.AePsNamedTuple(ae=ae_atan_logder_l, ps=ps_atan_logder_l)

    @lazy_property
    def ene_vs_ecut(self) -> dict:
        """Convergence of energy versus ecut for different l-values."""
        #convergence profiles, (ll=0,lmax)
        #!C     0    5.019345    0.010000
        #...
        #!C     1   19.469226    0.010000
        class ConvData(namedtuple("ConvData", "l energies values")):
            @property
            def to_dict(self):
                return dict(
                    l=self.l,
                    energies=list(self.energies),
                    values=list(self.values))

        conv_l = OrderedDict()

        for l in range(self.lmax + 1):
            data = self._grep(tag="!C     %d" % l).data
            conv_l[l] = ConvData(l=l, energies=data[:, 1], values=data[:, 2])

        return conv_l

    @lazy_property
    def hints(self):
        """Hints for the cutoff energy."""
        # Extract the hints
        hints = 3 * [-np.inf]
        ene_vs_ecut = self.ene_vs_ecut
        for i in range(3):
            for l in range(self.lmax + 1):
                hints[i] = max(hints[i], ene_vs_ecut[l].energies[-i-1])
        hints.reverse()

        # print("hints:", hints)
        # Truncate to the nearest int
        hints = [np.rint(h) for h in hints]

        hints = dict(
            low={"ecut": hints[0], "pawecutdg": hints[0]},
            normal={"ecut": hints[1], "pawecutdg": hints[1]},
            high={"ecut": hints[2], "pawecutdg": hints[2]})

        return hints

    def get_results(self):
        """"
        Return the most important results reported by the pp generator.
        Set the value of self.results
        """
        #if not self.run_completed:
        #    self.Results(info="Run is not completed")

        # Get the ecut needed to converge within ... TODO
        max_ecut = 0.0
        for l in range(self.lmax + 1):
            max_ecut = max(max_ecut, self.ene_vs_ecut[l].energies[-1])

        # Compute the l1 error in atag(logder)
        from scipy.integrate import cumtrapz
        max_l1err = 0.0
        for l in range(self.lmax + 1):
            f1, f2 = self.atan_logders.ae[l], self.atan_logders.ps[l]

            adiff = np.abs(f1.values - f2.values)
            integ = cumtrapz(adiff, x=f1.energies) / (f1.energies[-1] - f1.energies[0])
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

        self._results = self.Results(
            max_ecut=max_ecut, max_atan_logder_l1err=max_l1err,
            herm_err=herm_err, max_psexc_abserr=max_psexc_abserr)

        return self._results

    def find_string(self, s: str) -> int:
        """
        Returns the index of the first line containing string s.
        Raise self.Error if s cannot be found.
        """
        for i, line in enumerate(self.lines):
            if s in line:
                return i
        else:
            raise self.Error("Cannot find %s in lines" % s)

    def get_input_str(self) -> str:
        """String with the input file."""
        try:
            # oncvpsp 3.2.3
            i = self.find_string("<INPUT>")
            j = self.find_string("</INPUT>")
            return "\n".join(self.lines[i+1:j]) + "\n"
        except self.Error:
            #raise
            i = self.find_string("Reference configufation results")
            return "\n".join(self.lines[:i])

    def get_psp8_str(self) -> Union[str, None]:
        """
        Return string with the pseudopotential data in psp8 format.
        None if field is not present.
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
        None if field is not present.
        """
        start, stop = None, None
        for i, line in enumerate(self.lines):
            if "Begin PSP_UPF" in line: start = i
            if start is not None and 'END_PSP' in line:
                stop = i
                break
        if start is None and stop is None: return None
        return "\n".join(self.lines[start+1:stop])

    def get_pseudo_str(self):
        """Return string with the pseudopotential data."""

        psdata = self.get_psp8_str()
        if psdata is not None:
            return psdata
        psdata = self.get_upf_str()
        if psdata is not None:
            return psdata

        raise ValueError("Cannot find neither PSPCODE8 not PSP_UPF tag in output file")

    def make_plotter(self):
        """Builds an instance of :class:`PseudoGenDataPlotter`."""
        kwargs = {k: getattr(self, k) for k in self.Plotter.all_keys}
        return self.Plotter(**kwargs)

    @property
    def to_dict(self) -> dict:
        """
        Returns a dictionary with the radial functions and the other
        important results produced by ONCVPSP in JSON format.
        """
        # Dimensions and basic info.
        jdict = dict(
            lmax=self.lmax,
            ncore=self.nc,
            nvalence=self.nv,
            calc_type=self.calc_type
        )

        # List of radial wavefunctions (AE and PS)
        jdict["radial_wfs"] = d = {}
        d["ae"] = [wave.to_dict for wave in self.radial_wfs.ae.values()]
        d["ps"] = [wave.to_dict for wave in self.radial_wfs.ps.values()]

        # List of projectors
        jdict["projectors"] = [proj.to_dict for proj in self.projectors.values()]

        # Charge densities
        jdict["densities"] = dict(
            rhoV=self.densities["rhoV"].to_dict,
            rhoC=self.densities["rhoC"].to_dict,
            rhoM=self.densities["rhoM"].to_dict)

        # Logders (AE and PS)
        jdict["atan_logders"] = d = {}
        d["ae"] = [f.to_dict for f in self.atan_logders.ae.values()]
        d["ps"] = [f.to_dict for f in self.atan_logders.ps.values()]

        # Convergence of the different l-channels as function of ecut.
        jdict["ene_vs_ecut"] = [f.to_dict for f in self.ene_vs_ecut.values()]

        return jdict

    def _grep(self, tag: str, beg: int = 0):
        """
        Finds the first field in the file with the specified tag.
        beg gives the initial position in the file.

        Returns:
            :class:`GrepResult` object
        """
        data, stop, intag = [], None, -1

        if beg >= len(self.lines):
            raise ValueError("beg > len(lines)")

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
            return self.GrepResults(data=None, start=intag, stop=stop)
        else:
            return self.GrepResults(data=np.array(data), start=intag, stop=stop)

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
        print("Working in %s" % workdir)

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
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser
onc_parser = OncvOutputParser('%s')""" % outpath),

        nbf.new_code_cell("""\
# Parse the file and build the plotter
onc_parser.scan()
if not onc_parser.run_completed:
    raise RuntimeError("Cannot parse output file")

plotter = onc_parser.make_plotter()"""),

        nbf.new_markdown_cell(r"# AE and PS radial wavefunctions $\phi(r)$:"),
        nbf.new_code_cell("fig = plotter.plot_radial_wfs(show=False)"),

        nbf.new_markdown_cell("# Arctan of the logarithmic derivatives:"),
        nbf.new_code_cell("fig = plotter.plot_atan_logders(show=False)"),

        nbf.new_markdown_cell("# Convergence in $G$-space estimated by ONCVPSP:"),
        nbf.new_code_cell("fig = plotter.plot_ene_vs_ecut(show=False)"),

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


def psp8_get_densities(path, fc_file=None, ae_file=None, plot=False):
    """
    Extract AE-core, AE-valence and PS-valence from from a psp8 file.
    Optionally write `.fc` and `.AE` file in format suitable for AIM and cut3d Hirshfeld code

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
    from pseudo_dojo.core.pseudos import Pseudo
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
                print("For l=%s, npl=%s" % (l, npl), "wrong line", line)
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
