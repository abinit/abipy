# coding: utf-8
"""
Classes and functions for parsing the ONCVPSP output file and plotting the results.
"""
from __future__ import annotations

import io
import os
import tempfile
import numpy as np

from typing import Any, Union, List, Optional
from monty.functools import lazy_property
from monty.collections import AttrDict, dict2namedtuple
from monty.os.path import which
from monty.termcolor import cprint
from abipy.core.atom import l2char # NlkState, RadialFunction, RadialWaveFunction,
from abipy.core.mixins import NotebookWriter
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_visible, set_axlims
from abipy.tools.derivatives import finite_diff
from abipy.ppcodes.oncv_parser import OncvParser


class OncvPlotter(NotebookWriter):
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
        parser = OncvParser(filepath)
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
        ylabel = "ATAN(LogDer)"
        ylabel = r"$\phi(E) = \arctan(R * d \psi_E(r)/dr |_R)$"

        self.decorate_ax(ax, xlabel=xlabel, ylabel=ylabel, title="",
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

    def write_notebook(self, nbpath=None):
        raise NotImplementedError("write_notebooks should be tested")
        return oncv_make_open_notebook(self.parser.filepath)


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
from abipy.ppcodes.oncv_parser import OncvParser
onc_parser = OncvParser('%s')""" % outpath),

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


class MultiOncvPlotter(NotebookWriter):
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

        parser = OncvParser(filepath)
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

    def _get_ax_list(self, ax_list, sharex=True, ravel=True):

        num_plots, ncols, nrows = len(self), 1, len(self)

        # Build grid of plots.
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols, # figsize=(8, 8),
                                                sharex=sharex, sharey=False, squeeze=False)
        if ravel:
            ax_list = ax_list.ravel()

        return ax_list, fig, plt

    @add_fig_kwargs
    def plot_atan_logders(self, ax_list=None, with_xlabel=True, xlims=None, ylims=None, fontsize: int = 8, **kwargs):
        """
        Plot arctan of logder on axis ax.

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
        """
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=True)
        n = len(ax_list)

        for i, (ax, (label, plotter)) in enumerate(zip(ax_list, self.items())):
            plotter.plot_atan_logders(ax=ax, with_xlabel=with_xlabel, fontsize=fontsize, show=False)
            ax.set_title(label, fontsize=fontsize)
            #xlims = [-5, 5]
            set_axlims(ax, xlims, "x")
            set_axlims(ax, ylims, "y")
            if i != n - 1:
                set_visible(ax, False, "legend", "xlabel", "ylabel")

        return fig

    @add_fig_kwargs
    def plot_radial_wfs(self, ax_list=None, what="bound_states", fontsize: int = 8, **kwargs):
        """
        Plot AE and PS radial wavefunctions of ax_list.

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            what: "bound_states" or "scattering_states".
        """
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=True)

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
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=True)

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
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=False)

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
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=False)

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
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=False)

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
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=False)

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
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=True)

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
                                               sharex=True, sharey=False, squeeze=False)

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
        ax_list, fig, plt = self._get_ax_list(ax_list, sharex=False)

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
        yield self.plot_kene_vs_ecut(show=False)
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

    def write_notebook(self, nbpath=None):
        raise NotImplementedError("write_notebooks should be tested")
        #return oncv_make_open_notebook(self.parser.filepath)


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
