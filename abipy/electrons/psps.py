# coding: utf-8
"""
Interface to the PSPS.nc file containing form factors computed by ABINIT.
"""
from __future__ import annotations

import os
import numpy as np

from typing import List, Any
from collections import OrderedDict
from monty.bisect import find_gt
from monty.string import list_strings, marquee
from monty.functools import lazy_property
from abipy.iotools import ETSF_Reader
from abipy.core.structure import Structure
from abipy.core.mixins import AbinitNcFile, NotebookWriter
from abipy.abio.robots import Robot
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_visible


def _mklabel(fsym: str, der: int, arg: str) -> str:
    """mklabel(f, 2, x) --> $f''(x)$"""
    if der == 0:
        return "$%s(%s)$" % (fsym, arg)
    else:
        fsym = fsym + "^{" + (der * r"\prime") + "}"
        return "$%s(%s)$" % (fsym, arg)


def _rescale(arr, scale=1.0):
    if scale is None:
        return arr, 0.0

    amax = np.abs(arr).max()
    fact = scale / amax if amax != 0 else 1
    return fact * arr, fact


def dataframe_from_pseudos(pseudos, index=None):
    """
    Build pandas dataframe with the most important info associated to
    a list of pseudos or a list of objects that can be converted into pseudos.

    Args:
        pseudos: List of objects that can be converted to pseudos.
        index: Index of the dataframe.

    Return: pandas Dataframe.
    """
    from abipy.flowtk import PseudoTable
    pseudos = PseudoTable.as_table(pseudos)

    import pandas as pd
    attname = ["Z_val", "l_max", "l_local", "nlcc_radius", "xc", "supports_soc", "type"]
    rows = []
    for p in pseudos:
        row = OrderedDict([(k, getattr(p, k, None)) for k in attname])
        row["ecut_normal"], row["pawecutdg_normal"] = None, None
        if p.has_hints:
            hint = p.hint_for_accuracy(accuracy="normal")
            row["ecut_normal"] = hint.ecut
            if hint.pawecutdg: row["pawecutdg_normal"] = hint.pawecutdg
        rows.append(row)

    return pd.DataFrame(rows, index=index, columns=list(rows[0].keys()) if rows else None)


#def run_and_compare_pseudos_same_element(pseudos):
#    """
#    """
#    from abipy.abio.factories import gs_input
#    pseudos = for pseudo in pseudos
#
#    input= gs_input(structure: Structure, pseudos,
#                    kppa=None, ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
#                    smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None)

class PspsFile(AbinitNcFile, NotebookWriter):
    """
    Netcdf file with the tables used in Abinit to apply the pseudopotential part of the KS Hamiltonian.

    Usage example:

    .. code-block:: python

        with PspsFile("foo_PSPS.nc") as psps:
            psps.plot_tcore_rspace()
    """
    linestyles_der = ["-", "--", '-.', ':', ":", ":"]
    color_der = ["black", "red", "green", "orange", "cyan"]

    @classmethod
    def from_file(cls, filepath: str):
        """Initialize the object from a filepath."""
        return cls(filepath)

    @classmethod
    def from_abinit_run(cls, pseudo):
        """Initialize the object from a filepath."""
        from abipy.flowtk.pseudos import Pseudo
        pseudo = Pseudo.as_pseudo(pseudo)
        structure = Structure.boxed_atom(pseudo)

        from abipy.abio.factories import gs_input
        inp = gs_input(structure, pseudo,
                       kppa=1, ecut=None, pawecutdg=None, accuracy="normal", spin_mode="unpolarized",
                       smearing="fermi_dirac:0.1 eV", charge=0.0)

        inp["prtpsps"] = -1  # Print PSPS.nc and exit
        inp["ecut"] = 100  # Print PSPS.nc and exit
        print(inp)
        filepath = ""

        from abipy.flowtk import AbinitTask
        task = AbinitTask.temp_shell_task(inp) #, workdir=workdir, manager=manager)
        retcode = task.start_and_wait(autoparal=False) #, exec_args=["--dry-run"])
        print(retcode)
        print(task.workdir)

        filepath = os.path.join(task.outdir.path_in("out_PSPS.nc"))
        return cls(filepath)

    def __init__(self, filepath: str):
        super().__init__(filepath)
        self.reader = PspsReader(filepath)

    def close(self) -> None:
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self) -> dict:
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        return {}

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation."""
        lines = []; app = lines.append
        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")

        if verbose > 1:
            app("")
            #app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        return "\n".join(lines)


    @add_fig_kwargs
    def plot(self, **kwargs):
        """
        Driver routine to plot several quantities on the same graph.

        Args:
            ecut_ffnl: Max cutoff energy for ffnl plot (optional)

        Return: |matplotlib-Figure|
        """
        methods = [
            "plot_tcore_rspace",
            "plot_tcore_qspace",
            "plot_ffspl",
            "plot_vlocq",
        ]

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=2, ncols=2,
                                                sharex=False, sharey=False, squeeze=True)

        ecut_ffnl = kwargs.pop("ecut_ffnl", None)
        for m, ax in zip(methods, ax_list.ravel()):
            getattr(self, m)(ax=ax, ecut_ffnl=ecut_ffnl, show=False)

        return fig

    @add_fig_kwargs
    def plot_tcore_rspace(self, ax=None, ders=(0, 1, 2, 3), scale=1.0, rmax=3.0, **kwargs):
        """
        Plot the model core charge and its derivatives in real space.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ders: Tuple used to select the derivatives to be plotted.
            scale: 1.0 if all derivatives should be scaled to 1 else None.
            rmax: Max radius for plot in Bohr. None is full grid is wanted.

        Returns: |matplotlib-Figure|
        """
        if not isinstance(ders, (list, tuple)): ders = [ders]

        rmeshes, coresd = self.reader.read_coresd(rmax=rmax)
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        for rmesh, mcores in zip(rmeshes, coresd):
            for der, values in enumerate(mcores):
                if der not in ders: continue
                yvals, fact, = _rescale(values, scale=scale)
                ax.plot(rmesh, yvals,
                        color=kwargs.get("color", self.color_der[der]),
                        linewidth=kwargs.get("linewidth", 2.0),
                        linestyle=kwargs.get("linestyle", self.linestyles_der[der]),
                        label=_mklabel(r"\tilde{n}_c", der, "r") + " x %.4f" % fact
                       )

        ax.grid(True)
        ax.set_xlabel("r (Bohr)")
        ax.set_title("Model core in r-space")
        if kwargs.get("with_legend", False):
            ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_tcore_qspace(self, ax=None, ders=(0,), with_fact=True, with_qn=0, scale=1.0, **kwargs):
        """
        Plot the model core charge in q space

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ders: Tuple used to select the derivatives to be plotted.
            scale: 1.0 if all derivatives should be scaled to 1 else None.
            with_qn:

        Returns: |matplotlib-Figure|
        """
        if not isinstance(ders, (list, tuple)): ders = [ders]

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        color = kwargs.pop("color", "black")
        linewidth = kwargs.pop("linewidth", 2.0)

        qmesh, tcore_spl = self.reader.read_tcorespl()
        #print(qmesh, tcore_spl)
        ecuts = 2 * (np.pi * qmesh)**2
        lines = []
        for atype, tcore_atype in enumerate(tcore_spl):
            for der, values in enumerate(tcore_atype):
                #if der == 1: der = 2
                if der not in ders: continue
                yvals, fact = _rescale(values, scale=scale)

                label = _mklabel("\\tilde{n}_{c}", der, "q")
                if with_fact: label += " x %.4f" % fact

                line, = ax.plot(ecuts, yvals, color=color, linewidth=linewidth,
                                linestyle=self.linestyles_der[der], label=label)
                lines.append(line)

                if with_qn and der == 0:
                    yvals, fact = _rescale(qmesh * values, scale=scale)
                    line, ax.plot(ecuts, yvals, color=color, linewidth=linewidth,
                                  label=_mklabel("q f", der, "q") + " x %.4f" % fact)

                    lines.append(line)

        ax.grid(True)
        ax.set_xlabel("Ecut (Hartree)")
        ax.set_title("Model core in q-space")
        if kwargs.get("with_legend", False):
            ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_vlocq(self, ax=None, ders=(0,), with_qn=0, with_fact=True, scale=None, **kwargs):
        """
        Plot the local part of the pseudopotential in q space.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ders: Tuple used to select the derivatives to be plotted.
            scale: 1.0 if all derivatives should be scaled to 1 else None.
            with_qn:

        Returns: |matplotlib-Figure|
        """
        if not isinstance(ders, (list, tuple)): ders = [ders]

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        color = kwargs.pop("color", "black")
        linewidth = kwargs.pop("linewidth", 2.0)

        qmesh, vlspl = self.reader.read_vlspl()
        ecuts = 2 * (np.pi * qmesh)**2
        for atype, vl_atype in enumerate(vlspl):
            for der, values in enumerate(vl_atype):
                #if der == 1: der = 2
                if der not in ders: continue

                yvals, fact = _rescale(values, scale=scale)
                label = _mklabel("v_{loc}", der, "q")
                if with_fact: label += " x %.4f" % fact

                ax.plot(ecuts, yvals, color=color, linewidth=linewidth,
                        linestyle=self.linestyles_der[der], label=label)

                if with_qn and der == 0:
                    yvals, fact = _rescale(qmesh * values, scale=scale)
                    ax.plot(ecuts, yvals, color=color, linewidth=linewidth,
                            label="q*f(q) x %2.f" % fact)

        ax.grid(True)
        ax.set_xlabel("Ecut (Hartree)")
        ax.set_title("Vloc(q)")
        if kwargs.get("with_legend", False):
            ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_ffspl(self, ax=None, ecut_ffnl=None, ders=(0,), l_select=None,
                   with_qn=0, with_fact=False, scale=None, **kwargs):
        """
        Plot the nonlocal part of the pseudopotential in q-space.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ecut_ffnl: Max cutoff energy for ffnl plot (optional)
            scale: 1.0 if all derivatives should be scaled to 1 else None.
            ders: Tuple used to select the derivatives to be plotted.
            with_qn:

        Returns: |matplotlib-Figure|
        """
        if not isinstance(ders, (list, tuple)): ders = [ders]

        if l_select is not None:
            if not isinstance(l_select, (list, tuple)): l_select = [l_select]
            #print("l_select:", l_select)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        color = kwargs.pop("color", "black")
        linewidth = kwargs.pop("linewidth", 2.0)

        color_l = {-1: "black", 0: "red", 1: "blue", 2: "green", 3: "orange"}
        linestyles_n = ["solid", '-', '--', '-.', ":"]
        l_seen = set()

        # vlspl has shape [ntypat, 2, mqgrid_vl]
        qmesh, vlspl = self.reader.read_vlspl()
        all_projs = self.reader.read_projectors()

        for itypat, projs_type in enumerate(all_projs):
            # Loop over the projectors for this atom type.
            for p in projs_type:
                if l_select is not None and p.l not in l_select: continue
                print("Printing p.l:", p.l)

                for der, values in enumerate(p.data):
                    if der not in ders: continue
                    #yvals, fact = _rescale(values, scale=scale)

                    label = None
                    if p.l not in l_seen:
                        l_seen.add(p.l)
                        label = _mklabel("v_{nl}", der, "q") + ", l=%d" % p.l

                    stop = len(p.ecuts) +  1
                    if ecut_ffnl is not None:
                        stop = find_gt(p.ecuts, ecut_ffnl)

                    #values = p.ekb * p.values - vlspl[itypat, 0, :]
                    #values = vlspl[itypat, der] + p.sign_sqrtekb * p.values
                    #values = p.sign_sqrtekb * p.values
                    values = p.data[der]

                    #print(values.min(), values.max())
                    ax.plot(p.ecuts[:stop], values[:stop],
                            color=color_l[p.l], linewidth=linewidth,
                            linestyle=linestyles_n[p.n], label=label)

        ax.grid(True)
        ax.set_xlabel("Ecut (Hartree)")
        #ax.set_title("ffnl(q)")
        if kwargs.get("with_legend", False):
            ax.legend(loc="best")

        #ax.axhline(y=0, linewidth=linewidth, color='k', linestyle="solid")
        fig.tight_layout()

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot(show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("psps = abilab.abiopen(%s)" % self.filepath),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class PspsRobot(Robot):
    """
    This robot analyzes the results contained in multiple PSPS.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: PspsRobot
    """
    EXT = "PSPS"

    #@classmethod
    #def from_pseudos(cls, pseudos):
    #    """Initialize the object from a list of filepaths or Pseudos."""
    #
    #    PspsFile.from_pseudo(pseudo)
    #    filepaths = []
    #    return cls.from_files(filepaths)

    def _mkcolor(self, count):
        npseudos = len(self)
        if npseudos <= 2:
            return {0: "red", 1: "blue"}[count]
        else:
            cmap = plt.get_cmap("jet")
            return cmap(float(count) / (1 + len(others)))

    @add_fig_kwargs
    def plot_tcore_rspace(self, ders=(0, 1, 2, 3), with_qn=0, scale=None, fontsize=8, **kwargs):
        """
        Plot the model core charge and its derivatives in r-space.

        Args:
            fontsize: fontsize for subtitles.

        Returns: |matplotlib-Figure|
        """
        nrows, ncols = len(self), len(ders)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        fig.suptitle(f"Model core in r-space")
        for i, (label, psps) in enumerate(self.items()):
            kws = dict(color=self._mkcolor(i), show=False)
            for j, der in enumerate(ders):
                ax = ax_mat[i,j]
                psps.plot_tcore_rspace(ax=ax, ders=der, with_legend=False, scale=scale, **kws)
                ax.set_title(f"$rho_M^{der}(r)$", fontsize=fontsize) if i == 0 else ax.set_title("")
                if i != len(self) - 1: set_visible(ax, False, "xlabel")

        return fig

    @add_fig_kwargs
    def plot_tcore_qspace(self, ders=(0, 1), with_qn=0, scale=None, fontsize=8, **kwargs):
        """
        Plot the model core charge and its derivatives in q-space.

        Args:
            fontsize: fontsize for subtitles.

        Returns: |matplotlib-Figure|
        """
        nrows, ncols = len(self), len(ders)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        fig.suptitle(f"Model core in q-space")
        for i, (label, psps) in enumerate(self.items()):
            kws = dict(color=self._mkcolor(i), show=False)
            for j, der in enumerate(ders):
                ax = ax_mat[i,j]
                psps.plot_tcore_qspace(ax=ax, ders=der, with_qn=with_qn, scale=scale, **kws)
                ax.set_title(f"$rho_M^{der}(q)$", fontsize=fontsize) if i == 0 else ax.set_title("")
                if i != len(self) - 1: set_visible(ax, False, "xlabel")

        return fig

    @add_fig_kwargs
    def plot_vlocq(self, ders=(0, 1), with_qn=0, with_fact=True, scale=None, fontsize=8, **kwargs):
        """
        Plot the local part of the pseudopotential in q space.
        """
        nrows, ncols = len(self), len(ders)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        fig.suptitle(r"$V_{loc}(q)$")
        for i, (label, psps) in enumerate(self.items()):
            kws = dict(color=self._mkcolor(i), show=False)
            for j, der in enumerate(ders):
                ax = ax_mat[i,j]
                psps.plot_vlocq(ax=ax, ders=der, with_qn=with_qn, scale=scale, **kws)
                vloc = "v_{loc}"
                ax.set_title(f"${vloc}^{der}(q)$", fontsize=fontsize) if i == 0 else ax.set_title("")
                if i != len(self) - 1: set_visible(ax, False, "xlabel")

        return fig

    @add_fig_kwargs
    def plot_ffspl(self, ecut_ffnl=None, ders=(0, 1), with_qn=0, scale=None, fontsize=8, **kwargs):
        """
        Plot the nonlocal part of the pseudopotential in q-space.
        """
        l_select = [0, 1, 2]
        #nrows, ncols = len(self), len(ders)
        nrows, ncols = len(self), len(l_select)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                               sharex=True, sharey=False, squeeze=False)

        fig.suptitle(f"ffnl in q-space")
        for i, (label, psps) in enumerate(self.items()):
            kws = dict(color=self._mkcolor(i), show=False)
            for j, l in enumerate(l_select):
            #for j, der in enumerate(ders):
                ax = ax_mat[i,j]
                psps.plot_ffspl(ax=ax, ders=ders, l_select=l, with_qn=with_qn, scale=scale, **kws)
                #ax.set_title(f"$ff_{nl}{der}(q)$", fontsize=fontsize) if i == 0 else ax.set_title("")
                if i != len(self) - 1: set_visible(ax, False, "xlabel")

        #fig.suptitle(r"$V_{loc} in q-space$")

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        verbose = kwargs.get("verbose", 0)
        kws = dict(show=False, tight_layout=True)
        yield self.plot_ffspl(**kws)
        yield self.plot_vlocq(**kws)
        yield self.plot_tcore_rspace(**kws)
        yield self.plot_tcore_qspace(**kws)


    #def yield_plotly_figs(self, **kwargs):  # pragma: no cover
    #    """
    #    This function *generates* a predefined list of plotly figures with minimal input from the user.
    #    """
    #    verbose = kwargs.get("verbose", 0)
    #    for fig in self.yield_ebands_plotly_figs(**kwargs): yield fig
    #    if verbose:
    #        for fig in self.yield_structure_plotly_figs(**kwargs): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.PspsRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class PspsReader(ETSF_Reader):
    """
    This object reads the results stored in the PSPS file produced by ABINIT.
    It provides helper functions to access the most important quantities.
    """
    def __init__(self, filepath):
        super().__init__(filepath)

        # Get important quantities.
        self.usepaw, self.useylm = self.read_value("usepaw"), self.read_value("useylm")
        assert self.usepaw == 0 and self.useylm == 0
        self.ntypat = self.read_dimvalue("ntypat")
        self.lmnmax = self.read_dimvalue("lmnmax")
        self.indlmn = self.read_value("indlmn")

        self.znucl_typat = self.read_value("znucltypat")
        self.zion_typat = self.read_value("ziontypat")

    def read_coresd(self, rmax=None):
        """
        Read the core charges and derivatives for the different types of atoms.

        Args:
            rmax: Maximum radius in Bohr. If None, data on the full grid is returned.

        Returns:
            meshes: List of ntypat arrays. Each array contains the linear meshes in real space.
            coresd: List with nytpat arrays of shape [6, npts].

            (np.zeros. np.zeros) if core charge is not present

        xccc1d[ntypat, n1xccc*(1-usepaw)]

        Norm-conserving psps only
        The component xccc1d(n1xccc,1,ntypat) is the pseudo-core charge
        for each type of atom, on the radial grid. The components
        xccc1d(n1xccc,ideriv,ntypat) give the ideriv-th derivative of the
        pseudo-core charge with respect to the radial distance.
        """

        xcccrc = self.read_value("xcccrc")
        try:
            all_coresd = self.read_value("xccc1d")
        except self.Error:
            # model core may not be present!
            return self.ntypat * [np.linspace(0, 6, num=100)], self.ntypat * [np.zeros((2, 100))]

        npts = all_coresd.shape[-1]
        rmeshes, coresd = [], []
        for itypat, rc in enumerate(xcccrc):
            rvals, step = np.linspace(0, rc, num=npts, retstep=True)
            ir_stop = -1
            if rmax is not None:
                # Truncate the mesh
                ir_stop = min(int(rmax / step), npts) + 1
                #print(rmax, step, ir_stop, npts)

            rmeshes.append(rvals[:ir_stop])
            coresd.append(all_coresd[itypat, :, :ir_stop])

        return rmeshes, coresd

    def read_tcorespl(self):
        """
        Returns:
            qmesh: Linear q-mesh in q-space
            tcorespl: numpy array of shape [ntypat, 2, mqgrid_vl]
                with the pseudo core density in reciprocal space on a regular grid.
                Only if pseudo has_tcore
        """
        return self.read_value("qgrid_vl"), self.read_value("nc_tcorespl")

    def read_vlspl(self):
        """
        Returns:
            qmesh: Linear q-mesh in G-space
            vlspl: numpy array of shape [ntypat, two, mqgrid_vl]
            with the local part of each type of psp in q-space
        """
        return self.read_value("qgrid_vl"), self.read_value("vlspl")

    def read_projectors(self):
        """
        ffspl[ntypat, lnmax, 2, mqgrid_ff]
        Gives, on the radial grid, the {n,l} non-local projectors both for NC and PAW.
        """
        # ekb(dimekb,ntypat*(1-usepaw))
        ekb = self.read_value("ekb")
        qgrid_ff = self.read_value("qgrid_ff")
        #print("qgrid", qgrid_ff.min(), qgrid_ff.max())

        # ffspl[ntypat, lnmax, 2, mqgrid_ff]
        ffspl = self.read_value("ffspl")

        projs = self.ntypat * [None]
        for itypat in range(self.ntypat):
            projs_type = []
            ln_list = self.get_lnlist_for_type(itypat)
            for i, ln in enumerate(ln_list):
                #print(ffspl[itypat, i, :, :])
                p = VnlProjector(itypat, ln, ekb[itypat, i], qgrid_ff, ffspl[itypat, i])
                projs_type.append(p)

            projs[itypat] = projs_type

        return projs

    def get_lnlist_for_type(self, itypat):
        """
        Return a list of (l, n) indices for this atom type.
        """
        # indlmn(6,lmn_size,ntypat) = array giving l,m,n,lm,ln,s for i=lmn
        indlmn_type = self.indlmn[itypat,:,:]

        iln0 = 0; ln_list = []
        for ilmn in range(self.lmnmax):
            iln = indlmn_type[ilmn, 4]
            if iln > iln0:
                iln0 = iln
                l = indlmn_type[ilmn, 0]  # l
                n = indlmn_type[ilmn, 2]  # n
                ln_list.append((l, n))

        return ln_list


class VnlProjector:
    """
    Data and parameters associated to a non-local projector.
    """

    def __init__(self, itypat, ln, ekb, qmesh, data):
        """
        Args:
            itypat: Type atom index (C index >= 0)
            ln: Tuple with l and n.
            ekb: KB energy in Hartree.
            qmesh: Mesh of q-points.
            data: numpy array [2, nqpt]
        """
        self.ln = ln
        self.l, self.n, self.ekb = ln[0], ln[1], ekb
        self.qmesh, self.data = qmesh, data

        assert len(self.qmesh) == len(self.values)
        assert len(self.qmesh) == len(self.der2)

    @property
    def values(self):
        """Values of the projector in q-space."""
        return self.data[0,:]

    @property
    def der2(self):
        """Second order derivative."""
        return self.data[1,:]

    @property
    def ecuts(self):
        """List of cutoff energies in Ha corresponding to self.qmesh."""
        return 2 * (np.pi * self.qmesh)**2

    @property
    def sign_sqrtekb(self):
        return np.sign(self.ekb) * np.sqrt(np.abs(self.ekb))
