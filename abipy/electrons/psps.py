# coding: utf-8
"""PSPS file with tabulated data."""
import numpy as np

from collections import OrderedDict
from monty.bisect import find_gt
from monty.functools import lazy_property
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt
from abipy.iotools import ETSF_Reader
from abipy.core.mixins import AbinitNcFile

import logging
logger = logging.getLogger(__name__)


def mklabel(fsym, der, arg):
    """mklabel(f, 2, x) --> $f''(x)$"""
    if der == 0:
        return "$%s(%s)$" % (fsym, arg)
    else:
        fsym = fsym + "^{" + (der * r"\prime") + "}"
        return "$%s(%s)$" % (fsym, arg)


def rescale(arr, scale=1.0):
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


class PspsFile(AbinitNcFile):
    """
    Netcdf file with the tables used in Abinit to apply the
    pseudopotential part of the KS Hamiltonian.

    Usage example:

    .. code-block:: python

        with PspsFile("foo_PSPS.nc") as psps:
            psps.plot_tcore_rspace()
    """
    linestyles_der = ["-", "--", '-.', ':', ":", ":"]
    color_der = ["black", "red", "green", "orange", "cyan"]

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a Netcdf file"""
        return cls(filepath)

    def __init__(self, filepath):
        super().__init__(filepath)
        self.reader = r = PspsReader(filepath)

    def close(self):
        """Close the file."""
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        return {}

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
    def plot_tcore_rspace(self, ax=None, ders=(0, 1, 2, 3), rmax=3.0,  **kwargs):
        """
        Plot the model core and its derivatives in real space.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ders: Tuple used to select the derivatives to be plotted.
            rmax: Max radius for plot in Bohr. None is full grid is wanted.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        linewidth = kwargs.pop("linewidth", 2.0)
        rmeshes, coresd = self.reader.read_coresd(rmax=rmax)

        scale = None
        scale = 1.0
        for rmesh, mcores in zip(rmeshes, coresd):
            for der, values in enumerate(mcores):
                if der not in ders: continue
                yvals, fact, = rescale(values, scale=scale)
                ax.plot(rmesh, yvals, color=self.color_der[der], linewidth=linewidth,
                        linestyle=self.linestyles_der[der],
                        label=mklabel("\\tilde{n}_c", der, "r") + " x %.4f" % fact)

        ax.grid(True)
        ax.set_xlabel("r [Bohr]")
        ax.set_title("Model core in r-space")
        if kwargs.get("with_legend", False): ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_tcore_qspace(self, ax=None, ders=(0,), with_fact=True, with_qn=0, **kwargs):
        """
        Plot the model core in q space

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ders: Tuple used to select the derivatives to be plotted.
            with_qn:

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        color = kwargs.pop("color", "black")
        linewidth = kwargs.pop("linewidth", 2.0)

        qmesh, tcore_spl = self.reader.read_tcorespl()
        #print(qmesh, tcore_spl)
        ecuts = 2 * (np.pi * qmesh)**2
        lines = []
        scale = 1.0
        scale = None
        for atype, tcore_atype in enumerate(tcore_spl):
            for der, values in enumerate(tcore_atype):
                if der == 1: der = 2
                if der not in ders: continue
                yvals, fact = rescale(values, scale=scale)

                label = mklabel("\\tilde{n}_{c}", der, "q")
                if with_fact: label += " x %.4f" % fact

                line, = ax.plot(ecuts, yvals, color=color, linewidth=linewidth,
                                linestyle=self.linestyles_der[der], label=label)
                lines.append(line)

                if with_qn and der == 0:
                    yvals, fact = rescale(qmesh * values, scale=scale)
                    line, ax.plot(ecuts, yvals, color=color, linewidth=linewidth,
                                  label=mklabel("q f", der, "q") + " x %.4f" % fact)

                    lines.append(line)

        ax.grid(True)
        ax.set_xlabel("Ecut [Hartree]")
        ax.set_title("Model core in q-space")
        if kwargs.get("with_legend", False): ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_vlocq(self, ax=None, ders=(0,), with_qn=0, with_fact=True, **kwargs):
        """
        Plot the local part of the pseudopotential in q space.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ders: Tuple used to select the derivatives to be plotted.
            with_qn:

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        color = kwargs.pop("color", "black")
        linewidth = kwargs.pop("linewidth", 2.0)

        qmesh, vlspl = self.reader.read_vlspl()
        ecuts = 2 * (np.pi * qmesh)**2
        scale = 1.0
        scale = None
        for atype, vl_atype in enumerate(vlspl):
            for der, values in enumerate(vl_atype):
                if der == 1: der = 2
                if der not in ders: continue

                yvals, fact = rescale(values, scale=scale)
                label = mklabel("v_{loc}", der, "q")
                if with_fact: label += " x %.4f" % fact

                ax.plot(ecuts, yvals, color=color, linewidth=linewidth,
                        linestyle=self.linestyles_der[der], label=label)

                if with_qn and der == 0:
                    yvals, fact = rescale(qmesh * values, scale=scale)
                    ax.plot(ecuts, yvals, color=color, linewidth=linewidth,
                            label="q*f(q) x %2.f" % fact)

        ax.grid(True)
        ax.set_xlabel("Ecut [Hartree]")
        ax.set_title("Vloc(q)")
        if kwargs.get("with_legend", False): ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_ffspl(self, ax=None, ecut_ffnl=None, ders=(0,), with_qn=0, with_fact=False, **kwargs):
        """
        Plot the nonlocal part of the pseudopotential in q-space.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            ecut_ffnl: Max cutoff energy for ffnl plot (optional)
            ders: Tuple used to select the derivatives to be plotted.
            with_qn:

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        color = kwargs.pop("color", "black")
        linewidth = kwargs.pop("linewidth", 2.0)

        color_l = {-1: "black", 0: "red", 1: "blue", 2: "green", 3: "orange"}
        linestyles_n = ["solid", '-', '--', '-.', ":"]
        scale = None
        l_seen = set()

        qmesh, vlspl = self.reader.read_vlspl()

        all_projs = self.reader.read_projectors()
        for itypat, projs_type in enumerate(all_projs):
            # Loop over the projectors for this atom type.
            for p in projs_type:
                for der, values in enumerate(p.data):
                    if der == 1: der = 2
                    if der not in ders: continue
                    #yvals, fact = rescale(values, scale=scale)
                    label = None
                    if p.l not in l_seen:
                        l_seen.add(p.l)
                        label = mklabel("v_{nl}", der, "q") + ", l=%d" % p.l

                    stop = len(p.ecuts)
                    if ecut_ffnl is not None:
                        stop = find_gt(p.ecuts, ecut_ffnl)

                    #values = p.ekb * p.values - vlspl[itypat, 0, :]
                    values = vlspl[itypat, 0, :] + p.sign_sqrtekb * p.values

                    #print(values.min(), values.max())
                    ax.plot(p.ecuts[:stop], values[:stop], color=color_l[p.l], linewidth=linewidth,
                            linestyle=linestyles_n[p.n], label=label)

        ax.grid(True)
        ax.set_xlabel("Ecut [Hartree]")
        ax.set_title("ffnl(q)")
        if kwargs.get("with_legend", False): ax.legend(loc="best")

        ax.axhline(y=0, linewidth=linewidth, color='k', linestyle="solid")
        fig.tight_layout()

        return fig

    @add_fig_kwargs
    def compare(self, others, **kwargs):
        """Produce matplotlib plot comparing self with another list of pseudos ``others``."""
        if not isinstance(others, (list, tuple)):
            others = [others]

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=2, ncols=2,
                                                sharex=False, sharey=False, squeeze=True)
        ax_list = ax_list.ravel()
        #fig.suptitle("%s vs %s" % (self.basename, ", ".join(o.basename for o in others)))

        def mkcolor(count):
            npseudos = 1 + len(others)
            if npseudos <= 2:
                return {0: "red", 1: "blue"}[count]
            else:
                cmap = plt.get_cmap("jet")
                return cmap(float(count) / (1 + len(others)))

        ic = 0; ax = ax_list[ic]
        self.plot_tcore_rspace(ax=ax, color=mkcolor(0), show=False, with_legend=False)
        for count, other in enumerate(others):
            other.plot_tcore_rspace(ax=ax, color=mkcolor(count+1), show=False, with_legend=False)

        ic += 1; ax = ax_list[ic]
        self.plot_tcore_qspace(ax=ax, with_qn=0, color=mkcolor(0), show=False)
        for count, other in enumerate(others):
            other.plot_tcore_qspace(ax=ax, with_qn=0, color=mkcolor(count+1), show=False)

        ic += 1; ax = ax_list[ic]
        self.plot_vlocq(ax=ax, with_qn=0, color=mkcolor(0), show=False)
        for count, other in enumerate(others):
            other.plot_vlocq(ax=ax, with_qn=0, color=mkcolor(count+1), show=False)

        ic += 1; ax = ax_list[ic]
        self.plot_ffspl(ax=ax, with_qn=0, color=mkcolor(0), show=False)
        for count, other in enumerate(others):
            other.plot_ffspl(ax=ax, with_qn=0, color=mkcolor(count+1), show=False)

        return fig


class PspsReader(ETSF_Reader):
    """
    This object reads the results stored in the PSPS file produced by ABINIT.
    It provides helper function to access the most important quantities.
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

        # TODO
        #self.psps_files = []
        #for strng in r.read_value("filpsp"):
        #    s = "".join(strng)
        #    print(s)
        #    self.psps_files.append(s)
        #print(self.psps_files)

    def read_coresd(self, rmax=None):
        """
        Read the core charges and derivatives for the different types of atoms.

        Args:
            rmax: Maximum radius in Bohr. If None, data on the full grid is returned.

        Returns:
            meshes: List of ntypat arrays. Each array contains the linear meshes in real space.
            coresd: List with nytpat arrays of shape [6, npts].

            (np.zeros. np.zeros) if core charge is not present

        xccc1d[ntypat6,n1xccc*(1-usepaw)]

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
                # Truncate mesh
                ir_stop = min(int(rmax / step), npts) + 1
                #print(rmax, step, ir_stop, npts)

            rmeshes.append(rvals[:ir_stop])
            coresd.append(all_coresd[itypat, :, :ir_stop])

        return rmeshes, coresd

    def read_tcorespl(self):
        """
        Returns:
            qmesh: Linear q-mesh in G-space
            tcorespl:

        tcorespl[ntypat, 2, mqgrid_vl]
        Gives the pseudo core density in reciprocal space on a regular grid.
        Only if has_tcore
        """
        return self.read_value("qgrid_vl"), self.read_value("nc_tcorespl")

    def read_vlspl(self):
        """
        Returns:
            qmesh: Linear q-mesh in G-space
            vlspl:

        vlspl[2, ntypat, mqgrid_vl]
        Gives, on the radial grid, the local part of each type of psp.
        """
        return self.read_value("qgrid_vl"), self.read_value("vlspl")

    def read_projectors(self):
        """
        ffspl(ntypat, lnmax, 2, mqgrid_ff]
        Gives, on the radial grid, the different non-local projectors,
        in both the norm-conserving case, and the PAW case
        """
        # ekb(dimekb,ntypat*(1-usepaw))
        ekb = self.read_value("ekb")
        qgrid_ff = self.read_value("qgrid_ff")
        ffspl = self.read_value("ffspl")
        #print("qgrid", qgrid_ff.min(), qgrid_ff.max())

        projs = self.ntypat * [None]
        for itypat in range(self.ntypat):
            projs_type = []
            ln_list = self.get_lnlist_for_type(itypat)
            for i, ln in enumerate(ln_list):
                #print(ffspl[itypat, i, :, :])
                p = VnlProjector(itypat, ln, ekb[itypat, i], qgrid_ff, ffspl[itypat, i, :, :])
                projs_type.append(p)

            projs[itypat] = projs_type

        return projs

    def get_lnlist_for_type(self, itypat):
        """Return a list of (l, n) indices for this atom type."""
        # indlmn(6,lmn_size,ntypat)=array giving l,m,n,lm,ln,s for i=lmn
        indlmn_type = self.indlmn[itypat, :, :]

        iln0 = 0; ln_list = []
        for ilmn in range(self.lmnmax):
            iln = indlmn_type[ilmn, 4]
            if iln > iln0:
                iln0 = iln
                l = indlmn_type[ilmn, 0]  # l
                n = indlmn_type[ilmn, 2]  # n
                ln_list.append((l, n))

        return ln_list


class VnlProjector(object):
    """Data and parameters associated to a non-local projector."""
    def __init__(self, itypat, ln, ekb, qmesh, data):
        """
        Args:
            itypat:
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
        return self.data[0, :]

    @property
    def der2(self):
        """Second order derivative."""
        return self.data[1, :]

    @property
    def ecuts(self):
        """List of cutoff energies corresponding to self.qmesh."""
        return 2 * (np.pi * self.qmesh)**2

    @property
    def sign_sqrtekb(self):
        return np.sign(self.ekb) * np.sqrt(np.abs(self.ekb))
