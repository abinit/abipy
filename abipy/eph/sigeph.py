# coding: utf-8
"""
This module contains objects for the postprocessing of Sigma_eph calculations.

For a theoretical introduction see :cite:`Giustino2017`

Warning:

    Work in progress, DO NOT USE THIS CODE
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pandas as pd
import pymatgen.core.units as units
import abipy.core.abinit_units as abu

from collections import OrderedDict, namedtuple
from tabulate import tabulate
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import KpointList
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title)
from abipy.tools import duck
from abipy.electrons.ebands import ElectronsReader, RobotWithEbands
#from abipy.dfpt.phonons import PhononBands, RobotWithPhbands, factor_ev2units, unit_tag, dos_label_from_units
from abipy.abio.robots import Robot


# TODO QPState and QPList from electrons.gw (Define base abstract class?).
# __eq__ based on skb?
# broadening, adiabatic ??

class QpTempState(namedtuple("QpTempState", "tmesh e0 qpe ze0 spin kpoint band")):
    """
    Quasi-particle result for given (spin, kpoint, band).

    .. Attributes:

        spin: Spin index (C convention, i.e >= 0)
        kpoint: |Kpoint| object.
        band: band index. (C convention, i.e >= 0).
        tmesh: Temperature mesh in Kelvin.
        e0: Initial KS energy.
        qpe: Quasiparticle energy (complex) computed with the perturbative approach.
        ze0: Renormalization factor computed at e=e0.

    .. note::

        Energies are in eV.
    """

    @property
    def qpeme0(self):
        """E_QP[T] - E_0"""
        return self.qpe - self.e0

    @property
    def re_qpe(self):
        """Real part of the QP energy."""
        return self.qpe.real

    @property
    def imag_qpe(self):
        """Imaginay part of the QP energy."""
        return self.qpe.imag

    @property
    def skb(self):
        """Tuple with (spin, kpoint, band)"""
        return self.spin, self.kpoint, self.band

    #def copy(self):
    #    """Shallow copy."""
    #    d = {f: copy.copy(getattr(self, f)) for f in self._fields}
    #    return self.__class__(**d)

    @classmethod
    def get_fields(cls, exclude=()):
        fields = list(cls._fields) + ["qpeme0"]
        for e in exclude:
            fields.remove(e)
        return tuple(fields)

    def _repr_html_(self):
        """Integration with jupyter_ notebooks."""
        return self.get_dataframe()._repr_html_()

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0, title=None):
        """
        String representation with verbosity level ``verbose`` and optional ``title``.
        """
        s = str(self.get_dataframe())
        return "\n".join([marquee(title, mark="="), s]) if title is not None else s

    def get_dataframe(self, index=None, params=None):
        """
        Build pandas dataframe with QP data

        Args:
            index: dataframe index.
            params: Optional (Ordered) dictionary with extra parameters.

        Return: |pandas-DataFrame|
        """
        od = OrderedDict()
        for k in "tmesh e0 qpe qpeme0 ze0 spin kpoint band".split():
            if k in ("e0", "spin", "kpoint", "band"):
                od[k] = [getattr(self, k)] * len(self.tmesh)
            else:
                od[k] = getattr(self, k)

        if params is not None: od.update(params)

        return pd.DataFrame(od, index=index)

    @add_fig_kwargs
    def plot(self, with_fields="all", exclude_fields=None, ax_list=None, label=None, **kwargs):
        """
        Plot the QP results as function of temperature.

        Args:
            with_fields: The names of the QpTempState attributes to plot as function of e0.
                Accepts: List of strings or string with tokens separated by blanks.
                See :class:`QpTempState` for the list of available fields.
            exclude_fields: Similar to `with_field` but excludes fields.
            ax_list: List of matplotlib axes for plot. If None, new figure is produced.
            label: Label for plot.

        Returns: |matplotlib-Figure|
        """
        fields = _get_fields_for_plot("temp", with_fields, exclude_fields)
        if not fields: return None

        num_plots, ncols, nrows = len(fields), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        # Build plot grid.
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        linestyle = kwargs.pop("linestyle", "o")
        #kw_color = kwargs.pop("color", None)
        #kw_label = kwargs.pop("label", None)
        for i, (field, ax) in enumerate(zip(fields, ax_list)):
            irow, icol = divmod(i, ncols)
            ax.grid(True)
            if irow == nrows - 1: ax.set_xlabel("Temperature [K]")
            ax.set_ylabel(field)
            yy = getattr(self, field)
            lbl = label if i == 0 and label is not None else None

            # Handle complex arrays
            #if np.iscomplexobj(yy):
            #    ax.plot(self.tmesh, yy.real, linestyle, label=lbl, **kwargs)
            #    ax.plot(self.tmesh, yy.imag, linestyle, label=lbl, **kwargs)
            #else:
            ax.plot(self.tmesh, yy.real, linestyle, label=lbl, **kwargs)

        # Get around a bug in matplotlib
        if num_plots % ncols != 0: ax_list[-1].axis('off')

        if lbl is not None:
            ax_list[0].legend(loc="best")

        fig.tight_layout()
        return fig


def _get_fields_for_plot(vs, with_fields, exclude_fields):
    """
    Return list of QpTempState fields to plot from input arguments.
    vs in ["temp", "e0"]
    """
    if vs == "temp":
        all_fields = list(QpTempState.get_fields(exclude=["spin", "kpoint", "band", "e0", "tmesh"]))
    elif vs == "e0":
        all_fields = list(QpTempState.get_fields(exclude=["spin", "kpoint", "band", "e0", "tmesh"]))
    else:
        raise ValueError("Invalid vs: `%s`" % str(vs))

    # Initialize fields.
    if duck.is_string(with_fields) and with_fields == "all":
        fields = all_fields
    else:
        fields = list_strings(with_fields)
        for f in fields:
            if f not in all_fields:
                raise ValueError("Field %s not in allowed values %s" % (f, str(all_fields)))

    # Remove entries
    if exclude_fields:
        if duck.is_string(exclude_fields):
            exclude_fields = exclude_fields.split()
        for e in exclude_fields:
            fields.remove(e)

    return fields


class QpTempList(list):
    """
    A list of quasiparticle corrections (usually for a given spin).
    """
    def __init__(self, *args, **kwargs):
        super(QpTempList, self).__init__(*args)
        self.is_e0sorted = kwargs.get("is_e0sorted", False)

    @property
    def tmesh(self):
        """Temperature mesh in Kelvin."""
        if len(self):
            return self[0].tmesh
        return []

    @property
    def ntemp(self):
        """Number of temperatures."""
        return len(self.tmesh)

    def __repr__(self):
        return "<%s at %s, len=%d>" % (self.__class__.__name__, id(self), len(self))

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0, title=None):
        """String representation."""
        return " "

    #def copy(self):
    #    """Copy of self."""
    #    return self.__class__([qp.copy() for qp in self], is_e0sorted=self.is_e0sorted)

    def sort_by_e0(self):
        """Return a new :class:`QpTempList` object with the E0 energies sorted in ascending order."""
        return self.__class__(sorted(self, key=lambda qp: qp.e0), is_e0sorted=True)

    def get_e0mesh(self):
        """Return the E0 energies."""
        if not self.is_e0sorted:
            raise ValueError("QPState corrections are not sorted. Use sort_by_e0")

        return np.array([qp.e0 for qp in self])

    def get_field_itemp(self, field, itemp):
        """|numpy-array| containing the values of field at temperature ``itemp``"""
        return np.array([getattr(qp, field)[itemp] for qp in self])

    #def get_skb_field(self, skb, field):
    #    """Return the value of field for the given spin kp band tuple, None if not found"""
    #    for qp in self:
    #        if qp.skb == skb:
    #            return getattr(qp, field)
    #    return None

    #def get_qpenes(self):
    #    """Return an array with the :class:`QPState` energies."""
    #    return self.get_field("qpe")

    #def get_qpeme0(self):
    #    """Return an arrays with the :class:`QPState` corrections."""
    #    return self.get_field("qpeme0")

    def merge(self, other, copy=False):
        """
        Merge self with other. Return new :class:`QpTempList` object

        Raise: `ValueError` if merge cannot be done.
        """
        skb0_list = [qp.skb for qp in self]
        for qp in other:
            if qp.skb in skb0_list:
                raise ValueError("Found duplicated (s,b,k) indexes: %s" % str(qp.skb))

        # Test on tmesh?
        # TODO Why copy?
        qps = self.copy() + other.copy() if copy else self + other
        return self.__class__(qps)

    @add_fig_kwargs
    def plot_vs_e0(self, itemp_list="all", with_fields="all", exclude_fields=None, fermie=None, colormap="jet",
                   ax_list=None, xlims=None, fontsize=12, **kwargs):
        """
        Plot the QP results as a function of the initial KS energy.

        Args:
            itemp_list: "all" to plot all temperatures. List of integers to select a particular temperature.
            with_fields: The names of the qp attributes to plot as function of e0.
                Accepts: List of strings or string with tokens separated by blanks.
                See :class:`QPState` for the list of available fields.
            exclude_fields: Similar to `with_field` but excludes fields.
            fermie: Value of the Fermi level used in plot. None for absolute e0s.
            colormap: matplotlib color map.
            ax_list: List of |matplotlib-Axes| for plot. If None, new figure is produced.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            fontsize: Legend and title fontsize.
            kwargs: linestyle, color, label, marker

        Returns: |matplotlib-Figure|
        """
        fields = _get_fields_for_plot("e0", with_fields, exclude_fields)
        if not fields: return None

        num_plots, ncols, nrows = len(fields), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        # Build plot grid.
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()
        cmap = plt.get_cmap(colormap)

        # Get QpTempList and sort it.
        qps = self if self.is_e0sorted else self.sort_by_e0()
        e0mesh = qps.get_e0mesh()
        xlabel = r"$\epsilon_{KS}\;[eV]$"
        if fermie is not None:
            e0mesh -= fermie
            xlabel = r"$\epsilon_{KS}-\epsilon_F\;[eV]$"

        kw_linestyle = kwargs.pop("linestyle", "o")
        kw_color = kwargs.pop("color", None)
        kw_label = kwargs.pop("label", None)
        if duck.is_string(itemp_list) and itemp_list == "all":
            itemp_list = range(self.ntemp)

        for i, (field, ax) in enumerate(zip(fields, ax_list)):
            irow, icol = divmod(i, ncols)
            ax.grid(True)
            if irow == nrows - 1:
                ax.set_xlabel(xlabel)
            ax.set_ylabel(field, fontsize=fontsize)
            for itemp in itemp_list:
                yt = qps.get_field_itemp(field, itemp)
                # TODO real and imag?
                label = kw_label
                if kw_label is None:
                    label = "T = %.1f K" % self.tmesh[itemp] if i == 0 else None
                ax.plot(e0mesh, yt.real, kw_linestyle,
                        color=cmap(itemp / self.ntemp) if kw_color is None else kw_color,
                        label=label, **kwargs)
            set_axlims(ax, xlims, "x")

        # Get around a bug in matplotlib
        if num_plots % ncols != 0: ax_list[-1].axis('off')

        if kw_label:
            ax_list[0].legend(loc="best", fontsize=fontsize, shadow=True)

        return fig


class EphSelfEnergy(object):
    r"""
    Electron self-energy due to phonon interaction :math:`\Sigma_{nk}(\omega,T)`
    Actually these are the diagonal matrix elements in the KS basis set.
    """
    # Symbols used in matplotlib plots.
    latex_symbol = dict(
        re=r"$\Re{\Sigma(\omega)}$",
        im=r"$\Im{\Sigma(\omega)}$",
        spfunc=r"$A(\omega)}$",
    )

    def __init__(self, wmesh, qp, vals_e0ks, dvals_de0ks, dw_vals, vals_wr, spfunc_wr):
        """
        Args:
            wmesh: Frequency mesh in eV.
            qp: :class:`QpTempState` instance.
            vals_e0ks: complex [ntemp] array with Sigma_eph(omega=eKS, kT)
            dvals_de0ks: complex [ntemp] arrays with d Sigma_eph(omega, kT) / d omega (omega=eKS)
            dw_vals: [ntemp] array with Debye-Waller term (static)
            vals_wr: [ntemp, nwr] complex array with Sigma_eph(omega, kT). enk_KS corresponds to nwr//2 + 1.
            spfunc_wr: [ntemp, nwr] real array with spectral function.
        """
        self.qp = qp
        self.spin, self.kpoint, self.band = qp.spin, qp.kpoint, qp.band
        self.wmesh, self.tmesh = wmesh, qp.tmesh
        self.nwr, self.ntemp = len(self.wmesh), len(self.tmesh)

        # Get refs to values
        self.vals_e0ks = vals_e0ks
        assert self.vals_e0ks.shape == (self.ntemp,)
        self.dvals_de0ks = dvals_de0ks
        assert self.dvals_de0ks.shape == (self.ntemp,)
        self.dw_vals = dw_vals
        assert self.dw_vals.shape == (self.ntemp,)
        self.vals_wr = vals_wr
        assert self.vals_wr.shape == (self.ntemp, self.nwr)
        self.spfunc_wr = spfunc_wr
        assert self.spfunc_wr.shape == (self.ntemp, self.nwr)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0, title=None):
        """String representation."""
        lines = []; app = lines.append
        if title is not None: app(marquee(title, mark="="))
        app("K-point: %s, band: %d, spin: %d" % (repr(self.kpoint), self.band, self.spin))
        app("Number of temperatures: %d, from %.1f to %.1f [K]" % (self.ntemp, self.tmesh[0], self.tmesh[-1]))
        app("Number of frequencies: %d, from %.1f to %.1f [eV]" % (self.nwr, self.wmesh[0], self.wmesh[-1]))
        app(self.qp.to_string(verbose=verbose, title="QP data"))

        return "\n".join(lines)

    def _get_wmesh_xlabel(self, zero_energy):
        """Return (wmesh, xlabel) from zero_energy input argument."""
        if zero_energy is None:
            xx = self.wmesh
            xlabel = r"$\omega\;[eV]$"
        elif zero_energy == "e0":
            xx = self.wmesh - self.qp.e0
            xlabel = r"$\omega - \epsilon_{nk}\;[eV]$"
        # TODO: chemical potential? but then I have mu(T) to handle in plots!
        #elif zero_energy == "fermie":
        #    xx = self.wmesh - self.fermie
        #    xlabel = r"$\omega\;[eV]$"
        else:
            raise ValueError("Invalid value of zero_energy: `%s`" % str(zero_energy))

        return xx.copy(), xlabel

    def _get_ys_itemp(self, what, itemp):
        """Return array(T) to plot from what and itemp index."""
        return dict(
            re=self.vals_wr[itemp].real,
            im=self.vals_wr[itemp].imag,
            spfunc=self.spfunc_wr[itemp],
        )[what]

    def _get_itemps_labels(self, itemps):
        """Return list of temperature indices and labels from itemps."""
        if duck.is_string(itemps):
            if itemps == "all":
                itemps = list(range(self.ntemp))
            else:
                raise ValueError("Invalid value for itemps: `%s`" % str(itemps))
        else:
            itemps = np.array(itemps, dtype=np.int)
            itemps = [itemps] if itemps.size == 1 else itemps.tolist()
            if not all(self.ntemp > it >= 0 for it in itemps):
                raise ValueError("Invalid list of temperature indices. ntemp is %d, received itemps:\n\t%s" % (
                        self.ntemp, str(itemps)))

        return itemps, ["T=%.1f K" % self.tmesh[it] for it in itemps]

    @add_fig_kwargs
    def plot_tdep(self, itemps="all", zero_energy="e0", colormap="jet", ax_list=None,
                  what_list=("re", "im", "spfunc"), xlims=None, fontsize=8, **kwargs):
        """
        Plot the real/imaginary part of self-energy as well as the spectral function for
        the different temperatures with a color map.

        Args:
            itemps: List of temperature indices. "all" to plot'em all.
            zero_energy:
            colormap: matplotlib color map.
            ax_list: List of |matplotlib-Axes|. If None, new figure is produced.
            what_list: List of strings selecting the quantity to plot.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                or scalar e.g. ``left``. If left (right) is None, default values are used.
            fontsize: legend and label fontsize.
            kwargs: Keyword arguments passed to ax.plot

        Returns: |matplotlib-Figure|
        """
        # FIXME zero_energy or e0?
        what_list = list_strings(what_list)
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(what_list), ncols=1, sharex=True, sharey=False)
        cmap = plt.get_cmap(colormap)
        xs, xlabel = self._get_wmesh_xlabel(zero_energy)

        itemps, tlabels = self._get_itemps_labels(itemps)
        kw_color = kwargs.get("color", None)
        kw_label = kwargs.get("label", None)
        for i, (what, ax) in enumerate(zip(what_list, ax_list)):
            ax.grid(True)
            ax.set_ylabel(self.latex_symbol[what])
            if (i == len(ax_list) - 1): ax.set_xlabel(xlabel)
            for itemp in itemps:
                ax.plot(xs, self._get_ys_itemp(what, itemp),
                        color=cmap(itemp / self.ntemp) if kw_color is None else kw_color,
                        label=tlabels[itemp] if (i == 0 and kw_label is None) else kw_label,
                )
            if i == 0: ax.legend(loc="best", shadow=True, fontsize=fontsize)
            set_axlims(ax, xlims, "x")

        if "title" not in kwargs:
            title = "K-point: %s, band: %d, spin: %d" % (repr(self.kpoint), self.band, self.spin)
            fig.suptitle(title, fontsize=fontsize)

        return fig


class SigEPhFile(AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter):
    """
    This file contains the Fan-Migdal Debye-Waller self-energy, the |ElectronBands| on the k-mesh.
    Provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        with SigEPhFile("out_SIGEPH.nc") as ncfile:
            print(ncfile)
            ncfile.ebands.plot()

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: SigEPhFile
    """
    # Markers used for up/down bands.
    marker_spin = {0: "^", 1: "v"}
    color_spin = {0: "k", 1: "r"}

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath):
        super(SigEPhFile, self).__init__(filepath)
        self.reader = r = SigmaPhReader(filepath)

        # Get important dimensions.
        self.nkcalc = r.nkcalc
        self.ntemp = r.ntemp
        self.nqbz = r.nqbz
        self.nqibz = r.nqibz
        self.ngqpt = r.ngqpt

        self.symsigma = r.read_value("symsigma")
        # TODO zcut?
        self.zcut = r.read_value("eta")
        self.nbsum = int(r.read_value("nbsum"))

        self.bstart_sk = self.reader.bstart_sk
        self.nbcalc_sk = self.reader.nbcalc_sk
        self.bstop_sk = self.reader.bstop_sk

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
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="KS Electronic Bands"))
        app("")
        # SigmaEPh section.
        app(marquee("SigmaEPh calculation", mark="="))
        app("Number of k-points computed: %d" % (self.nkcalc))
        app("Q-mesh: nqibz: %s, nqbz: %s, ngqpt: %s" % (self.nqibz, self.nqbz, str(self.ngqpt)))
        app("Number of bands included in self-energy: %d" % (self.nbsum))
        app("zcut: %.3f [Ha], %.3f [eV]" % (self.zcut, self.zcut * units.Ha_to_eV))
        app("Number of temperatures: %d, from %.1f to %.1f [K]" % (self.ntemp, self.tmesh[0], self.tmesh[-1]))
        app("symsigma: %s" % (self.symsigma))
        app("Has Spectral function: %s" % (self.has_spectral_function))

        # Build Table with direct gaps. Only the results for the first and the last T are shown if not verbose.
        if verbose:
            it_list = list(range(self.ntemp))
        else:
            it_list = [0, -1] if self.ntemp != 1 else [0]

        for it in it_list:
            app("\nKS and QP direct gaps for T = %.1f K:" % self.tmesh[it])
            data = []
            for ik, kpoint in enumerate(self.sigma_kpoints):
                for spin in range(self.nsppol):
                    ks_gap = self.ks_dirgaps[spin, ik]
                    qp_gap = self.qp_dirgaps_t[spin, ik, it]
                    data.append([spin, repr(kpoint), ks_gap, qp_gap, qp_gap - ks_gap])
            app(str(tabulate(data, headers=["Spin", "K-point", "KS_gap", "QP_gap", "QP - KS"], floatfmt=".3f")))
            app("")

        if verbose > 1:
            app("K-points and bands included in self-energy corrections:")
            for spin in range(self.nsppol):
                for ik, kpoint in enumerate(self.sigma_kpoints):
                    post = "ik: %d" % (ik if self.nsppol == 1 else "ik: %d, spin: %d" % (ik, spin))
                    app("\t%s: bstart: %d, bstop: %d, %s" % (
                        repr(kpoint), self.bstart_sk[spin, ik], self.bstop_sk[spin, ik], post))

        return "\n".join(lines)

    @lazy_property
    def ebands(self):
        """|ElectronBands| object."""
        return self.reader.read_ebands()

    @property
    def structure(self):
        """|Structure| object."""
        return self.ebands.structure

    def close(self):
        """Close the file."""
        self.reader.close()

    @property
    def has_spectral_function(self):
        """True if file contains spectral function data."""
        return self.reader.nwr != 0

    @property
    def sigma_kpoints(self):
        """The K-points where QP corrections have been calculated."""
        return self.reader.sigma_kpoints

    @property
    def tmesh(self):
        """Temperature mesh in Kelvin."""
        return self.reader.tmesh

    @lazy_property
    def ks_dirgaps(self):
        """
        |numpy-array| of shape [nsppol, nkcalc] with the KS gaps in eV ordered as kcalc.
        """
        return self.reader.read_value("ks_gaps") * units.Ha_to_eV

    @lazy_property
    def qp_dirgaps_t(self):
        """
        |numpy-array| of shape [nsppol, nkcalc, ntemp] with the QP gaps in eV ordered as kcalc.
        """
        return self.reader.read_value("qp_gaps") * units.Ha_to_eV

    @lazy_property
    def mu_e(self):
        """mu_e[ntemp] chemical potential (eV) of electrons for the different temperatures."""
        return self.reader.read_value("mu_e") * units.Ha_to_eV

    # TODO
    #integer,allocatable :: kcalc2ibz(:,:)
    #!kcalc2ibz (nkcalc, 6))
    #! Mapping kcalc --> ibz as reported by listkk.

    def sigkpt2index(self, kpoint):
        """
        Returns the index of the self-energy k-point in sigma_kpoints
        Used to access data in the arrays that are dimensioned with [0:nkcalc]
        """
        return self.reader.sigkpt2index(kpoint)

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with the convergence parameters, e.g. ``nbsum``."""
        return OrderedDict([
            ("nbsum", self.nbsum),
            ("zcut", self.zcut),
            ("symsigma", self.symsigma),
            ("nqbz", self.reader.nqbz),
            ("nqibz", self.reader.nqibz),
        ])

    def get_dataframe(self, with_params=True, ignore_imag=False):
        """
        Returns |pandas-Dataframe| with QP results for all k-points, bands and spins
        included in the calculation.

        Args:
            with_params: False to exclude calculation parameters from the dataframe.
            ignore_imag: only real part is returned if ``ignore_imag``.
        """
        df_list = []; app = df_list.append
        for spin in range(self.nsppol):
            for ik, kpoint in enumerate(self.sigma_kpoints):
                app(self.get_dataframe_sk(spin, ik, with_params=with_params, ignore_imag=ignore_imag))

        return pd.concat(df_list)

    def get_dataframe_sk(self, spin, kpoint, itemp=None, index=None, with_params=True, ignore_imag=False):
        """
        Returns |pandas-DataFrame| with QP results for the given (spin, k-point).

        Args:
            spin: Spin index.
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            itemp: Temperature index, if None all temperatures are returned.
            index: dataframe index.
            with_params: False to exclude calculation parameters from the dataframe.
            ignore_imag: Only real part is returned if ``ignore_imag``.
        """
        ik = self.sigkpt2index(kpoint)
        rows = []
        for band in range(self.bstart_sk[spin, ik], self.bstop_sk[spin, ik]):
            # Read QP data.
            qp = self.reader.read_qp(spin, ik, band, ignore_imag=ignore_imag)
            # Convert to dataframe and add other entries useful when comparing different calculations.
            rows.append(qp.get_dataframe(params=self.params if with_params else None))

        df = pd.concat(rows)
        if itemp is not None: df = df[df["tmesh"] == self.tmesh[itemp]]
        return df

    #def get_dirgaps_dataframe(self):

    @add_fig_kwargs
    def plot_qpgaps_t(self, ax_list=None, plot_qpmks=False, fontsize=8, **kwargs):
        """
        Plot the KS and the QP(T) direct gaps for all the k-points available on file.

        Args:
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            plot_qpmks: If True, plot QP_gap - KS_gap
            fontsize: legend and title fontsize.
            kwargs: Passed to ax.plot method except for marker.

        Returns: |matplotlib-Figure|
        """
        # Build grid plot.
        nrows, ncols = len(self.sigma_kpoints), 1
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=True)
        ax_list = ax_list.ravel()
        label = kwargs.pop("label", None)

        for ik, (kpt, ax) in enumerate(zip(self.sigma_kpoints, ax_list)):
            for spin in range(self.nsppol):
                if not plot_qpmks:
                    # Plot QP_{spin,kpt}(T)
                    ax.plot(self.tmesh, self.qp_dirgaps_t[spin, ik], marker=self.marker_spin[spin],
                            label=label, **kwargs) #label="QP gap")
                    # Add KS gap (assumed at T=0).
                    ax.scatter(0, self.ks_dirgaps[spin, ik]) #, label="KS gap %s" % label)
                else:
                    # Plot QP_{spin,kpt}(T) - KS_gap
                    ax.plot(self.tmesh, self.qp_dirgaps_t[spin, ik] - self.ks_dirgaps[spin, ik],
                            marker=self.marker_spin[spin]) #, label="QP-KS gap %s")

            ax.grid(True)
            if ik == len(self.sigma_kpoints) - 1:
                ax.set_xlabel("Temperature [K]")
            if plot_qpmks:
                ax.set_ylabel("QP-KS gap [eV]")
            else:
                ax.set_ylabel("QP direct gap [eV]")
            ax.set_title("k:%s" % (repr(kpt)), fontsize=fontsize)
            if label:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_qpdata_t(self, spin, kpoint, band_list=None, fontsize=12, **kwargs):
        """
        Plot the QP results as function T for a given (spin, k-point) and all bands.

        Args:
            spin: Spin index
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            band_list: List of band indices to be included. If None, all bands are shown.
            fontsize: legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        # TODO: Add more quantities DW, Fan(0)
        # Quantities to plot.
        what_list = ["re_qpe", "imag_qpe", "ze0"]

        # Build grid plot.
        nrows, ncols = len(what_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # Read all QPs for this (spin, kpoint) and all bands.
        qp_list = self.reader.read_qplist_sk(spin, kpoint)

        for i, (ax, what) in enumerate(zip(ax_list, what_list)):
            # Plot QP(T)
            for qp in qp_list:
                if band_list is not None and qp.band not in band_list: continue
                yvals = getattr(qp, what)
                ax.plot(qp.tmesh, yvals, marker=self.marker_spin[spin],
                        label="band: %s" % qp.band)

            ax.grid(True)
            ax.set_ylabel(what)
            if i == len(what_list) - 1: ax.set_xlabel("Temperature [K]")
            if i == 0: ax.legend(loc="best", fontsize=fontsize, shadow=True)

        if "title" not in kwargs:
            title = "QP results spin:%s, k:%s" % (spin, repr(qp_list[0].kpoint))
            fig.suptitle(title, fontsize=fontsize)

        return fig

    @lazy_property
    def qplist_spin(self):
        """Tuple of :class:`QpTempList` objects indexed by spin."""
        return self.reader.read_allqps()

    @add_fig_kwargs
    def plot_qps_vs_e0(self, itemp_list="all", with_fields="all", exclude_fields=None, e0="fermie",
                       colormap="jet", xlims=None, ax_list=None, fontsize=8, **kwargs):
        """
        Plot the QP results in the SIGEPH file as function of the initial KS energy.

        Args:
            itemp_list: "all" to plot all temperatures. List of integers to select a particular temperature.
            with_fields: The names of the qp attributes to plot as function of e0.
                Accepts: List of strings or string with tokens separated by blanks.
                See :class:`QPState` for the list of available fields.
            exclude_fields: Similar to ``with_field`` but excludes fields.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            ax_list: List of |matplotlib-Axes| for plot. If None, new figure is produced.
            colormap: matplotlib color map.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        fermie = self.ebands.get_e0(e0)
        for spin in range(self.nsppol):
            fig = self.qplist_spin[spin].plot_vs_e0(itemp_list=itemp_list,
                with_fields=with_fields, exclude_fields=exclude_fields, fermie=fermie,
                colormap=colormap, xlims=xlims, ax_list=ax_list, fontsize=fontsize, marker=self.marker_spin[spin],
                show=False, **kwargs)
            ax_list = fig.axes

        return fig

    def write_notebook(self, nbpath=None, title=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=title)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
            nbv.new_code_cell("ncfile.ebands.plot();"),
            #nbv.new_code_cell("ncfile.get_dirgaps_dataframe()"),
            #nbv.new_code_cell("alldata = ncfile.get_dataframe()\nalldata"),
            nbv.new_code_cell("ncfile.plot_qpgaps_t();"),
            nbv.new_code_cell("#ncfile.plot_qpgaps_t(plot_qpmks=True);"),
            nbv.new_code_cell("ncfile.plot_qps_vs_e0();"),
            nbv.new_code_cell("""\
for spin in range(ncfile.nsppol):
    for sigma_kpoint in ncfile.sigma_kpoints:
        ncfile.plot_qpdata_t(spin, sigma_kpoint, band_list=None, fontsize=12);"""),
            nbv.new_code_cell("#df = ncfile.get_dataframe_sk(spin=0, kpoint=(0, 0, 0))"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class SigEPhRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple SIGEPH.nc files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: SigEPhRobot
    """
    # Try to have API similar to SigresRobot
    EXT = "SIGEPH"

    def __init__(self, *args):
        super(SigEPhRobot, self).__init__(*args)
        if len(self.abifiles) in (0, 1): return

        # Check dimensions and self-energy states and issue warning.
        warns = []; wapp = warns.append
        nc0 = self.abifiles[0]
        same_nsppol, same_nkcalc = True, True
        if any(nc.nsppol != nc0.nsppol for nc in self.abifiles):
            same_nsppol = False
            wapp("Comparing ncfiles with different values of nsppol.")
        if any(nc.nkcalc != nc0.nkcalc for nc in self.abifiles):
            same_nkcalc = False
            wapp("Comparing ncfiles with different number of k-points in self-energy. Doh!")

        if same_nsppol and same_nkcalc:
            # FIXME
            # Different values of bstart_ks are difficult to handle
            # Because the high-level API assumes an absolute global index
            # Should decide how to treat this case: either raise or interpret band as an absolute band index.
            if any(np.any(nc.bstart_sk != nc0.bstart_sk) for nc in self.abifiles):
                wapp("Comparing ncfiles with different values of bstart_sk")
            if any(np.any(nc.bstop_sk != nc0.bstop_sk) for nc in self.abifiles):
                wapp("Comparing ncfiles with different values of bstop_sk")

        if warns:
            for w in warns:
                cprint(w, color="yellow")

    def _check_dims_and_params(self):
        """Test that nsppol, sigma_kpoints, tlist are consistent."""
        if not len(self.abifiles) > 1:
            return

        nc0 = self.abifiles[0]
        errors = []
        eapp = errors.append

        if any(nc.nsppol != nc0.nsppol for nc in self.abifiles[1:]):
            eapp("Files with different values of `nsppol`")

        if any(nc.nkcalc != nc0.nkcalc for nc in self.abifiles[1:]):
            eapp("Files with different values of `nkcalc`")
        else:
            for nc in self.abifiles[1:]:
                for k0, k1 in zip(nc0.sigma_kpoints, nc.sigma_kpoints):
                    if k0 != k1:
                        eapp("Files with different values of `sigma_kpoints`")

        if any(not np.allclose(nc.tmesh, nc0.tmesh) for nc in self.abifiles[1:]):
            eapp("Files with different tmesh")

        if errors:
            raise ValueError("Cannot compare multiple SIGEPH.nc files. Reason:\n %s" % "\n".join(errors))

    def get_dataframe_sk(self, spin, kpoint, with_params=True, ignore_imag=False):
        """
        Return |pandas-Dataframe| with qp results for this spin, k-point

        Args:
            spin: Spin index
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            with_params:
            ignore_imag: only real part is returned if ``ignore_imag``.
        """
        df_list = []; app = df_list.append
        for label, ncfile in self.items():
            df = ncfile.get_dataframe_sk(spin, kpoint, index=None,
                                         with_params=with_params, ignore_imag=ignore_imag)
        return pd.concat(df_list)

    def get_dataframe(self, with_params=True, ignore_imag=False):
        """
        Return |pandas-Dataframe| with QP results for all k-points, bands and spins
        present in the files treated by the robot.

        Args:
            with_params:
            ignore_imag: only real part is returned if ``ignore_imag``.
        """
        df_list = []; app = df_list.append
        for label, ncfile in self.items():
            for spin in range(ncfile.nsppol):
                for ik, kpoint in enumerate(ncfile.sigma_kpoints):
                    app(ncfile.get_dataframe_sk(spin, ik, with_params=with_params, ignore_imag=ignore_imag))
        return pd.concat(df_list)

    @add_fig_kwargs
    def plot_selfenergy_conv(self, spin, kpoint, band, itemp=0, sortby=None, hue=None,
                             colormap="jet", xlims=None, fontsize=8, **kwargs):
        """
        Plot the convergence of the EPH self-energy wrt to the ``sortby`` parameter.
        Values can be optionally grouped by `hue`.

        Args:
            spin: Spin index.
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            band: Band index.
            itemp: Temperature index.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            colormap: matplotlib color map.
            xlims: Set the data limits for the x-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        # Make sure that nsppol, sigma_kpoints, and tmesh are consistent.
        self._check_dims_and_params()
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(colormap)

        if hue is None:
            ax_list = None
            lnp_list = self.sortby(sortby)
            for i, (label, ncfile, param) in enumerate(lnp_list):
                sigma = ncfile.reader.read_sigeph_skb(spin, kpoint, band)
                fig = sigma.plot_tdep(itemps=itemp, ax_list=ax_list,
                    label=label, color=cmap(i/len(lnp_list)), show=False)
                ax_list = fig.axes
        else:
            # group_and_sortby and build (3, ngroups) subplots
            groups = self.group_and_sortby(hue, sortby)
            nrows, ncols = 3, len(groups)
            ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                    sharex=True, sharey=True, squeeze=False)
            for ig, g in enumerate(groups):
                subtitle = "%s: %s" % (self._get_label(hue), g.hvalue)
                ax_mat[0, ig].set_title(subtitle, fontsize=fontsize)
                for i, (nclabel, ncfile, param) in enumerate(g):
                    sigma = ncfile.reader.read_sigeph_skb(spin, kpoint, band)
                    fig = sigma.plot_tdep(itemps=itemp, ax_list=ax_mat[:, ig],
                        label="%s: %s" % (self._get_label(sortby), param),
                        color=cmap(i/len(g)), show=False)

            if ig != 0:
                for ax in ax_mat[:, ig]:
                    set_visible(ax, False, "ylabel")

            for ax in ax_mat.ravel():
                set_axlims(ax, xlims, "x")

        return fig

    @add_fig_kwargs
    def plot_qpgaps_t(self, plot_qpmks=False, sortby=None, hue=None, fontsize=8, **kwargs):
        """
        Compare the QP(T) direct gaps for all the k-points available in the robot.

        Args:
            plot_qpmks: If True, plot QP_gap - KS_gap
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        if hue is None:
            ax_list = None
            for i, ((label, ncfile, param), lineopt) in enumerate(zip(self.sortby(sortby), self.iter_lineopt())):
                fig = ncfile.plot_qpgaps_t(ax_list=ax_list, plot_qpmks=plot_qpmks,
                    label=label if i == 0 else None, show=False, fontsize=fontsize, **lineopt)
                ax_list = fig.axes
        else:
            # (nkcalc, ngroups) subplots
            groups = self.group_and_sortby(hue, sortby)
            nrows, ncols = self.abifiles[0].nkcalc, len(groups)
            ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                   sharex=True, sharey=False, squeeze=False)
            for ig, g in enumerate(groups):
                for i, (nclabel, ncfile, param) in enumerate(g):
                    fig = ncfile.plot_qpgaps_t(ax_list=ax_mat[:, ig], plot_qpmks=plot_qpmks,
                            label="%s: %s" % (self._get_label(sortby), param),
                            fontsize=fontsize, show=False) #, **lineopt)

                # Add label with hue to previous title with k-point info.
                ax_append_title(ax_mat[0, ig], "\n%s: %s" % (self._get_label(hue), g.hvalue), fontsize=fontsize)
                if ig != 0:
                    for ax in ax_mat[:, ig]:
                        set_visible(ax, False, "ylabel")

            # Hide legends except first one
            if nrows > 1:
                for ax in ax_mat[1:, :].ravel():
                    set_visible(ax, False, "legend")

        return fig

    @add_fig_kwargs
    def plot_qpgaps_convergence(self, itemp=0, sortby=None, hue=None, plot_qpmks=False, fontsize=8, **kwargs):
        """
        Plot the convergence of the direct QP gaps at given temperature
        for all the k-points treated by the robot.

        Args:
            itemp: Temperature index.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            plot_qpmks: If True, plot QP_gap - KS_gap
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        # Make sure that nsppol, sigma_kpoints, and tmesh are consistent.
        self._check_dims_and_params()

        nc0 = self.abifiles[0]
        nsppol, sigma_kpoints = nc0.nsppol, nc0.sigma_kpoints
        # Build grid with (nkpt, 1) plots.
        ncols, nrows = 1, len(sigma_kpoints)
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        if hue is None:
            labels, ncfiles, params = self.sortby(sortby, unpack=True)
        else:
            groups = self.group_and_sortby(hue, sortby)

        name = "QP dirgap" if not plot_qpmks else "QP-KS dirgap"
        for ik, (kpt, ax) in enumerate(zip(sigma_kpoints, ax_list)):
            for spin in range(nsppol):
                ax.set_title("%s k:%s, T = %.1f K" % (
                    name, repr(kpt), nc0.tmesh[itemp]), fontsize=fontsize)

                # Extract QP dirgap for [spin, kpt, itemp]
                if hue is None:
                    yvals = [ncfile.qp_dirgaps_t[spin, ik, itemp] for ncfile in ncfiles]
                    if plot_qpmks:
                        yvals = np.array(yvals) - np.array([ncfile.ks_dirgaps[spin, ik] for ncfile in ncfiles])
                    ax.plot(params, yvals, marker=nc0.marker_spin[spin])
                else:
                    for g in groups:
                        yvals = [ncfile.qp_dirgaps_t[spin, ik, itemp] for ncfile in g.abifiles]
                        if plot_qpmks:
                            yvals = np.array(yvals) - np.array([ncfile.ks_dirgaps[spin, ik] for ncfile in g.abifiles])
                        label = "%s: %s" % (self._get_label(hue), g.hvalue)
                        ax.plot(g.xvalues, yvals, marker=nc0.marker_spin[spin], label=label)

            ax.grid(True)
            ax.set_ylabel("%s [eV]" % name)
            if ik == len(sigma_kpoints) - 1:
                ax.set_xlabel("%s" % self._get_label(sortby))
                if sortby is None: rotate_ticklabels(ax, 15)
            if hue is not None:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_qpdata_conv_skb(self, spin, kpoint, band,
                             itemp=0, sortby=None, hue=None, fontsize=8, **kwargs):
        """
        Plot the convergence of the QP results at the given temperature for given (spin, kpoint, band)

        Args:
            spin: Spin index.
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            band: Band index.
            itemp: Temperature index.
            sortby: Define the convergence parameter, sort files and produce plot labels.
                Can be None, string or function. If None, no sorting is performed.
                If string and not empty it's assumed that the abifile has an attribute
                with the same name and `getattr` is invoked.
                If callable, the output of sortby(abifile) is used.
            hue: Variable that define subsets of the data, which will be drawn on separate lines.
                Accepts callable or string
                If string, it's assumed that the abifile has an attribute with the same name and getattr is invoked.
                If callable, the output of hue(abifile) is used.
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        # Make sure that nsppol, sigma_kpoints, and tmesh are consistent.
        self._check_dims_and_params()

        # TODO: Add more quantities DW, Fan(0)
        # TODO: Decide how to treat complex quantities, avoid annoying ComplexWarning
        # TODO: Format for g.hvalue
        # Treat fundamental gaps
        # Quantities to plot.
        what_list = ["re_qpe", "imag_qpe", "ze0"]

        # Build grid plot.
        nrows, ncols = len(what_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        nc0 = self.abifiles[0]
        ik = nc0.sigkpt2index(kpoint)
        kpoint = nc0.sigma_kpoints[ik]

        # Sort and read QP data.
        if hue is None:
            labels, ncfiles, params = self.sortby(sortby, unpack=True)
            qplist = [ncfile.reader.read_qp(spin, kpoint, band) for ncfile in ncfiles]
        else:
            groups = self.group_and_sortby(hue, sortby)
            qplist_group = []
            for g in groups:
                lst = [ncfile.reader.read_qp(spin, kpoint, band) for ncfile in g.abifiles]
                qplist_group.append(lst)

        for i, (ax, what) in enumerate(zip(ax_list, what_list)):
            if hue is None:
                # Extract QP data.
                yvals = [getattr(qp, what)[itemp] for qp in qplist]
                ax.plot(params, yvals, marker=nc0.marker_spin[spin])
            else:
                for g, qplist in zip(groups, qplist_group):
                    # Extract QP data.
                    yvals = [getattr(qp, what)[itemp] for qp in qplist]
                    label = "%s: %s" % (self._get_label(hue), g.hvalue)
                    ax.plot(g.xvalues, yvals, marker=nc0.marker_spin[spin],
                            label=label if i == 0 else None)

            ax.grid(True)
            ax.set_ylabel(what)
            if i == len(what_list) - 1:
                ax.set_xlabel("%s" % self._get_label(sortby))
                if sortby is None: rotate_ticklabels(ax, 15)
            if i == 0 and hue is not None:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        if "title" not in kwargs:
            title = "QP results spin: %s, k:%s, band: %s, T = %.1f K" % (
                    spin, repr(kpoint), band, nc0.tmesh[itemp])
            fig.suptitle(title, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_qpfield_vs_e0(self, field, itemp=0, sortby=None, hue=None, fontsize=8,
                           colormap="jet", e0="fermie", **kwargs):
        """
        For each file in the robot, plot one of the attributes of :class:`QpTempState`
        at temperature `itemp` as a function of the KS energy.

        Args:
            field (str): String defining the attribute to plot.
            itemp (int): Temperature index.

        .. note::

            For the meaning of the other arguments, see other robot methods.

        Returns: |matplotlib-Figure|
        """
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(colormap)

        if hue is None:
            ax_list = None
            lnp_list = self.sortby(sortby)
            for i, (label, ncfile, param) in enumerate(lnp_list):
                if sortby is not None:
                    label = "%s: %s" % (self._get_label(sortby), param)
                fig = ncfile.plot_qps_vs_e0(itemp_list=[itemp], with_fields=list_strings(field),
                    e0=e0, ax_list=ax_list, color=cmap(i/ len(lnp_list)), fontsize=fontsize,
                    label=label, show=False)
                ax_list = fig.axes
        else:
            # group_and_sortby and build (ngroups,) subplots
            groups = self.group_and_sortby(hue, sortby)
            nrows, ncols = 1, len(groups)
            ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                   sharex=True, sharey=True, squeeze=False)
            for ig, g in enumerate(groups):
                subtitle = "%s: %s" % (self._get_label(hue), g.hvalue)
                ax_mat[0, ig].set_title(subtitle, fontsize=fontsize)
                for i, (nclabel, ncfile, param) in enumerate(g):
                    fig = ncfile.plot_qps_vs_e0(itemp_list=[itemp], with_fields=list_strings(field),
                        e0=e0, ax_list=ax_mat[:, ig], color=cmap(i/ len(g)), fontsize=fontsize,
                        label="%s: %s" % (self._get_label(sortby), param), show=False)

                if ig != 0:
                    for ax in ax_mat[:, ig]:
                        set_visible(ax, False, "ylabel")

        return fig

    def write_notebook(self, nbpath=None, title=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=title)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.SigEPhRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("robot.get_params_dataframe()"),
            nbv.new_code_cell("# data = robot.get_dataframe()\ndata"),
            nbv.new_code_cell("robot.plot_qpgaps_convergence(itemp=0, sortby=None, hue=None);"),
            #nbv.new_code_cell("robot.plot_qpgaps_t(sortby=None);"),
            nbv.new_code_cell("""\
nc0 = robot.abifiles[0]
for spin in range(nc0.nsppol):
    for ik, sigma_kpoint in enumerate(nc0.sigma_kpoints):
        for band in range(nc0.bstart_sk[spin, ik], nc0.bstop_sk[spin, ik]):
            robot.plot_qpdata_conv_skb(spin, sigma_kpoint, band, itemp=0, sortby=None, hue=None);"""),

            nbv.new_code_cell("""\
#nc0 = robot.abifiles[0]
#for spin in range(nc0.nsppol):
#    for ik, sigma_kpoint in enumerate(nc0.sigma_kpoints):
#        for band in range(nc0.bstart_sk[spin, ik], nc0.bstop_sk[spin, ik]):
#           robot.plot_selfenergy_conv(spin, sigma_kpoint, band, itemp=0, sortby=None);"),"""),
        ])

        # Mixins.
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)


class SigmaPhReader(ElectronsReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: SigmaPhReader
    """
    def __init__(self, path):
        super(SigmaPhReader, self).__init__(path)

        self.nsppol = self.read_dimvalue("nsppol")

        # Get important dimensions.
        self.nkcalc = self.read_dimvalue("nkcalc")
        self.ntemp = self.read_dimvalue("ntemp")
        self.nqbz = self.read_dimvalue("nqbz")
        self.nqibz = self.read_dimvalue("nqibz")
        self.ngqpt = self.read_value("ngqpt")

        # T and frequency meshes.
        self.ktmesh = self.read_value("kTmesh")
        self.tmesh = self.ktmesh / abu.kb_HaK
        self.ktmesh *= units.Ha_to_eV

        # The K-points where QP corrections have been calculated.
        structure = self.read_structure()
        self.sigma_kpoints = KpointList(structure.reciprocal_lattice, self.read_value("kcalc"))
        for kpoint in self.sigma_kpoints:
            name = structure.findname_in_hsym_stars(kpoint)
            kpoint.set_name(name)

        # [nsppol, nkcalc] arrays with index of KS bands computed.
        # Note conversion between Fortran and python convention.
        self.bstart_sk = self.read_value("bstart_ks") - 1
        self.nbcalc_sk = self.read_value("nbcalc_ks")
        self.bstop_sk = self.bstart_sk + self.nbcalc_sk

        # Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
        # This dimension is optional, its presence signals that we have Sigma(w)
        self.nwr = self.read_dimvalue("nwr", default=0)

    def get_sigma_skb_kpoint(self, spin, kpoint, band):
        """
        Check k-point, band and spin index. Raise ValueError if invalid.
        """
        ik = self.sigkpt2index(kpoint)

        if not (len(self.sigma_kpoints) > ik >= 0):
            raise ValueError("Invalid k-point index %d. should be in [0, %d[" % (ik, len(self.sigma_kpoints)))
        if not (self.nsppol > spin >= 0):
            raise ValueError("Invalid spin index %d. should be in [0, %d[" % (ik, self.nsppol))
        if not (self.bstop_sk[spin, ik] > band >= self.bstart_sk[spin, ik]):
            raise ValueError("Invalid band index %d. should be in [%d, %d[" % (
                band, self.bstart_sk[spin, ik], self.bstop_sk[spin, ik]))

        return spin, ik, band - self.bstart_sk[spin, ik], self.sigma_kpoints[ik]

    def sigkpt2index(self, sigma_kpoint):
        """
        Returns the index of the self-energy k-point in sigma_kpoints
        Used to access data in the arrays that are dimensioned [0:nkcalc]
        sigma_kpoint can be either an integer or list with reduced coordinates.
        """
        if duck.is_intlike(sigma_kpoint):
            ik = int(sigma_kpoint)
            if self.nkcalc > ik >= 0: return ik
            raise ValueError("kpoint index should be in [0, %d] but received: %d" % (self.nkcalc, ik))
        else:
            return self.sigma_kpoints.index(sigma_kpoint)

    def read_qplist_sk(self, spin, kpoint, ignore_imag=False):
        """
        Read and return :class:`QpTempList` object for the given spin, kpoint.

        Args:
            spin: Spin index.
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            ignore_imag: Only real part is returned if ``ignore_imag``.
        """
        ik = self.sigkpt2index(kpoint)
        bstart, bstop = self.bstart_sk[spin, ik], self.bstop_sk[spin, ik]

        return QpTempList([self.read_qp(spin, ik, band, ignore_imag=ignore_imag)
                           for band in range(bstart, bstop)])

    def read_sigeph_skb(self, spin, kpoint, band):
        """
        Returns the e-ph self energy for the given (spin, k-point, band).

        Args:
            spin: Spin index
            kpoint: K-point in self-energy. Accepts |Kpoint|, vector or index.
            band: band index.

        Return: :class:`EphSelfEnergy` object.
        """
        if self.nwr == 0:
            raise ValueError("%s does not contain spectral function data." % self.path)

        spin, ik, ib, kpoint = self.get_sigma_skb_kpoint(spin, kpoint, band)

        # Abinit fortran (Ha units)
        # wrmesh_b(nwr, max_nbcalc, nkcalc, nsppol)
        # Frequency mesh along the real axis (Ha units) used for the different bands
        #print(spin, ik, ib, self.read_variable("wrmesh_b").shape)
        wmesh = self.read_variable("wrmesh_b")[spin, ik, ib, :] * units.Ha_to_eV

        # complex(dpc) :: vals_e0ks(ntemp, max_nbcalc, nkcalc, nsppol)
        # Sigma_eph(omega=eKS, kT, band)
        vals_e0ks = self.read_variable("vals_e0ks")[spin, ik, ib, :, :] * units.Ha_to_eV
        vals_e0ks = vals_e0ks[:, 0] + 1j * vals_e0ks[:, 1]

        # complex(dpc) :: dvals_de0ks(ntemp, max_nbcalc, nkcalc, nsppol)
        # d Sigma_eph(omega, kT, band, kcalc, spin) / d omega (omega=eKS)
        dvals_de0ks = self.read_variable("dvals_de0ks")[spin, ik, ib, :, :] * units.Ha_to_eV
        dvals_de0ks = dvals_de0ks[:, 0] + 1j * dvals_de0ks[:, 1]

        # real(dp) :: dw_vals(ntemp, max_nbcalc, nkcalc, nsppol)
        # Debye-Waller term (static).
        dw_vals = self.read_variable("dw_vals")[spin, ik, ib, :] * units.Ha_to_eV

        # complex(dpc) :: vals_wr(nwr, ntemp, max_nbcalc, nkcalc, nsppol)
        # Sigma_eph(omega, kT, band) for given (k, spin).
        # Note: enk_KS corresponds to nwr/2 + 1.
        vals_wr = self.read_variable("vals_wr")[spin, ik, ib, :, :, :] * units.Ha_to_eV
        vals_wr = vals_wr[:, :, 0] + 1j * vals_wr[:, :, 1]

        # Spectral function
        # nctkarr_t("spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol")
        spfunc_wr = self.read_variable("spfunc_wr")[spin, ik, ib, :, :] / units.Ha_to_eV

        # Read QP data. Note band instead of ib index.
        qp = self.read_qp(spin, ik, band)

        return EphSelfEnergy(wmesh, qp, vals_e0ks, dvals_de0ks, dw_vals, vals_wr, spfunc_wr)

    def read_qp(self, spin, kpoint, band, ignore_imag=False):
        """
        Return :class:`QpTempState` for the given (spin, kpoint, band).
        Only real part is returned if ``ignore_imag``.
        """
        spin, ik, ib, kpoint = self.get_sigma_skb_kpoint(spin, kpoint, band)

        def ri(a):
            return np.real(a) if ignore_imag else a

        # complex(dpc),allocatable :: qp_enes(:,:)
        # qp_enes(ntemp, max_nbcalc)
        # (Complex) QP energies computed with the non-adiabatic formalism.
        # nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol")
        var = self.read_variable("qp_enes")
        qpe = (var[spin, ik, ib, :, 0] + 1j * var[spin, ik, ib, :, 1]) * units.Ha_to_eV

        #var = self.read_variable("qpadb_enes")
        #qpe_adb = var[spin, ik, ib, :] * units.Ha_to_eV

        # nctkarr_t("dw_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol"),
        # Debye-Waller term (static).

        # complex(dpc),allocatable :: vals_e0ks(:,:)
        # vals_e0ks(ntemp, max_nbcalc)
        # Sigma_eph(omega=eKS, kT, band) for fixed (kcalc, spin).
        # TODO: Add Fan0 instead of computing Sigm - DW?

        e0 = self.read_variable("ks_enes")[spin, ik, ib] * units.Ha_to_eV
        ze0 = self.read_variable("ze0_vals")[spin, ik, ib]

        return QpTempState(
            spin=spin,
            kpoint=kpoint,
            band=band,
            tmesh=self.tmesh,
            e0=e0,
            qpe=ri(qpe),
            #sigxme=self._sigxme[spin, ik, ib_file],
            #sigcmee0=ri(self._sigcmee0[spin, ik, ib_file]),
            ze0=ze0,
        )

    def read_allqps(self, ignore_imag=False):
        """
        Return list with ``nsppol`` items. Each item is a :class:`QpTempList` with the QP results.

        Args:
            ignore_imag: Only real part is returned if ``ignore_imag``.
        """
        qps_spin = self.nsppol * [None]

        for spin in range(self.nsppol):
            qps = []
            for ik, kpoint in enumerate(self.sigma_kpoints):
                for band in range(self.bstart_sk[spin, ik], self.bstop_sk[spin, ik]):
                    qps.append(self.read_qp(spin, ik, band, ignore_imag=ignore_imag))
            qps_spin[spin] = QpTempList(qps)

        return tuple(qps_spin)
