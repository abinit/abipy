# coding: utf-8
"""
This module contains objects for the postprocessing of Sigma_eph calculations.

For a theoretical introduction see :cite:`Giustino2017`

Warning:

    Work in progress, DO NOT USE THIS CODE
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import tempfile
import pickle
import numpy as np
import pandas as pd
import pymatgen.core.units as units
import abipy.core.abinit_units as abu

from collections import OrderedDict, namedtuple, Iterable
from tabulate import tabulate
from monty.string import marquee, list_strings
from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core.mixins import AbinitNcFile, Has_Structure, Has_ElectronBands, NotebookWriter
from abipy.core.kpoints import Kpoint, KpointList, Kpath, IrredZone, has_timrev_from_kptopt
from abipy.tools.plotting import (add_fig_kwargs, get_ax_fig_plt, get_axarray_fig_plt, set_axlims, set_visible,
    rotate_ticklabels, ax_append_title)
from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands, RobotWithEbands, ElectronBandsPlotter, ElectronDosPlotter
#from abipy.dfpt.phonons import PhononBands, RobotWithPhbands, factor_ev2units, unit_tag, dos_label_from_units
from abipy.abio.robots import Robot
from abipy.eph.common import BaseEphReader


# TODO QPState and QPList from electrons.gw (Define base abstract class?).
# __eq__ based on skb?
# broadening, adiabatic ??

class QpTempState(namedtuple("QpTempState", "spin kpoint band tmesh e0 qpe ze0 fan0 dw qpe_adb")):
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
        fan0: Fan term (complex) at e_KS
        dw: Debye-Waller
        qpe_adb: Quasiparticle energy (real) in the adiabatic approximation.

    .. note::

        Energies are in eV.
    """

    @lazy_property
    def qpeme0(self):
        """E_QP[T] - E_0"""
        return self.qpe - self.e0

    @lazy_property
    def re_qpe(self):
        """Real part of the QP energy."""
        return self.qpe.real

    @lazy_property
    def imag_qpe(self):
        """Imaginay part of the QP energy."""
        return self.qpe.imag

    @property
    def re_fan0(self):
        """Real part of the Fan term at KS."""
        return self.fan0.real

    @property
    def imag_fan0(self):
        """Imaginary part of the Fan term at KS."""
        return self.fan0.imag

    @lazy_property
    def skb(self):
        """Tuple with (spin, kpoint, band)"""
        return self.spin, self.kpoint, self.band

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

    @classmethod
    def get_fields_for_plot(cls, vs, with_fields, exclude_fields):
        """
        Return list of QpTempState fields to plot from input arguments.

        Args:
            vs in ["temp", "e0"]
            with_fields:
            exclude_fields:
        """
        if vs == "temp":
            all_fields = list(cls.get_fields(exclude=["spin", "kpoint", "band", "e0", "tmesh"]))
        elif vs == "e0":
            all_fields = list(cls.get_fields(exclude=["spin", "kpoint", "band", "e0", "tmesh"]))
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

    @add_fig_kwargs
    def plot(self, with_fields="all", exclude_fields=None, ax_list=None, label=None, fontsize=12, **kwargs):
        """
        Plot the QP results as function of temperature.

        Args:
            with_fields: The names of the QpTempState attributes to plot as function of e0.
                Accepts: List of strings or string with tokens separated by blanks.
                See :class:`QpTempState` for the list of available fields.
            exclude_fields: Similar to `with_field` but excludes fields.
            ax_list: List of matplotlib axes for plot. If None, new figure is produced.
            label: Label for plot.
            fontsize: Fontsize for legend and title.

        Returns: |matplotlib-Figure|
        """
        fields = self.get_fields_for_plot("temp", with_fields, exclude_fields)
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
        for ix, (field, ax) in enumerate(zip(fields, ax_list)):
            irow, icol = divmod(ix, ncols)
            ax.grid(True)
            if irow == nrows - 1: ax.set_xlabel("Temperature [K]")
            ax.set_ylabel(field)
            yy = getattr(self, field)
            lbl = label if ix == 0 and label is not None else None

            # Handle complex arrays
            #if np.iscomplexobj(yy):
            #    ax.plot(self.tmesh, yy.real, linestyle, label=lbl, **kwargs)
            #    ax.plot(self.tmesh, yy.imag, linestyle, label=lbl, **kwargs)
            #else:
            ax.plot(self.tmesh, yy.real, linestyle, label=lbl, **kwargs)

        # Get around a bug in matplotlib
        if num_plots % ncols != 0: ax_list[-1].axis('off')

        if lbl is not None:
            ax_list[0].legend(loc="best", fontsize=fontsize, shadow=True)

        #fig.tight_layout()
        return fig

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

    # TODO: Linewidths
    @add_fig_kwargs
    def plot_vs_e0(self, itemp_list=None, with_fields="all", exclude_fields=None, fermie=None,
                   colormap="jet", ax_list=None, xlims=None, fontsize=12, **kwargs):
        """
        Plot the QP results as a function of the initial KS energy.

        Args:
            itemp_list: List of integers to select a particular temperature. None for all
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
        fields = QpTempState.get_fields_for_plot("e0", with_fields, exclude_fields)
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
        xlabel = r"$\epsilon_{KS}\;(eV)$"
        if fermie is not None:
            e0mesh -= fermie
            xlabel = r"$\epsilon_{KS}-\epsilon_F\;(eV)$"

        kw_linestyle = kwargs.pop("linestyle", "o")
        kw_color = kwargs.pop("color", None)
        kw_label = kwargs.pop("label", None)

        itemp_list = list(range(self.ntemp)) if itemp_list is None else duck.list_ints(itemp_list)
        for ix, (field, ax) in enumerate(zip(fields, ax_list)):
            irow, icol = divmod(ix, ncols)
            ax.grid(True)
            if irow == nrows - 1:
                ax.set_xlabel(xlabel)
            ax.set_ylabel(field, fontsize=fontsize)
            has_legend = False
            # Plot different temperatures.
            for itemp in itemp_list:
                yt = qps.get_field_itemp(field, itemp)
                # TODO real and imag?
                label = kw_label
                if kw_label is None:
                    label = "T = %.1f K" % self.tmesh[itemp] if ix == 0 else None
                has_legend = has_legend or bool(label)
                ax.plot(e0mesh, yt.real, kw_linestyle,
                        color=cmap(itemp / self.ntemp) if kw_color is None else kw_color,
                        label=label, **kwargs)

            set_axlims(ax, xlims, "x")
            #print("ix", ix, "has", has_legend, "label", label)
            if ix == 0 and has_legend:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        # Get around a bug in matplotlib
        if num_plots % ncols != 0: ax_list[-1].axis('off')

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
        app("Number of temperatures: %d, from %.1f to %.1f (K)" % (self.ntemp, self.tmesh[0], self.tmesh[-1]))
        app("Number of frequencies: %d, from %.1f to %.1f (eV)" % (self.nwr, self.wmesh[0], self.wmesh[-1]))
        app(self.qp.to_string(verbose=verbose, title="QP data"))

        return "\n".join(lines)

    def _get_wmesh_xlabel(self, zero_energy):
        """Return (wmesh, xlabel) from zero_energy input argument."""
        if zero_energy is None:
            xx = self.wmesh
            xlabel = r"$\omega\;(eV)$"
        elif zero_energy == "e0":
            xx = self.wmesh - self.qp.e0
            xlabel = r"$\omega - \epsilon_{nk}\;(eV)$"
        # TODO: chemical potential? but then I have mu(T) to handle in plots!
        #elif zero_energy == "fermie":
        #    xx = self.wmesh - self.fermie
        #    xlabel = r"$\omega\;(eV)$"
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
        for ix, (what, ax) in enumerate(zip(what_list, ax_list)):
            ax.grid(True)
            ax.set_ylabel(self.latex_symbol[what])
            if (ix == len(ax_list) - 1): ax.set_xlabel(xlabel)
            for itemp in itemps:
                ax.plot(xs, self._get_ys_itemp(what, itemp),
                        color=cmap(itemp / self.ntemp) if kw_color is None else kw_color,
                        label=tlabels[itemp] if (ix == 0 and kw_label is None) else kw_label,
                )
            if ix == 0: ax.legend(loc="best", shadow=True, fontsize=fontsize)
            set_axlims(ax, xlims, "x")

        if "title" not in kwargs:
            title = "K-point: %s, band: %d, spin: %d" % (repr(self.kpoint), self.band, self.spin)
            fig.suptitle(title, fontsize=fontsize)

        return fig

    plot = plot_tdep


class _MyQpkindsList(list):
    """Returned by find_qpkinds."""


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

    """
    def get_fundamental_gaps(self):

        ib_lumo = self.ebands.nelect // 2
        ib_homo = ib_lumo - 1

        # nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol")
        qpes = self.reader.read_value("qp_enes", cmode="c").real * units.Ha_to_eV
        for spin in range(self.nsppol):
            for ikc, kpoint in enumerate(self.sigma_kpoints):
                qpes[spin, ikc, :, :]

        def difference_matrix(a):
            x = np.reshape(a, (len(a), 1))
            return x - x.transpose()

        for spin, kset in enumerate(self.ebands.fundamental_gaps):
            ks_fgap = kset.energy
            # Find index in nkcalc
            ik_homo = self.sigkpt2index(kset.in_state.kpoint)
            ik_lumo = self.sigkpt2index(kset.out_state.kpoint)
            # nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol")
            qpes = self.reader.read_value("qp_enes", cmode="c") * units.Ha_to_eV
            qpes[spin, ik_homo, ib_homo].real
            qpes[spin, ik_lumo, ib_lumo].real
    """

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
        app("Max bstart: %d, min bstop: %d" % (self.reader.max_bstart, self.reader.min_bstop))
        app("Q-mesh: nqibz: %s, nqbz: %s, ngqpt: %s" % (self.nqibz, self.nqbz, str(self.ngqpt)))
        app("K-mesh for electrons:")
        app(self.ebands.kpoints.ksampling.to_string(verbose=verbose))
        app("Number of bands included in self-energy: %d" % (self.nbsum))
        app("zcut: %.3f [Ha], %.3f (eV)" % (self.zcut, self.zcut * units.Ha_to_eV))
        app("Number of temperatures: %d, from %.1f to %.1f (K)" % (self.ntemp, self.tmesh[0], self.tmesh[-1]))
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
            for ikc, kpoint in enumerate(self.sigma_kpoints):
                for spin in range(self.nsppol):
                    ks_gap = self.ks_dirgaps[spin, ikc]
                    qp_gap = self.qp_dirgaps_t[spin, ikc, it]
                    #qp_gap_adb = self.qp_dirgaps_adb_t[spin, ikc, it]
                    data.append([spin, repr(kpoint), ks_gap, qp_gap, qp_gap - ks_gap])
            app(str(tabulate(data, headers=["Spin", "k-point", "KS_gap", "QP_gap", "QP - KS"], floatfmt=".3f")))
            app("")

        if verbose > 1:
            app("K-points and bands included in self-energy corrections:")
            for spin in range(self.nsppol):
                for ikc, kpoint in enumerate(self.sigma_kpoints):
                    post = "ikc: %d" % (ikc if self.nsppol == 1 else "ikc: %d, spin: %d" % (ikc, spin))
                    app("\t%s: bstart: %d, bstop: %d, %s" % (
                        repr(kpoint), self.bstart_sk[spin, ikc], self.bstop_sk[spin, ikc], post))

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

    def sigkpt2index(self, kpoint):
        """
        Returns the index of the self-energy k-point in sigma_kpoints
        Used to access data in the arrays that are dimensioned with [0:nkcalc]
        """
        return self.reader.sigkpt2index(kpoint)

    def find_qpkinds(self, qp_kpoints):
        """
        Find kpoints for QP corrections from user input.
        Return list of (kpt, ikcalc) tuples where kpt is a |Kpoint| and
        ikcalc is the index in the nkcalc array..
        """
        if isinstance(qp_kpoints, _MyQpkindsList):
            return qp_kpoints

        if isinstance(qp_kpoints, Kpoint):
            qp_kpoints = [qp_kpoints]

        if qp_kpoints is None or (duck.is_string(qp_kpoints) and qp_kpoints == "all"):
            # qp_kpoints in (None, "all")
            items = self.sigma_kpoints, list(range(self.nkpt))

        elif duck.is_intlike(qp_kpoints):
            # qp_kpoints = 1
            ikc = int(qp_kpoints)
            items = [self.sigma_kpoints[ikc]], [ikc]

        elif isinstance(qp_kpoints, Iterable):
            # either [0, 1] or [[0, 0.5, 0]]
            # note possible ambiguity with [0, 0, 0] that is treated as list of integers.
            if duck.is_intlike(qp_kpoints[0]):
                ik_list = duck.list_ints(qp_kpoints)
                items = [self.sigma_kpoints[ikc] for ikc in ik_list], ik_list
            else:
                ik_list = [self.reader.sigkpt2index(kpt) for kpt in qp_kpoints]
                qp_kpoints = [self.sigma_kpoints[ikc] for ikc in ik_list]
                items = qp_kpoints, ik_list
        else:
            raise TypeError("Dont' know how to interpret `%s`" % (type(qp_kpoints)))

        # Check indices
        errors = []
        eapp = errors.append
        for ikc in items[1]:
            if ikc >= self.nkcalc:
                eapp("K-point index %d >= nkcalc %d, check input qp_kpoints" % (ikc, self.nkcalc))
        if errors:
            raise ValueError("\n".join(errors))

        return _MyQpkindsList(zip(items[0], items[1]))

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with the convergence parameters, e.g. ``nbsum``."""
        od = OrderedDict([
            ("nbsum", self.nbsum),
            ("zcut", self.zcut),
            ("symsigma", self.symsigma),
            # TODO: Remove?
            ("nqbz", self.reader.nqbz),
            ("nqibz", self.reader.nqibz),
        ])
        # Add EPH parameters.
        od.update(self.reader.common_eph_params)

        return od

    def get_sigeph_skb(self, spin, kpoint, band):
        return self.reader.read_sigeph_skb(spin, kpoint, band)

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
            for ikc, kpoint in enumerate(self.sigma_kpoints):
                app(self.get_dataframe_sk(spin, ikc, with_params=with_params, ignore_imag=ignore_imag))

        return pd.concat(df_list)

    def get_dataframe_sk(self, spin, kpoint, itemp=None, index=None, with_params=False, ignore_imag=False):
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
        ikc = self.sigkpt2index(kpoint)
        rows = []
        for band in range(self.bstart_sk[spin, ikc], self.bstop_sk[spin, ikc]):
            # Read QP data.
            qp = self.reader.read_qp(spin, ikc, band, ignore_imag=ignore_imag)
            # Convert to dataframe and add other entries useful when comparing different calculations.
            rows.append(qp.get_dataframe(params=self.params if with_params else None))

        df = pd.concat(rows)
        if itemp is not None: df = df[df["tmesh"] == self.tmesh[itemp]]
        return df

    #def get_dirgaps_dataframe(self):

    def interpolate(self, itemp_list=None, lpratio=5, ks_ebands_kpath=None, ks_ebands_kmesh=None, ks_degatol=1e-4,
                    vertices_names=None, line_density=20, filter_params=None, only_corrections=False, verbose=0): # pragma: no cover
        """
        Interpolated the self-energy corrections in k-space on a k-path and, optionally, on a k-mesh.

        Args:
            itemp_list: List of temperature indices to interpolate. None for all.
            lpratio: Ratio between the number of star functions and the number of ab-initio k-points.
                The default should be OK in many systems, larger values may be required for accurate derivatives.
            ks_ebands_kpath: KS |ElectronBands| on a k-path. If present,
                the routine interpolates the QP corrections and apply them on top of the KS band structure
                This is the recommended option because QP corrections are usually smoother than the
                QP energies and therefore easier to interpolate. If None, the QP energies are interpolated
                along the path defined by ``vertices_names`` and ``line_density``.
            ks_ebands_kmesh: KS |ElectronBands| on a homogeneous k-mesh. If present, the routine
                interpolates the corrections on the k-mesh (used to compute the QP DOS)
            ks_degatol: Energy tolerance in eV. Used when either ``ks_ebands_kpath`` or ``ks_ebands_kmesh`` are given.
                KS energies are assumed to be degenerate if they differ by less than this value.
                The interpolator may break band degeneracies (the error is usually smaller if QP corrections
                are interpolated instead of QP energies). This problem can be partly solved by averaging
                the interpolated values over the set of KS degenerate states.
                A negative value disables this ad-hoc symmetrization.
            vertices_names: Used to specify the k-path for the interpolated QP band structure
                when ``ks_ebands_kpath`` is None.
                It's a list of tuple, each tuple is of the form (kfrac_coords, kname) where
                kfrac_coords are the reduced coordinates of the k-point and kname is a string with the name of
                the k-point. Each point represents a vertex of the k-path. ``line_density`` defines
                the density of the sampling. If None, the k-path is automatically generated according
                to the point group of the system.
            line_density: Number of points in the smallest segment of the k-path. Used with ``vertices_names``.
            filter_params: TO BE DESCRIBED
            only_corrections: If True, the output contains the interpolated QP corrections instead of the QP energies.
                Available only if ks_ebands_kpath and/or ks_ebands_kmesh are used.
            verbose: Verbosity level

        Returns: class:`TdepElectronBands`.
        """
        # TODO: Consistency check.
        errlines = []
        eapp = errlines.append
        if len(self.sigma_kpoints) != len(self.ebands.kpoints):
            eapp("QP energies should be computed for all k-points in the IBZ but nkibz != nkptgw")
        if len(self.sigma_kpoints) == 1:
            eapp("QP Interpolation requires nkptgw > 1.")
        #if (np.any(self.bstop_sk[0, 0] != self.gwbstop_sk):
        #    cprint("Highest bdgw band is not constant over k-points. QP Bands will be interpolated up to...")
        #if (np.any(self.gwbstart_sk[0, 0] != self.gwbstart_sk):
        #if (np.any(self.gwbstart_sk[0, 0] != 0):
        if errlines:
            raise ValueError("\n".join(errlines))

        # Get symmetries from abinit spacegroup (read from file).
        abispg = self.structure.abi_spacegroup
        fm_symrel = [s for (s, afm) in zip(abispg.symrel, abispg.symafm) if afm == 1]

        if ks_ebands_kpath is None:
            # Generate k-points for interpolation. Will interpolate all bands available in the sigeph file.
            bstart, bstop = 0, -1
            if vertices_names is None:
                vertices_names = [(k.frac_coords, k.name) for k in self.structure.hsym_kpoints]
            kpath = Kpath.from_vertices_and_names(self.structure, vertices_names, line_density=line_density)
            kfrac_coords, knames = kpath.frac_coords, kpath.names

        else:
            # Use list of k-points from ks_ebands_kpath.
            ks_ebands_kpath = ElectronBands.as_ebands(ks_ebands_kpath)
            kfrac_coords = [k.frac_coords for k in ks_ebands_kpath.kpoints]
            knames = [k.name for k in ks_ebands_kpath.kpoints]

            # Find the band range for the interpolation.
            bstart, bstop = 0, ks_ebands_kpath.nband
            # FIXME what about bstart?
            bstop = min(bstop, self.reader.min_bstop)
            if ks_ebands_kpath.nband < self.reader.min_bstop:
                cprint("Number of bands in KS band structure smaller than the number of bands in GW corrections", "red")
                cprint("Highest GW bands will be ignored", "red")

        # Interpolate QP energies if ks_ebands_kpath is None else interpolate QP corrections
        # and re-apply them on top of the KS band structure.
        gw_kcoords = [k.frac_coords for k in self.sigma_kpoints]

        # Read QP energies from file (real + imag part) and compute corrections if ks_ebands_kpath.
        # nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol")
        qpes = self.reader.read_value("qp_enes", cmode="c") * units.Ha_to_eV

        if ks_ebands_kpath is not None:
            if ks_ebands_kpath.structure != self.structure:
                cprint("sigres.structure and ks_ebands_kpath.structures differ. Check your files!", "red")
            # nctkarr_t("ks_enes", "dp", "max_nbcalc, nkcalc, nsppol")
            ks_enes = self.reader.read_value("ks_enes") * units.Ha_to_eV
            for itemp in range(self.ntemp):
                qpes[:, :, :, itemp] -= ks_enes

        # Note there's no guarantee that the sigma_kpoints and the corrections have the same k-point index.
        # Be careful because the order of the k-points and the band range stored in the SIGRES file may differ ...
        # FIXME: This part is not clear to me now !
        #qpdata = np.empty(qpes.shape)
        #for gwk in self.sigma_kpoints:
        #    ik_ibz = self.reader.kpt2fileindex(gwk)
        #    for spin in range(self.nsppol):
        #        qpdata[spin, ik_ibz, :] = egw_rarr[spin, ik_ibz, :]

        # Build interpolator for QP corrections.
        from abipy.core.skw import SkwInterpolator
        cell = (self.structure.lattice.matrix, self.structure.frac_coords, self.structure.atomic_numbers)
        has_timrev = has_timrev_from_kptopt(self.reader.read_value("kptopt"))

        qp_ebands_kpath_t, qp_ebands_kmesh_t, interpolators_t = [], [], []
        itemp_list = list(range(self.ntemp)) if itemp_list is None else duck.list_ints(itemp_list)
        for itemp in itemp_list:
            skw_reim = []
            for reim in ("real", "imag"):
                qpdata = qpes[:, :, bstart:bstop, itemp]
                qpdata = getattr(qpdata, reim).copy()

                skw = SkwInterpolator(lpratio, gw_kcoords, qpdata, self.ebands.fermie, self.ebands.nelect,
                                      cell, fm_symrel, has_timrev,
                                      filter_params=filter_params, verbose=verbose)
                skw_reim.append(skw)

                if ks_ebands_kpath is None:
                    # Interpolate QP energies.
                    if reim == "real":
                        eigens_kpath = skw.interp_kpts(kfrac_coords).eigens
                    else:
                        lw_kpath = skw.interp_kpts(kfrac_coords).eigens
                else:
                    # Interpolate QP energies corrections and add them to KS.
                    ref_eigens = ks_ebands_kpath.eigens[:, :, bstart:bstop]
                    qp_corrs = skw.interp_kpts_and_enforce_degs(kfrac_coords, ref_eigens, atol=ks_degatol).eigens
                    if reim == "real":
                        eigens_kpath = qp_corrs if only_corrections else ref_eigens + qp_corrs
                    else:
                        lw_kpath = qp_corrs

            interpolators_t.append(skw_reim)

            # Build new ebands object with k-path.
            kpts_kpath = Kpath(self.structure.reciprocal_lattice, kfrac_coords, weights=None, names=knames)
            occfacts_kpath = np.zeros(eigens_kpath.shape)

            # Finding the new Fermi level of the interpolated bands is not trivial, in particular if metallic.
            # because one should first interpolate the QP bands on a mesh. Here I align the QP bands
            # at the HOMO of the KS bands.
            homos = ks_ebands_kpath.homos if ks_ebands_kpath is not None else self.ebands.homos
            qp_fermie = max([eigens_kpath[e.spin, e.kidx, e.band] for e in homos])
            #qp_fermie = self.ebands.fermie
            #qp_fermie = 0.0

            newt = ElectronBands(self.structure, kpts_kpath, eigens_kpath, qp_fermie, occfacts_kpath,
                                 self.ebands.nelect, self.ebands.nspinor, self.ebands.nspden,
                                 linewidths=lw_kpath)
            qp_ebands_kpath_t.append(newt)

        # TODO
        qp_ebands_kmesh = None
        if ks_ebands_kmesh is not None:
            raise NotImplementedError()
            # Interpolate QP corrections on the same k-mesh as the one used in the KS run.
            ks_ebands_kmesh = ElectronBands.as_ebands(ks_ebands_kmesh)
            if bstop > ks_ebands_kmesh.nband:
                raise ValueError("Not enough bands in ks_ebands_kmesh, found %s, minimum expected %d\n" % (
                    ks_ebands_kmesh%nband, bstop))
            if ks_ebands_kpath.structure != self.structure:
                raise ValueError("sigres.structure and ks_ebands_kmesh.structures differ. Check your files!")

            # K-points and weight for DOS are taken from ks_ebands_kmesh
            dos_kcoords = [k.frac_coords for k in ks_ebands_kmesh.kpoints]
            dos_weights = [k.weight for k in ks_ebands_kmesh.kpoints]

            # Interpolate QP corrections from bstart to bstop
            ref_eigens = ks_ebands_kmesh.eigens[:, :, bstart:bstop]
            qp_corrs = skw.interp_kpts_and_enforce_degs(dos_kcoords, ref_eigens, atol=ks_degatol).eigens
            eigens_kmesh = qp_corrs if only_corrections else ref_eigens + qp_corrs

            # Build new ebands object with k-mesh
            kpts_kmesh = IrredZone(self.structure.reciprocal_lattice, dos_kcoords, weights=dos_weights,
                                   names=None, ksampling=ks_ebands_kmesh.kpoints.ksampling)
            occfacts_kmesh = np.zeros(eigens_kmesh.shape)
            linewidths_kmesh = np.ones(eigens_kmesh.shape) * 0.2
            newt = ElectronBands(self.structure, kpts_kmesh, eigens_kmesh, qp_fermie, occfacts_kmesh,
                                 self.ebands.nelect, self.ebands.nspinor, self.ebands.nspden,
                                 linewidths=linewidths_kmesh)
            qp_ebands_kmesh_t.append(newt)

        return TdepElectronBands(self.tmesh[itemp_list], ks_ebands_kpath, qp_ebands_kpath_t,
                                 ks_ebands_kmesh, qp_ebands_kmesh_t, interpolators_t)

    @add_fig_kwargs
    def plot_qpgaps_t(self, qp_kpoints=0, ax_list=None, plot_qpmks=True, fontsize=8, **kwargs):
        """
        Plot the KS and the QP(T) direct gaps for all the k-points available in the SIGEPH file.

        Args:
            qp_kpoints: List of k-points in self-energy. Accept integers (list or scalars), list of vectors,
                or None to plot all k-points.
            ax_list: List of |matplotlib-Axes| or None if a new figure should be created.
            plot_qpmks: If False, plot QP_gap, KS_gap else (QP_gap - KS_gap)
            fontsize: legend and title fontsize.
            kwargs: Passed to ax.plot method except for marker.

        Returns: |matplotlib-Figure|
        """
        qpkinds = self.find_qpkinds(qp_kpoints)
        # Build grid plot.
        nrows, ncols = len(qpkinds), 1
        ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=True)
        ax_list = np.array(ax_list).ravel()
        label = kwargs.pop("label", None)

        for ix, ((kpt, ikc), ax) in enumerate(zip(qpkinds, ax_list)):
            for spin in range(self.nsppol):
                if not plot_qpmks:
                    # Plot QP_{spin,kpt}(T)
                    ax.plot(self.tmesh, self.qp_dirgaps_t[spin, ikc], marker=self.marker_spin[spin],
                            label=label, **kwargs) #label="QP gap")
                    # Add KS gap (assumed at T=0).
                    ax.scatter(0, self.ks_dirgaps[spin, ikc]) #, label="KS gap %s" % label)
                else:
                    # Plot QP_{spin,kpt}(T) - KS_gap
                    ax.plot(self.tmesh, self.qp_dirgaps_t[spin, ikc] - self.ks_dirgaps[spin, ikc],
                            marker=self.marker_spin[spin], label=label)

            ax.grid(True)
            if ix == len(qpkinds) - 1:
                ax.set_xlabel("Temperature (K)")
                if plot_qpmks:
                    ax.set_ylabel("QP-KS gap (eV)")
                else:
                    ax.set_ylabel("QP direct gap (eV)")
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
        what_list = ["re_qpe", "imag_qpe", "ze0"]   # "re_fan0", "imag_fan0", "dw"

        # Build grid plot.
        nrows, ncols = len(what_list), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        # Read all QPs for this (spin, kpoint) and all bands.
        qp_list = self.reader.read_qplist_sk(spin, kpoint)

        for ix, (ax, what) in enumerate(zip(ax_list, what_list)):
            # Plot QP(T)
            for qp in qp_list:
                if band_list is not None and qp.band not in band_list: continue
                yvals = getattr(qp, what)
                ax.plot(qp.tmesh, yvals, marker=self.marker_spin[spin],
                        label="band: %s" % qp.band)

            ax.grid(True)
            ax.set_ylabel(what)
            if ix == len(what_list) - 1: ax.set_xlabel("Temperature [K]")
            if ix == 0: ax.legend(loc="best", fontsize=fontsize, shadow=True)

        if "title" not in kwargs:
            title = "QP results spin:%s, k:%s" % (spin, repr(qp_list[0].kpoint))
            fig.suptitle(title, fontsize=fontsize)

        return fig

    @lazy_property
    def qplist_spin(self):
        """Tuple of :class:`QpTempList` objects indexed by spin."""
        return self.reader.read_allqps()

    @add_fig_kwargs
    def plot_qps_vs_e0(self, itemp_list=None, with_fields="all", exclude_fields=None, e0="fermie",
                       colormap="jet", xlims=None, ax_list=None, fontsize=8, **kwargs):
        """
        Plot the QP results in the SIGEPH file as function of the initial KS energy.

        Args:
            itemp_list: List of integers to select a particular temperature. None means all
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

        #for ix, ax in enumerate(ax_list):
        #    if ix != 0:
        #        set_visible(ax, False, "legend")
        #    else:
        #        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot_qpbands_ibzt(self, itemp_list=None, e0="fermie", colormap="jet", ylims=None, fontsize=8, **kwargs):
        r"""
        Plot the KS band structure in the IBZ with the QP(T) energies.

        Args:
            itemp_list: List of integers to select a particular temperature. None for all
            e0: Option used to define the zero of energy in the band structure plot.
            colormap: matplotlib color map.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used.
            fontsize: Legend and title fontsize.

        Returns: |matplotlib-Figure|
        """
        itemp_list = list(range(self.ntemp)) if itemp_list is None else duck.list_ints(itemp_list)

        # Map sigma_kpoints to ebands.kpoints
        kcalc2ibz = np.empty(self.nkcalc, dtype=np.int)
        for ikc, sigkpt in enumerate(self.sigma_kpoints):
            kcalc2ibz[ikc] = self.ebands.kpoints.index(sigkpt)

        # TODO: It seems there's a minor issue with fermie if SCF band structure.
        e0 = self.ebands.get_e0(e0)
        #print("e0",e0, self.ebands.fermie)

        # Build grid with (1, nsppol) plots.
        nrows, ncols = 1, self.nsppol
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = np.array(ax_list).ravel()
        cmap = plt.get_cmap(colormap)

        # Read QP energies: nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol")
        qpes = self.reader.read_value("qp_enes", cmode="c") * units.Ha_to_eV
        band_range = (self.reader.max_bstart, self.reader.min_bstop)

        for spin, ax in zip(range(self.nsppol), ax_list):
            # Plot KS bands in the band range included in self-energy calculation.
            self.ebands.plot(ax=ax, e0=e0, spin=spin, band_range=band_range, show=False)
            # Add (scattered) QP(T) energies for the calculated k-points.
            for itemp in itemp_list:
                yvals = qpes[spin, :, :, itemp].real - e0
                for band in range(self.reader.min_bstop):
                    ax.scatter(kcalc2ibz, yvals[:, band],
                        label="T = %.1f K" % self.tmesh[itemp] if band == 0 else None,
                        color=cmap(itemp / self.ntemp), alpha=0.6, marker="o", s=20,
                    )

            set_axlims(ax, ylims, "y")
            if spin == 0:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        verbose = kwargs.pop("verbose", 0)
        yield self.plot_qpbands_ibzt(show=False)
        #yield self.plot_qpgaps_t(qp_kpoints=0, show=False)
        for i, qp_kpt in enumerate(self.sigma_kpoints):
            if i > 5 and not verbose:
                print("File contains more than 5 k-points. Only the first five k-points are displayed.")
                break
            yield self.plot_qpgaps_t(qp_kpoints=qp_kpt, show=False)
        yield self.plot_qps_vs_e0(show=False)

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
                for ikc, kpoint in enumerate(ncfile.sigma_kpoints):
                    app(ncfile.get_dataframe_sk(spin, ikc, with_params=with_params, ignore_imag=ignore_imag))
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
            for ix, (label, ncfile, param) in enumerate(lnp_list):
                sigma = ncfile.reader.read_sigeph_skb(spin, kpoint, band)
                fig = sigma.plot_tdep(itemps=itemp, ax_list=ax_list,
                    label=label, color=cmap(ix / len(lnp_list)), show=False)
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
                for ix, (nclabel, ncfile, param) in enumerate(g):
                    sigma = ncfile.reader.read_sigeph_skb(spin, kpoint, band)
                    fig = sigma.plot_tdep(itemps=itemp, ax_list=ax_mat[:, ig],
                        label="%s: %s" % (self._get_label(sortby), param),
                        color=cmap(ix / len(g)), show=False)

            if ig != 0:
                for ax in ax_mat[:, ig]:
                    set_visible(ax, False, "ylabel")

            for ax in ax_mat.ravel():
                set_axlims(ax, xlims, "x")

        return fig

    @add_fig_kwargs
    def plot_qpgaps_t(self, qp_kpoints=0, plot_qpmks=True, sortby=None, hue=None, fontsize=8, **kwargs):
        """
        Compare the QP(T) direct gaps for all the k-points available in the robot.

        Args:
            qp_kpoints: List of k-points in self-energy. Accept integers (list or scalars), list of vectors,
                or None to plot all k-points.
            plot_qpmks: If False, plot QP_gap, KS_gap else (QP_gap - KS_gap)
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
            for ix, ((label, ncfile, param), lineopt) in enumerate(zip(self.sortby(sortby), self.iter_lineopt())):
                fig = ncfile.plot_qpgaps_t(qp_kpoints=qp_kpoints, ax_list=ax_list,
                    plot_qpmks=plot_qpmks, label=label, show=False, fontsize=fontsize, **lineopt)
                    #label=label if ix == 0 else None, show=False, fontsize=fontsize, **lineopt)
                ax_list = fig.axes
        else:
            # Need to know number of k-points here to build grid
            qpkinds = self.abifiles[0].find_qpkinds(qp_kpoints)
            # (nkpt, ngroups) subplots
            groups = self.group_and_sortby(hue, sortby)
            nrows, ncols = len(qpkinds), len(groups)
            ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                   sharex=True, sharey=False, squeeze=False)
            for ig, g in enumerate(groups):
                for ix, (nclabel, ncfile, param) in enumerate(g):
                    fig = ncfile.plot_qpgaps_t(qp_kpoints=qpkinds, ax_list=ax_mat[:, ig],
                            plot_qpmks=plot_qpmks, label="%s: %s" % (self._get_label(sortby), param),
                            fontsize=fontsize, show=False) #, **lineopt)

                # Add label with hue to previous title with k-point info.
                ax_append_title(ax_mat[0, ig], "\n%s: %s" % (self._get_label(hue), g.hvalue), fontsize=fontsize)
                if ig != 0:
                    for ax in ax_mat[:, ig]:
                        set_visible(ax, False, "ylabel")

            # Hide legends except the first one
            if nrows > 1:
                for ax in ax_mat[1:, :].ravel():
                    set_visible(ax, False, "legend")

        return fig

    @add_fig_kwargs
    def plot_qpgaps_convergence(self, qp_kpoints=0, itemp=0, sortby=None, hue=None,
                                plot_qpmks=True, fontsize=8, **kwargs):
        """
        Plot the convergence of the direct QP gaps at given temperature
        for all the k-points and spins treated by the robot.

        Args:
            qp_kpoints: List of k-points in self-energy. Accept integers (list or scalars), list of vectors,
                or None to plot all k-points.
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
            plot_qpmks: If False, plot QP_gap, KS_gap else (QP_gap - KS_gap)
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
        """
        # Make sure that nsppol, sigma_kpoints, and tmesh are consistent.
        self._check_dims_and_params()

        nc0 = self.abifiles[0]
        nsppol = nc0.nsppol
        qpkinds = nc0.find_qpkinds(qp_kpoints)

        # Build grid with (nkpt, 1) plots.
        nrows, ncols = len(qpkinds), 1
        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        if hue is None:
            labels, ncfiles, params = self.sortby(sortby, unpack=True)
        else:
            groups = self.group_and_sortby(hue, sortby)

        name = "QP dirgap" if not plot_qpmks else "QP-KS dirgap"
        for ix, ((kpt, ikc), ax) in enumerate(zip(qpkinds, ax_list)):
            for spin in range(nsppol):
                ax.set_title("%s k:%s, T = %.1f K" % (
                    name, repr(kpt), nc0.tmesh[itemp]), fontsize=fontsize)

                # Extract QP dirgap for [spin, kpt, itemp]
                if hue is None:
                    yvals = [ncfile.qp_dirgaps_t[spin, ikc, itemp] for ncfile in ncfiles]
                    if plot_qpmks:
                        yvals = np.array(yvals) - np.array([ncfile.ks_dirgaps[spin, ikc] for ncfile in ncfiles])

                    if not duck.is_string(params[0]):
                        ax.plot(params, yvals, marker=nc0.marker_spin[spin])
                    else:
                        # Must handle list of strings in a different way.
                        xn = range(len(params))
                        ax.plot(xn, yvals, marker=nc0.marker_spin[spin])
                        ax.set_xticks(xn)
                        ax.set_xticklabels(params, fontsize=fontsize)

                else:
                    for g in groups:
                        yvals = [ncfile.qp_dirgaps_t[spin, ikc, itemp] for ncfile in g.abifiles]
                        if plot_qpmks:
                            yvals = np.array(yvals) - np.array([ncfile.ks_dirgaps[spin, ikc] for ncfile in g.abifiles])
                        label = "%s: %s" % (self._get_label(hue), g.hvalue)
                        ax.plot(g.xvalues, yvals, marker=nc0.marker_spin[spin], label=label)

            ax.grid(True)
            if ix == len(qpkinds) - 1:
                ax.set_ylabel("%s (eV)" % name)
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
        ax_list = np.array(ax_list).ravel()

        nc0 = self.abifiles[0]
        ikc = nc0.sigkpt2index(kpoint)
        kpoint = nc0.sigma_kpoints[ikc]

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

        for ix, (ax, what) in enumerate(zip(ax_list, what_list)):
            if hue is None:
                # Extract QP data.
                yvals = [getattr(qp, what)[itemp] for qp in qplist]
                if not duck.is_string(params[0]):
                    ax.plot(params, yvals, marker=nc0.marker_spin[spin])
                else:
                    # Must handle list of strings in a different way.
                    xn = range(len(params))
                    ax.plot(xn, yvals, marker=nc0.marker_spin[spin])
                    ax.set_xticks(xn)
                    ax.set_xticklabels(params, fontsize=fontsize)
            else:
                for g, qplist in zip(groups, qplist_group):
                    # Extract QP data.
                    yvals = [getattr(qp, what)[itemp] for qp in qplist]
                    label = "%s: %s" % (self._get_label(hue), g.hvalue)
                    ax.plot(g.xvalues, yvals, marker=nc0.marker_spin[spin],
                            label=label if ix == 0 else None)

            ax.grid(True)
            ax.set_ylabel(what)
            if ix == len(what_list) - 1:
                ax.set_xlabel("%s" % self._get_label(sortby))
                if sortby is None: rotate_ticklabels(ax, 15)
            if ix == 0 and hue is not None:
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
            fontsize: legend and label fontsize.
            e0: Option used to define the zero of energy in the band structure plot.

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
                    e0=e0, ax_list=ax_list, color=cmap(i / len(lnp_list)), fontsize=fontsize,
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
                        e0=e0, ax_list=ax_mat[:, ig], color=cmap(i / len(g)), fontsize=fontsize,
                        label="%s: %s" % (self._get_label(sortby), param), show=False)

                if ig != 0:
                    for ax in ax_mat[:, ig]:
                        set_visible(ax, False, "ylabel")

        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        verbose = kwargs.pop("verbose", 0)
        itemp = 0
        for i, qp_kpt in enumerate(self.abifiles[0].sigma_kpoints):
            if i > 5 and not verbose:
                print("File contains more than 5 k-points. Only the first five k-points are displayed.")
                break
            yield self.plot_qpgaps_convergence(qp_kpoints=qp_kpt, itemp=itemp, show=False)
            yield self.plot_qpgaps_t(qp_kpoints=qp_kpt, show=False)

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
    for ikc, sigma_kpoint in enumerate(nc0.sigma_kpoints):
        for band in range(nc0.bstart_sk[spin, ikc], nc0.bstop_sk[spin, ikc]):
            robot.plot_qpdata_conv_skb(spin, sigma_kpoint, band, itemp=0, sortby=None, hue=None);"""),

            nbv.new_code_cell("""\
#nc0 = robot.abifiles[0]
#for spin in range(nc0.nsppol):
#    for ikc, sigma_kpoint in enumerate(nc0.sigma_kpoints):
#        for band in range(nc0.bstart_sk[spin, ikc], nc0.bstop_sk[spin, ikc]):
#           robot.plot_selfenergy_conv(spin, sigma_kpoint, band, itemp=0, sortby=None);"),"""),
        ])

        # Mixins.
        nb.cells.extend(self.get_baserobot_code_cells())
        nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)


class TdepElectronBands(object): # pragma: no cover
    """
    A list of |ElectronBands| (T) with a real part renormalized by
    the E-PH sel-energy. Imaginary part is stored in a specialized array.
    This object is not supposed to be instantiated directly.
    It is usually constructed in SigEPhFile by interpolating the ab-initio results

    Provides methods to plot renormalized band structures with linewidths.
    """
    def __init__(self, tmesh, ks_ebands_kpath, qp_ebands_kpath_t,
                 ks_ebands_kmesh, qp_ebands_kmesh_t, interpolators_t):
        """
        Args:
            tmesh: Array of temperatures in K.
            ks_ebands_kpath: KS bands on k-path (None if not available).
            qp_ebands_kpath_t: List of QP(T) bands on k-path (empty if not available).
            ks_ebands_kmesh: KS bands on k-mesh (None if not available).
            qp_ebands_kmesh_t: List of QP(T) bands on k-mesh (empty if not available).
            interpolators_t: [ntemp, 2] interpolators for real/imaginary part.
        """
        self.tmesh = np.array(tmesh)
        self.ntemp = self.tmesh.size
        self.interpolators_t = np.reshape(interpolators_t, (self.ntemp, 2))
        assert len(self.interpolators_t) == self.ntemp

        self.ks_ebands_kpath = ks_ebands_kpath
        self.qp_ebands_kpath_t = qp_ebands_kpath_t
        if qp_ebands_kpath_t:
            assert len(self.qp_ebands_kpath_t) == self.ntemp

        self.ks_ebands_kmesh = ks_ebands_kmesh
        self.qp_ebands_kmesh_t = qp_ebands_kmesh_t
        if qp_ebands_kmesh_t:
            assert len(self.qp_ebands_kmesh_t) == self.ntemp

    @lazy_property
    def has_kpath(self):
        """True if interpolated bands on the k-path are available."""
        return bool(self.qp_ebands_kpath_t)

    @lazy_property
    def has_kmesh(self):
        """True if interpolated bands on the k-mesh are available."""
        return bool(self.qp_ebands_kmesh_t)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation with verbosiy level ``verbose``."""
        lines = []
        app = lines.append
        app("Number of temperatures: %d, tmesh[0]: %.1f, tmesh[-1] %.1f [K]" % (
            self.ntemp, self.tmesh[0], self.tmesh[-1]))
        app("Has qp_bands on k-path: %s" % self.has_kpath)
        app("Has ks_bands on k-path: %s" % (self.ks_ebands_kpath is not None))
        app("Has qp_bands on k-mesh: %s" % self.has_kmesh)
        app("Has ks_bands on k-mesh: %s" % (self.ks_ebands_kmesh is not None))

        return "\n".join(lines)

    @classmethod
    def pickle_load(cls, filepath):
        """Loads the object from a pickle file."""
        with open(filepath, "rb") as fh:
            new = pickle.load(fh)
            #assert cls is new.__class__
            return new

    def pickle_dump(self, filepath=None):
        """
        Save the status of the object in pickle format.
        If filepath is None, a temporary file is created.

        Return: name of the pickle file.
        """
        if filepath is None:
            _, filepath = tempfile.mkstemp(suffix='.pickle')

        with open(filepath, "wb") as fh:
            pickle.dump(self, fh)
            return filepath

    @add_fig_kwargs
    def plot_itemp(self, itemp, ax=None, e0="fermie", ylims=None, fontsize=12, **kwargs):
        """
        Plot band structures with linewidth at temperature ``itemp``.

        Args:
            itemp: Temperature index.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: Fontsize for title.

        Return: |matplotlib-Figure|
        """
        #if not self.has_kpath: return None
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        # Plot KS bands
        if self.ks_ebands_kpath is not None:
            ks_opts = dict(color="black", lw=2, label="KS")
            self.ks_ebands_kpath.plot_ax(ax, e0, **ks_opts)

        # Plot QP(T) bands with linewidths
        title = "T = %.1f K" % self.tmesh[itemp]
        qp_opts = dict(color="red", lw=2, label="QP")
        self.qp_ebands_kpath_t[itemp].plot_ax(ax, e0, **qp_opts)
        self.qp_ebands_kpath_t[itemp].decorate_ax(ax, fontsize=fontsize, title=title)

        set_axlims(ax, ylims, "y")
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    @add_fig_kwargs
    def plot(self, e0="fermie", ylims=None, fontsize=8, **kwargs):
        """
        Plot grid of band structures with linewidth (one plot for each temperature).

        Args:

            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            fontsize: Fontsize for title.

        Return: |matplotlib-Figure|
        """
        num_plots, ncols, nrows = self.ntemp, 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=True, sharey=True, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        #e0 = 0
        for itemp, ax in enumerate(ax_list):
            fig = self.plot_itemp(itemp, ax=ax, e0=e0, ylims=ylims, fontsize=fontsize, show=False)
            if itemp != 0:
                set_visible(ax, False, "ylabel", "legend")

        return fig

    @add_fig_kwargs
    def plot_lws_vs_e0(self, itemp_list=None, ax=None, e0="fermie", exchange_xy=False,
                       colormap="jet", xlims=None, ylims=None, fontsize=8, **kwargs):
        r"""
        Plot the electronic linewidths vs KS energy at temperature ``itemp``

        Args:
            itemp_list: List of temperature indices to interpolate. None for all.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            e0: Option used to define the zero of energy in the band structure plot. Possible values:
                - ``fermie``: shift all eigenvalues to have zero energy at the Fermi energy (``self.fermie``).
                -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
                -  None: Don't shift energies, equivalent to e0=0
            exchange_xy: True to exchange x-y axis.
            colormap: matplotlib color map.
            xlims, ylims: Set the data limits for the x-axis or the y-axis. Accept tuple e.g. ``(left, right)``
                   or scalar e.g. ``left``. If left (right) is None, default values are used
            fontsize: fontsize for titles and legend.

        Returns: |matplotlib-Figure|
        """
        itemp_list = list(range(self.ntemp)) if itemp_list is None else duck.list_ints(itemp_list)

        ax, fig, plt = get_ax_fig_plt(ax=ax)
        cmap = plt.get_cmap(colormap)

        #kw_linestyle = kwargs.pop("linestyle", "o")
        #kw_color = kwargs.pop("color", None)
        #kw_label = kwargs.pop("label", None)

        for itemp in itemp_list:
            # Select kmesh or kpath
            if self.has_kmesh:
                qp_ebands = self.qp_ebands_mesh_t[itemp]
            else:
                qp_ebands = self.qp_ebands_kpath_t[itemp]

            fig = qp_ebands.plot_lws_vs_e0(ax=ax, e0=e0, exchange_xy=False,
                xlims=xlims, ylims=ylims, fontsize=fontsize,
                label="T = %.1f K" % self.tmesh[itemp],
                color=cmap(itemp / len(itemp_list)), # if kw_color is None else kw_color,
                show=False)

        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    def get_ebands_plotter(self, edos_kwargs=None, with_edos=True):
        """
        Build and return |ElectronBandsPlotter| with KS and QP(T) results

        Args:
            edos_kwargs: optional dictionary with the options passed to ``get_edos`` to compute the electron DOS.
            with_edos: False if DOSes are not wanted
        """
        if not self.has_kpath:
            raise ValueError("QP bands on k-path are required.")

        ebands_plotter = ElectronBandsPlotter()
        if self.ks_ebands_kpath is not None:
            ebands_plotter.add_ebands("KS", self.ks_ebands_kpath,
                edos=self.ks_ebands_kmesh if with_edos else None,
                edos_kwargs=edos_kwargs)

        for itemp in range(self.ntemp):
            label="T = %.1f K" % self.tmesh[itemp]
            ebands_plotter.add_ebands(label, self.qp_ebands_kpath_t[itemp],
                edos=self.qp_ebands_kmesh_t[itemp] if (self.has_kmesh and with_edos) else None,
                edos_kwargs=edos_kwargs)

        return ebands_plotter

    def get_edos_plotter(self, edos_kwargs=None):
        """
        Build and return |ElectronDosPlotter| with KS and QP(T) results.

        Args:
            edos_kwargs: optional dictionary with the options passed to ``get_edos`` to compute the electron DOS.
        """
        if not self.has_kmesh:
            raise ValueError("QP bands on k-mesh are required.")

        edos_plotter = ElectronDosPlotter()
        if self.ks_ebands_kmesh is not None:
            edos_plotter.add_edos("KS", edos=self.ks_ebands_kmesh, edos_kwargs=edos_kwargs)

        for itemp in range(self.ntemp):
            label="T = %.1f K" % self.tmesh[itemp]
            edos_plotter.add_edos(label, edos=self.qp_ebands_kmesh_t[itemp], edos_kwargs=edos_kwargs)

        return edos_plotter


class SigmaPhReader(BaseEphReader):
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
        self.max_bstart = self.bstart_sk.max()
        self.min_bstop = self.bstop_sk.min()

        # Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
        # This dimension is optional, its presence signals that we have Sigma(w)
        self.nwr = self.read_dimvalue("nwr", default=0)

    def get_sigma_skb_kpoint(self, spin, kpoint, band):
        """
        Check k-point, band and spin index. Raise ValueError if invalid.
        """
        ikc = self.sigkpt2index(kpoint)

        if not (len(self.sigma_kpoints) > ikc >= 0):
            raise ValueError("Invalid k-point index %d. should be in [0, %d[" % (ikc, len(self.sigma_kpoints)))
        if not (self.nsppol > spin >= 0):
            raise ValueError("Invalid spin index %d. should be in [0, %d[" % (ikc, self.nsppol))
        if not (self.bstop_sk[spin, ikc] > band >= self.bstart_sk[spin, ikc]):
            raise ValueError("Invalid band index %d. should be in [%d, %d[" % (
                band, self.bstart_sk[spin, ikc], self.bstop_sk[spin, ikc]))

        return spin, ikc, band - self.bstart_sk[spin, ikc], self.sigma_kpoints[ikc]

    def sigkpt2index(self, sigma_kpoint):
        """
        Returns the index of the self-energy k-point in sigma_kpoints
        Used to access data in the arrays that are dimensioned [0:nkcalc]
        sigma_kpoint can be either an integer or list with reduced coordinates.
        """
        if duck.is_intlike(sigma_kpoint):
            ikc = int(sigma_kpoint)
            if self.nkcalc > ikc >= 0: return ikc
            raise ValueError("kpoint index should be in [0, %d] but received: %d" % (self.nkcalc, ikc))
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
        ikc = self.sigkpt2index(kpoint)
        bstart, bstop = self.bstart_sk[spin, ikc], self.bstop_sk[spin, ikc]

        return QpTempList([self.read_qp(spin, ikc, band, ignore_imag=ignore_imag)
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

        spin, ikc, ib, kpoint = self.get_sigma_skb_kpoint(spin, kpoint, band)

        # Abinit fortran (Ha units)
        # wrmesh_b(nwr, max_nbcalc, nkcalc, nsppol)
        # Frequency mesh along the real axis (Ha units) used for the different bands
        #print(spin, ikc, ib, self.read_variable("wrmesh_b").shape)
        wmesh = self.read_variable("wrmesh_b")[spin, ikc, ib, :] * units.Ha_to_eV

        # complex(dpc) :: vals_e0ks(ntemp, max_nbcalc, nkcalc, nsppol)
        # Sigma_eph(omega=eKS, kT, band)
        vals_e0ks = self.read_variable("vals_e0ks")[spin, ikc, ib, :, :] * units.Ha_to_eV
        vals_e0ks = vals_e0ks[:, 0] + 1j * vals_e0ks[:, 1]

        # complex(dpc) :: dvals_de0ks(ntemp, max_nbcalc, nkcalc, nsppol)
        # d Sigma_eph(omega, kT, band, kcalc, spin) / d omega (omega=eKS)
        dvals_de0ks = self.read_variable("dvals_de0ks")[spin, ikc, ib, :, :] * units.Ha_to_eV
        dvals_de0ks = dvals_de0ks[:, 0] + 1j * dvals_de0ks[:, 1]

        # real(dp) :: dw_vals(ntemp, max_nbcalc, nkcalc, nsppol)
        # Debye-Waller term (static).
        dw_vals = self.read_variable("dw_vals")[spin, ikc, ib, :] * units.Ha_to_eV

        # complex(dpc) :: vals_wr(nwr, ntemp, max_nbcalc, nkcalc, nsppol)
        # Sigma_eph(omega, kT, band) for given (k, spin).
        # Note: enk_KS corresponds to nwr/2 + 1.
        vals_wr = self.read_variable("vals_wr")[spin, ikc, ib, :, :, :] * units.Ha_to_eV
        vals_wr = vals_wr[:, :, 0] + 1j * vals_wr[:, :, 1]

        # Spectral function
        # nctkarr_t("spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol")
        spfunc_wr = self.read_variable("spfunc_wr")[spin, ikc, ib, :, :] / units.Ha_to_eV

        # Read QP data. Note band instead of ib index.
        qp = self.read_qp(spin, ikc, band)

        return EphSelfEnergy(wmesh, qp, vals_e0ks, dvals_de0ks, dw_vals, vals_wr, spfunc_wr)

    def read_qp(self, spin, kpoint, band, ignore_imag=False):
        """
        Return :class:`QpTempState` for the given (spin, kpoint, band).
        Only real part is returned if ``ignore_imag``.
        """
        spin, ikc, ibc, kpoint = self.get_sigma_skb_kpoint(spin, kpoint, band)

        def ri(a):
            return np.real(a) if ignore_imag else a

        # (Complex) QP energies computed with the dynamic formalism.
        # nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol")
        var = self.read_variable("qp_enes")
        qpe = (var[spin, ikc, ibc, :, 0] + 1j * var[spin, ikc, ibc, :, 1]) * units.Ha_to_eV

        # Adiabatic energies.
        # nctkarr_t("qpadb_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol")
        var = self.read_variable("qpadb_enes")
        qpe_adb = var[spin, ikc, ibc, :, 0] * units.Ha_to_eV

        # Debye-Waller term (static).
        # nctkarr_t("dw_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol"),
        var = self.read_variable("dw_vals")
        dw = var[spin, ikc, ibc, :] * units.Ha_to_eV

        # Sigma_eph(omega=eKS, kT, band, ikcalc, spin)
        # nctkarr_t("vals_e0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol")
        # TODO: Add Fan0 instead of computing Sigma - DW?
        var = self.read_variable("vals_e0ks")
        sigc = (var[spin, ikc, ibc, :, 0] + 1j * var[spin, ikc, ibc, :, 1]) * units.Ha_to_eV
        fan0 = sigc - dw

        # nctkarr_t("ks_enes", "dp", "max_nbcalc, nkcalc, nsppol")
        e0 = self.read_variable("ks_enes")[spin, ikc, ibc] * units.Ha_to_eV

        # nctkarr_t("ze0_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol")
        ze0 = self.read_variable("ze0_vals")[spin, ikc, ibc]

        return QpTempState(
            spin=spin,
            kpoint=kpoint,
            band=band,
            tmesh=self.tmesh,
            e0=e0,
            qpe=ri(qpe),
            ze0=ze0,
            fan0=ri(fan0),
            dw=dw,
            qpe_adb=qpe_adb,
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
            for ikc, kpoint in enumerate(self.sigma_kpoints):
                for band in range(self.bstart_sk[spin, ikc], self.bstop_sk[spin, ikc]):
                    qps.append(self.read_qp(spin, ikc, band, ignore_imag=ignore_imag))
            qps_spin[spin] = QpTempList(qps)

        return tuple(qps_spin)
