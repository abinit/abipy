"""
This module contains objects for postprocessing e-ph calculations
using the results stored in the GSTORE.nc file.
"""

from __future__ import annotations

import dataclasses
import itertools

# import abipy.core.abinit_units as abu
from functools import cached_property, lru_cache

import numpy as np
import pandas as pd
from monty.string import marquee  # , list_strings
from monty.termcolor import cprint

from abipy.abio.robots import Robot
from abipy.core.kpoints import kpoints_indices
from abipy.core.mixins import AbinitNcFile, Has_ElectronBands, Has_Header, Has_Structure  # , NotebookWriter
from abipy.core.structure import Structure

# from abipy.tools import duck
from abipy.electrons.ebands import ElectronBands, RobotWithEbands
from abipy.eph.common import BaseEphReader
from abipy.tools.numtools import BzRegularGridInterpolator, nparr_to_df
from abipy.tools.plotting import (
    add_fig_kwargs,
    get_ax_fig_plt,
    get_axarray_fig_plt,
    set_grid_legend,
)
from abipy.tools.typing import Figure, PathLike

GSTORE_KQ_MISSING = 0        # (k, q, spin) has not been computed.
GSTORE_KQ_COMPUTED = 1       # (k, q, spin) has been computed.
GSTORE_KQ_SYMMETRIZED = 2    # (k, q, spin) has been reconstructed by symmetry.


def _allclose(arr_name, array1, array2, verbose: int, rtol=1e-5, atol=1e-8) -> bool:
    """
    Wraps numpy allclose.
    """
    if np.allclose(array1, array2, rtol=rtol, atol=atol):
        if verbose:
            cprint(f"The arrays for {arr_name} are almost equal within the tolerances {rtol=}, {atol=}", color="green")
        return True

    if verbose:
        cprint(f"The arrays for {arr_name} are not almost equal within the tolerances {rtol=}, {atol=}", color="red")

    # differing_indices = np.where(~np.isclose(array1, array2, atol=atol))
    # for index in zip(*differing_indices):
    #    print(f"Difference at index {index}: array1 = {array1[index]}, array2 = {array2[index]}, difference = {abs(array1[index] - array2[index])}")
    return False


def _linreg_stats(x, y, wrap: bool = False) -> dict:
    """
    Least-squares linear-fit descriptors for a parity plot (y vs x): slope and intercept
    of y = slope*x + intercept, the Pearson correlation coefficient, and the RMSE.

    Args:
        wrap: True for angle data (in degrees), which is only defined modulo 360: realigns
            y to the branch closest to x (y -> x + wrap_to_180(y - x), which does not change
            y's physical value) before fitting, so a point straddling the +-180 branch cut
            isn't mistaken for a large disagreement.
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if wrap:
        y = x + (((y - x) + 180.0) % 360.0 - 180.0)

    if x.size > 1 and x.max() > x.min():
        slope, intercept = np.polyfit(x, y, 1)
        r = np.corrcoef(x, y)[0, 1]
    else:
        slope, intercept, r = np.nan, np.nan, np.nan

    rmse = np.sqrt(np.mean((x - y) ** 2)) if x.size > 0 else np.nan

    return dict(slope=slope, intercept=intercept, r=r, rmse=rmse, n=x.size, y=y)


def _add_fit_annotation(ax, x, y, xlim, wrap: bool = False, fontsize: int = 8, color: str = "C2") -> dict:
    """
    Draw the least-squares linear-fit line (see :func:`_linreg_stats`) across ``xlim`` and
    annotate ``ax`` with a text box reporting slope, intercept, Pearson r, RMSE, and the
    number of points. ``xlim`` is taken as an explicit argument (rather than ``ax.get_xlim()``)
    so the fit line spans the same range as the data regardless of plotting order/autoscaling.

    Returns the stats dict from :func:`_linreg_stats` (with the branch-realigned y, if
    ``wrap=True``) in case the caller wants to reuse it.
    """
    stats = _linreg_stats(x, y, wrap=wrap)

    if np.isfinite(stats["slope"]):
        xs = np.array(xlim, dtype=float)
        ax.plot(xs, stats["slope"] * xs + stats["intercept"], "-.", color=color, lw=1.0, zorder=4,
                label=f"fit: y = {stats['slope']:.3f}x {stats['intercept']:+.3g}")

    ax.text(0.03, 0.03,
            f"N = {stats['n']:,}\n"
            f"slope = {stats['slope']:.4f}\n"
            f"intercept = {stats['intercept']:.3g}\n"
            f"Pearson r = {stats['r']:.5f}\n"
            f"RMSE = {stats['rmse']:.3e}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=fontsize,
            bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.85), zorder=6)

    return stats


class GstoreFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands):
    """
    This file stores the e-ph matrix elements produced by the EPH code of Abinit
    and provides methods to analyze and plot results.

    Usage example:

    .. code-block:: python

        from abipy.eph.gstore import GstoreFile
        with GstoreFile("out_GSTORE.nc") as gstore:
            print(gstore)

            for spin in range(gstore.nsppol):
                # Extract the object storing the g for this spin.
                gqk = gstore.gqk_spin[spin]
                print(gqk)

                # Get a Dataframe with g(k, q) for all modes and bands.
                df = gqk.get_gdf_at_qpt_kpt([1/2, 0, 0], [0, 0, 0])
                print(df)

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreFile
    """

    @classmethod
    def from_file(cls, filepath: PathLike) -> GstoreFile:
        """Initialize the object from a netcdf file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        """
        Args:
            filepath: Path to the netcdf file.
        """
        super().__init__(filepath)
        self.r = GstoreReader(filepath)

    @cached_property
    def ebands(self) -> ElectronBands:
        """|ElectronBands| object."""
        return self.r.read_ebands()

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.ebands.structure

    def close(self) -> None:
        """Close the file."""
        self.r.close()

    @cached_property
    def gqk_spin(self) -> list:
        """List of |Gqk| objects, one for each spin."""
        return [Gqk.from_gstore(self, spin) for spin in range(self.nsppol)]

    @cached_property
    def has_gwpt(self) -> bool:
        """True if GSTORE contains GWPT matrix elements."""
        return self.r.gtype == "gwpt"

    @cached_property
    def params(self) -> dict:
        """Dict with the convergence parameters, e.g. ``nbsum``."""
        # od = OrderedDict([
        #    ("nbsum", self.nbsum),
        #    ("zcut", self.zcut),
        #    ("symsigma", self.symsigma),
        #    ("nqbz", self.r.nqbz),
        #    ("nqibz", self.r.nqibz),
        # ])
        ## Add EPH parameters.
        # od.update(self.r.common_eph_params)

        od = {}
        return od

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []
        app = lines.append

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))

        app("")
        app(self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"))
        if verbose > 1:
            app("")
            app(self.hdr.to_string(verbose=verbose, title="Abinit Header"))

        app(marquee("Gstore parameters", mark="="))
        app(f"nsppol: {self.r.nsppol}")
        app(f"gstore_completed: {bool(self.r.completed)}")
        app(f"kzone: {self.r.kzone}")
        app(f"kfilter: {self.r.kfilter}")
        app(f"gtype: {self.r.gtype}")
        app(f"qzone: {self.r.qzone}")
        app(f"with_vk: {self.r.with_vk}")
        app(f"kptopt: {self.r.kptopt}")
        app(f"qptopt: {self.r.qptopt}")
        app(f"use_lgk: {self.r.use_lgk}")
        app(f"use_lgq: {self.r.use_lgq}")

        for spin in range(self.r.nsppol):
            app(f"For {spin= }")
            app(f"\tbrange_k: {self.r.brange_k_spin[spin]}")
            app(f"\tbrange_kq: {self.r.brange_kq_spin[spin]}")
            app(f"\terange_spin: {self.r.erange_spin[spin]}")
            app(f"\tglob_spin_nq: {self.r.glob_spin_nq[spin]}")

        return "\n".join(lines)

    def check_unfilled_entries_in_gvals(self):
        """ """
        r = self.r
        for spin in range(self.nsppol):
            # nctkarr_t("gvals", "dp", "gstore_cplex, nb_kq, nb_k, natom3, glob_nk, glob_nq)
            variable = r.read_variable("gvals", path=f"gqk_spin{spin + 1}")
            fill_value = variable._FillValue
            # Read the data
            data = variable[:]
            missing_entries = np.where(data == fill_value)
            # Print the indices of missing entries
            print("Missing entries found at indices:", missing_entries)

            if self.r.kfilter == "none":
                raise ValueError("when kfilter == 'none' all the entries in gvals should have been written!")

    @lru_cache
    def get_gwpt_label_data(self, what: str, spin: int) -> tuple[str, np.array]:
        """Return label and numpy array to analyze according to `what`."""
        gqk = self.gqk_spin[spin]
        # TODO: Handle 1/0 and little group with huge values

        where = (np.abs(gqk.gvals_ks) > 1e-10) & (np.abs(gqk.gvals_ks) < 1e10)

        if what == "ratio":
            label = r"Ratio $|g^{\text{GWPT}}|/|g^{\text{KS}}|$"
            # data = np.abs(gqk.gvals) / np.abs(gqk.gvals_ks)

            data = np.ones_like(gqk.gvals, dtype=float)
            np.divide(
                np.abs(gqk.gvals),
                np.abs(gqk.gvals_ks),
                out=data,
                where=where,
            )

        elif what == "gwpt":
            label = r"$|g|^{\text{GWPT}}$"
            data = np.abs(gqk.gvals)

        elif what == "gks":
            label = r"$|g|^\text{KS}}$"
            data = np.abs(gqk.gvals_ks)

        else:
            raise ValueError(f"Invalid value for {what=}")

        return label, data

    @add_fig_kwargs
    def plot_gwpt_hist(
        self,
        what: str = "ratio",
        spin: int = 0,
        ax=None,
        hist_kwargs: dict | None = None,
        ratio_min: float = 0.0,
        ratio_max: float = 3.0,
        **kwargs,
    ) -> Figure:
        """
        Plot histograms of the GWPT and KS e-ph matrix elements and of their ratio.

        Args:
            what:
            spin: spin index
            ax: |matplotlib-Axes| or None if a new figure should be created.
            hist_kwargs:
            ratio_min: lower bound used to clip the ratio |g^GWPT|/|g^KS|
                before histogramming. Defaults to 0.0.
            ratio_max: upper bound used to clip the ratio |g^GWPT|/|g^KS|
                before histogramming. Defaults to 3.0. Values outside
                [ratio_min, ratio_max] are dropped (they are usually numerical
                artifacts from |g^KS| being close to zero).
            fontsize: legend and label fontsize.
        """
        if not self.has_gwpt:
            raise ValueError("GSTORE does not contain GWPT matrix elements.")

        if ratio_min >= ratio_max:
            raise ValueError(f"ratio_min ({ratio_min}) must be < ratio_max ({ratio_max}).")

        what_list = ("gwpt", "gks", "ratio")
        xlabels = {
            "gwpt": r"$|g^{GW}|$",
            "gks": r"$|g^{KS}|$",
            "ratio": r"Ratio $|g^{GW}|/|g^{KS}|$",
        }

        ax_list = None
        ax_list, fig, plt = get_axarray_fig_plt(
            ax_list, nrows=len(what_list), ncols=1, sharex=False, sharey=True, squeeze=True
        )

        for what, ax in zip(what_list, ax_list, strict=False):
            _, data = self.get_gwpt_label_data(what, spin)
            data_flat = data.flatten()

            if what == "ratio":
                # Drop non-finite values and values outside [ratio_min, ratio_max].
                # The ratio can blow up when |g^KS| is close to zero, which would
                # otherwise stretch the histogram x-axis and hide the bulk of the
                # distribution (typically centered near 1).
                finite = np.isfinite(data_flat)
                in_range = (data_flat >= ratio_min) & (data_flat <= ratio_max)
                data_flat = data_flat[finite & in_range]

            hist_kwargs_ = hist_kwargs or {}
            ax.hist(data_flat, **hist_kwargs_)
            ax.set_xlabel(xlabels[what])
            ax.set_ylabel("Count")

            if what == "ratio":
                ax.set_xlim(ratio_min, ratio_max)
            else:
                # |g| is non-negative; anchor the left edge at 0 so the three
                # subplots align visually at x=0 (matplotlib's auto-scale
                # otherwise pads the left side with a small negative margin).
                ax.set_xlim(left=0)

        # Add vertical spacing so the xlabel of each subplot is not hidden
        # behind the next subplot's frame.
        fig.subplots_adjust(hspace=0.45)

        return fig

    @add_fig_kwargs
    def plot_gwpt_heatmap_kq(self, kpoint, qpoint, what: str = "ratio", spin: int = 0, ax=None, **kwargs) -> Figure:
        """
        Plot heatmap with the ratio between the GWPT and the KS e-ph matrix elements.

        Args:
            kpoint:
            qpoint:
            what:
            spin: spin index
            ax: |matplotlib-Axes| or None if a new figure should be created.
            hist_kwargs:
            fontsize: legend and label fontsize.
        """
        if not self.has_gwpt:
            raise ValueError("GSTORE does not contain GWPT matrix elements.")

        ik_glob, kpoint = self.r.find_ik_glob_kpoint(kpoint, spin)
        iq_glob, qpoint = self.r.find_iq_glob_qpoint(qpoint, spin)

        label, data = self.get_gwpt_label_data(what, spin)

        natom = len(self.structure)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=natom, ncols=3, sharex=True, sharey=True, squeeze=False)

        for iat, idir in itertools.product(range(natom), range(3)):
            # (glob_nq, glob_nk, natom3, nb_kq, nb_k)
            ipc = idir + iat * 3
            grid_mn = data[iq_glob, ik_glob, ipc]

            ax = ax_mat[iat, idir]
            ax.imshow(grid_mn, aspect="auto", origin="lower")
            # ax.set_colorbar(label=label)
            ax.set_xlabel("m band index")
            ax.set_ylabel("n band index")
            # ax.set_title("A/B averaged over bands")

        return fig

    @add_fig_kwargs
    def plot_gwpt_vs_qpts(
        self,
        kpoint,
        band_k: int,
        frac_bounds,
        band_kq_range: range | list,
        spin: int = 0,
        dist_tol: float = 1e-12,
        fontsize: int = 8,
        **kwargs,
    ) -> Figure:
        """ """
        if not self.has_gwpt:
            raise ValueError("GSTORE does not contain GWPT matrix elements.")

        gqk = self.gqk_spin[spin]
        in_k = band_k - gqk.bstart_k
        ik_glob, kpoint = self.r.find_ik_glob_kpoint(kpoint, spin)

        reciprocal_lattice = self.structure.reciprocal_lattice
        cart_bounds = reciprocal_lattice.get_cartesian_coords(frac_bounds)

        qpt_cart_coords = []
        for iq in range(gqk.glob_nq):
            iq_bz = self.r.qglob2bz[spin, iq]
            qpt_cart_coords.append(reciprocal_lattice.get_cartesian_coords(self.r.qbz[iq_bz]))

        from abipy.core.kpoints import find_points_along_path

        p = find_points_along_path(cart_bounds, qpt_cart_coords, dist_tol)
        if not len(p.ikfound):
            raise ValueError(
                f"Cannot find q-points with {dist_tol=}. Check input boundaries or try to increase dist_tol."
            )
        # p.ikfound, p.dist_list, p.path_ticks)

        natom = len(self.structure)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=natom, ncols=3, sharex=False, sharey=False, squeeze=False)

        # shape: (glob_nq, natom3, nb_kq)
        g_qpm = gqk.gvals[:, ik_glob, :, :, in_k]
        gks_qpm = gqk.gvals_ks[:, ik_glob, :, :, in_k]

        xs = list(range(len(p.ikfound)))
        for iat, idir in itertools.product(range(natom), range(3)):
            ax = ax_mat[iat, idir]
            ipert = iat + idir
            # (glob_nq, nb_kq) --> (nb_kq, glob_nq).
            gwpt_bq, ks_bq = g_qpm[:, ipert, :].T.copy(), gks_qpm[:, ipert, :].T.copy()

            gwpt_style = dict(marker="o", color="red")
            ks_style = dict(marker=".", color="blue")
            ratio_style = dict(ls="--", color="k")
            # band_kq_range = [0, 4]

            for ib_kq in range(gqk.nb_kq):
                if band_kq_range is not None and ib_kq not in band_kq_range:
                    continue
                gwpt_ys = np.abs(gwpt_bq[ib_kq, p.ikfound])
                ks_ys = np.abs(ks_bq[ib_kq, p.ikfound])
                ax.plot(xs, gwpt_ys, **gwpt_style)
                ax.plot(xs, ks_ys, **ks_style)

                # Plot enhancement ratio.
                ratio_ax = ax.twinx()
                ratio = gwpt_ys / ks_ys
                ratio_ax.plot(xs, ratio, **ratio_style)

            # set_grid_legend(ax, fontsize, xlabel=r"band index (kq)")

        return fig

    @add_fig_kwargs
    def plot_gwpt_vs_ks(
        self, kpoint, band_k: int, spin: int = 0, fontsize: int = 8, colormap: str = "jet", **kwargs
    ) -> Figure:
        """ """
        if not self.has_gwpt:
            raise ValueError("GSTORE does not contain GWPT matrix elements.")

        natom = len(self.structure)
        ax_mat, fig, plt = get_axarray_fig_plt(None, nrows=natom, ncols=3, sharex=False, sharey=False, squeeze=False)

        gqk = self.gqk_spin[spin]
        in_k = band_k - gqk.bstart_k
        ik_glob, kpoint = self.r.find_ik_glob_kpoint(kpoint, spin)

        # shape: (glob_nq, natom3, nb_kq)
        g_qpm = gqk.gvals[:, ik_glob, :, :, in_k]
        gks_qpm = gqk.gvals_ks[:, ik_glob, :, :, in_k]

        cmap = plt.get_cmap(colormap)
        colors = [cmap(iq / gqk.glob_nq) for iq in range(gqk.glob_nq)]
        xs = list(range(gqk.nb_kq))

        for iat, idir in itertools.product(range(natom), range(3)):
            ax = ax_mat[iat, idir]
            ipert = iat + idir
            # (glob_nq, nb_kq)
            gwpt_ys, ks_ys = g_qpm[:, ipert, :], gks_qpm[:, ipert, :]

            for iq in range(gqk.glob_nq):
                numerator = np.abs(gwpt_ys[iq])
                denominator = np.abs(ks_ys[iq])
                # ratio = numerator / denominator
                # ratio = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=np.abs(denominator) > 1e-2)
                # ax.plot(xs, ratio, color=colors[iq])
                ax.plot(xs, np.abs(gwpt_ys[iq]), color=colors[iq], ls="-")
                ax.plot(xs, np.abs(ks_ys[iq]), color=colors[iq], ls="--")

            set_grid_legend(ax, fontsize, xlabel=r"band index (kq)")

        return fig

    @add_fig_kwargs
    def plot_gwpt_vs_ks_scatter(
        self,
        spin: int = 0,
        ratio_min: float | None = None,
        ratio_max: float | None = None,
        ks_tol: float = 1e-7,
        fit_intercept: bool = False,
        colormap: str = "viridis",
        zoom_factor: float = 2.0,
        with_inset: bool = True,
        inset_loc: str = "lower right",
        scatter_kwargs: dict | None = None,
        ax=None,
        fontsize: int = 8,
        **kwargs,
    ) -> Figure:
        """
        Scatter plot of |g^GW| vs |g^KS| over all matrix elements for a given spin.

        Each point is one (q, k, band_kq, perturbation, band_k) tuple. Points are
        colored by the ratio |g^GW|/|g^KS|, and a linear least-squares fit is
        overlaid. An optional inset shows the full data range so outliers are
        visible without disturbing the zoomed main axes.

        Args:
            spin: spin index.
            ratio_min: optional lower bound of the ratio window. Points with
                ratio < ratio_min are dropped. If None (default), no lower
                bound is applied.
            ratio_max: optional upper bound of the ratio window. Points with
                ratio > ratio_max are dropped. If None (default), no upper
                bound is applied.
            ks_tol: |g^KS| values <= ks_tol are discarded to avoid huge
                ratios from a near-zero denominator. Default 1e-7.
            fit_intercept: if True, fit y = a*x + b; otherwise fit y = a*x
                (forced through the origin).
            colormap: matplotlib colormap used to color points by ratio.
                When ratio_min/ratio_max are None the colormap is clipped to
                the 1st/99th percentile of the ratio distribution so a few
                outliers don't wash out the rest of the colors.
            zoom_factor: the main axes go from 0 to
                zoom_factor * max(median(|g^KS|), median(|g^GW|)) on both axes.
            with_inset: if True, add an inset showing the full data range
                with a rectangle indicating the zoom window.
            inset_loc: matplotlib location string for the inset.
            scatter_kwargs: extra kwargs forwarded to ``ax.scatter``
                (e.g. ``{"s": 6, "alpha": 0.6}``).
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: legend fontsize.
        """
        if not self.has_gwpt:
            raise ValueError("GSTORE does not contain GWPT matrix elements.")
        if (ratio_min is not None and ratio_max is not None
                and ratio_min >= ratio_max):
            raise ValueError(f"ratio_min ({ratio_min}) must be < ratio_max ({ratio_max}).")

        gqk = self.gqk_spin[spin]
        g_gw = np.abs(np.asarray(gqk.gvals)).ravel()
        g_ks = np.abs(np.asarray(gqk.gvals_ks)).ravel()

        # Keep finite values with |g^KS| above the safety threshold.
        # ks_tol protects against numerically pathological ratios where the
        # denominator is essentially zero; it is *not* a user-facing
        # visualization filter (use ratio_min/ratio_max for that).
        valid = np.isfinite(g_gw) & np.isfinite(g_ks) & (g_ks > ks_tol)
        x = g_ks[valid]
        y = g_gw[valid]
        ratio = y / x

        # Optional ratio window. Sides default to ±inf (no filtering).
        lo = -np.inf if ratio_min is None else ratio_min
        hi = +np.inf if ratio_max is None else ratio_max
        if not (np.isneginf(lo) and np.isposinf(hi)):
            window = (ratio >= lo) & (ratio <= hi)
            x, y, ratio = x[window], y[window], ratio[window]

        if x.size < 2:
            raise RuntimeError(
                f"Only {x.size} point(s) remain after filtering; cannot fit."
            )

        # Linear fit.
        if fit_intercept:
            slope, intercept = np.polyfit(x, y, 1)
            fit_label = f"fit: y = {slope:.4g} x + {intercept:.4g}"
        else:
            # Closed-form least squares for y = a*x: a = (x.y) / (x.x).
            slope = float(np.dot(x, y) / np.dot(x, x))
            intercept = 0.0
            fit_label = f"fit: y = {slope:.4g} x"

        ax, fig, _ = get_ax_fig_plt(ax=ax)

        scatter_kwargs_ = {"s": 6, "alpha": 0.6, "edgecolors": "none"}
        if scatter_kwargs:
            scatter_kwargs_.update(scatter_kwargs)

        # Colormap range: user-supplied bounds win; otherwise clip to the
        # 1st/99th percentile so a handful of outliers don't wash out the
        # rest of the color resolution.
        if ratio_min is not None:
            cmap_vmin = ratio_min
        else:
            cmap_vmin = float(np.percentile(ratio, 1))
        if ratio_max is not None:
            cmap_vmax = ratio_max
        else:
            cmap_vmax = float(np.percentile(ratio, 99))

        sc = ax.scatter(
            x, y, c=ratio, cmap=colormap, vmin=cmap_vmin, vmax=cmap_vmax,
            **scatter_kwargs_,
        )
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label(r"$|g^{GW}|/|g^{KS}|$")

        # Zoomed main axes derived from medians (robust to outliers).
        med_x = float(np.median(x))
        med_y = float(np.median(y))
        main_max = zoom_factor * max(med_x, med_y)

        x_line = np.linspace(0.0, main_max, 200)
        ax.plot(x_line, slope * x_line + intercept,
                color="red", linewidth=2.0, label=fit_label)

        ax.set_xlabel(r"$|g^{KS}|$")
        ax.set_ylabel(r"$|g^{GW}|$")
        ax.set_xlim(0.0, main_max)
        ax.set_ylim(0.0, main_max)
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="upper left", fontsize=fontsize)

        if with_inset:
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            from matplotlib.patches import Rectangle

            full_max = max(float(x.max()), float(y.max())) * 1.02

            axins = inset_axes(
                ax, width="40%", height="40%", loc=inset_loc, borderpad=1.5,
            )
            axins.scatter(
                x, y, c=ratio, cmap=colormap,
                vmin=cmap_vmin, vmax=cmap_vmax, **scatter_kwargs_,
            )
            x_line_f = np.linspace(0.0, full_max, 200)
            axins.plot(x_line_f, slope * x_line_f + intercept,
                       color="red", linewidth=1.5)
            axins.set_xlim(0.0, full_max)
            axins.set_ylim(0.0, full_max)
            axins.set_aspect("equal", adjustable="box")
            axins.set_title("full range", fontsize=fontsize + 1)
            axins.tick_params(axis="both", labelsize=fontsize)

            # Rectangle on the inset outlining the zoom window of the main axes.
            if full_max > main_max:
                axins.add_patch(Rectangle(
                    (0.0, 0.0), main_max, main_max,
                    fill=False, edgecolor="gray", linewidth=0.8,
                ))

        return fig


@dataclasses.dataclass(kw_only=True)
class Gqk:
    """
    This object stores the e-ph matrix elements (g or g^2) and the matrix elements
    of the velocity operator for a given spin.
    """

    spin: int  # Spin index.
    nb_k: int  # Number of bands at k.
    nb_kq: int  # Number of bands at k+q.
    bstart_k: int  # Initial band at k.
    bstart_kq: int  # Initial band at k+q.
    glob_nk: int  # Total number of k/q points in global matrix.
    glob_nq: int  # Note that k-points/q-points can be filtered.

    gstore: GstoreFile

    gvals: np.ndarray  # Array of shape (glob_nq, glob_nk, natom3, nb_kq, nb_k)
    # storing complex g(k,q) in the atom representation.

    gvals_ks: np.ndarray | None  # Same as gvals but for KS if we are in GWPT mode.

    vk_cart_ibz: np.ndarray | None
    vkmat_cart_ibz: np.ndarray | None

    def __post_init__(self):
        """Implement consistency check"""
        natom3 = len(self.structure) * 3
        expected_shape = (self.glob_nq, self.glob_nk, natom3, self.nb_kq, self.nb_k)
        if self.gvals.shape != expected_shape:
            raise ValueError(f"{self.gvals.shape=} != {expected_shape=}")
        if self.gvals_ks is not None and self.gvals_ks.shape != expected_shape:
            raise ValueError(f"{self.gvals_ks.shape=} != {expected_shape=}")

    @classmethod
    def from_gstore(cls, gstore: GstoreFile, spin: int) -> Gqk:
        """
        Build an instance from a GstoreFile and the spin index.
        """
        ncr = gstore.r
        path = f"gqk_spin{spin + 1}"
        nb_k = ncr.read_dimvalue("nb_k", path=path)
        nb_kq = ncr.read_dimvalue("nb_kq", path=path)
        glob_nk = ncr.read_dimvalue("glob_nk", path=path)
        glob_nq = ncr.read_dimvalue("glob_nq", path=path)

        # Read e-ph matrix elements
        # nctkarr_t("gvals", "dp", "gstore_cplex, nb_kq, nb_k, natom3, glob_nk, glob_nq)
        # Have to transpose the (nb_kq, nb_k) submatrix written by Fortran.
        # Remember that gvals on disk are always complex, in Hartree units and in the atomic representation.

        gvals = ncr.read_value("gvals", path=path).transpose(0, 1, 2, 4, 3, 5).copy()
        gvals = gvals[..., 0] + 1j * gvals[..., 1]

        # Try to read KS gvals (produced by GWPT code)
        gvals_ks = None
        if gstore.has_gwpt:
            gvals_ks = ncr.read_value("gvals_ks", path=path).transpose(0, 1, 2, 4, 3, 5).copy()
            gvals_ks = gvals_ks[..., 0] + 1j * gvals_ks[..., 1]

        vk_cart_ibz, vkmat_cart_ibz = None, None
        if ncr.with_vk == 1:
            # nctk_def_arrays(spin_ncid, nctkarr_t("vk_cart_ibz", "dp", "three, nb_k, gstore_nkibz"))
            vk_cart_ibz = ncr.read_value("vk_cart_ibz", path=path)

        if ncr.with_vk == 2:
            # Full (nb_k x nb_k) matrix.
            # Have to transpose (nb_kq, nb_k) submatrix written by Fortran.
            # nctk_def_arrays(spin_ncid, nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb_k, nb_k, gstore_nkibz"))
            vkmat_cart_ibz = ncr.read_value("vkmat_cart_ibz", path=path).transpose(0, 1, 3, 2, 4).copy()
            vkmat_cart_ibz = vkmat_cart_ibz[..., 0] + 1j * vkmat_cart_ibz[..., 1]

        # Note conversion between Fortran and python indexing.
        bstart_k = ncr.read_value("bstart_k", path=path) - 1
        bstart_kq = ncr.read_value("bstart_kq", path=path) - 1

        data = locals()
        return cls(**{k: data[k] for k in [field.name for field in dataclasses.fields(Gqk)]})

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosiy level ``verbose``."""
        lines = []
        app = lines.append

        app(marquee(f"Gqk for spin: {self.spin}", mark="="))
        app(f"bstart_k: {self.bstart_k}")
        app(f"nb_k: {self.nb_k}")
        app(f"bstart_kq: {self.bstart_kq}")
        app(f"nb_kq: {self.nb_kq}")
        app(f"glob_nk: {self.glob_nk}")
        app(f"glob_nq: {self.glob_nq}")

        return "\n".join(lines)

    @property
    def structure(self) -> Structure:
        """|Structure| object."""
        return self.gstore.structure

    @cached_property
    def g2(self) -> np.ndarray:
        """g2 in the atomic representation in Ha^2."""
        return np.abs(self.gvals) ** 2

    # @cached_property
    # def g_nu(self) -> np.ndarray:
    #   """g in the phonon representation in Ha."""
    #   self.gvals

    @cached_property
    def g2_ks(self) -> np.ndarray | None:
        """KS g2 in the atomic representation in Ha^2."""
        if self.gvals_ks is None:
            return None
        return np.abs(self.gvals_ks) ** 2

    # @cached_property
    # def g_ks_nu(self) -> np.ndarray:
    #   """g in the phonon representation in Ha."""
    #   if self.gvals_ks is None:
    #       return None

    def get_dataframe(self, what: str = "g2") -> pd.DataFrame:
        """
        Build and return a dataframe with all the |g(k,q)|^2 if what == "g2" or
        all |v_nk|^2 if what == "v2".
        """
        if what == "g2":
            df = nparr_to_df("g2", self.g2, ["iq", "ik", "imode", "m_kq", "n_k"])

        elif what == "v2":
            if self.vk_cart_ibz is None:
                raise ValueError("vk_cart_ibz is not available in GSTORE!")
            # Compute the squared norm of each vector
            v2 = np.sum(self.vk_cart_ibz**2, axis=2)
            df = nparr_to_df("v2", v2, ["ik", "n_k"])

        else:
            raise ValueError(f"Invalid {what=}")

        # Shift band indices.
        df["m_kq"] += self.bstart_kq
        df["n_k"] += self.bstart_k

        return df

    def get_g2q_interpolator_kpoint(self, kpoint, method="linear", check_mesh=1) -> BzRegularGridInterpolator:
        r"""
        Build and return an interpolator that can be used to interpolate g^2(q)

        NB: Invoking the interpolation with an arbitrary q-point returns a numpy array
        of shape (nb_kq, nb_k, natom3) with g_{m_kq n_k, \nu}(q)
        """
        r = self.gstore.r

        # Find the index of the kpoint.
        ik_g, kpoint = r.find_ik_glob_kpoint(kpoint, self.spin)

        # Compute indices of qpoints in the ngqpt mesh.
        ngqpt, shifts = r.ngqpt, [0, 0, 0]
        q_indices = kpoints_indices(r.qbz, ngqpt, shifts, check_mesh=check_mesh)

        natom3 = 3 * len(self.structure)
        nb_k = self.nb_k
        nb_kq = self.nb_kq
        nx, ny, nz = ngqpt
        assert nb_k == nb_kq

        # (glob_nq, glob_nk, natom3, m_kq, n_k)
        g2_qph_mn = self.g2[:, ik_g]

        # Insert g2 in g2_grid
        g2_grid = np.empty((nb_k, nb_kq, natom3, nx, ny, nz))
        for nu in range(natom3):
            for g2_mn, q_inds in zip(g2_qph_mn[:, nu], q_indices, strict=False):
                ix, iy, iz = q_inds
                g2_grid[:, :, nu, ix, iy, iz] = g2_mn

        return BzRegularGridInterpolator(self.structure, shifts, g2_grid, method=method)

    def get_g_qpt_kpt(self, qpoint, kpoint, what) -> np.ndarray:
        """
        Return numpy array with e-ph matrix elements for the given (qpoint, kpoint) pair.

        Args:
            what="g2" for |g(k,q)|^2, "g" for g(k,q)
        """
        # Find the internal indices of (qpoint, kpoint)
        iq_g, qpoint = self.gstore.r.find_iq_glob_qpoint(qpoint, self.spin)
        ik_g, kpoint = self.gstore.r.find_ik_glob_kpoint(kpoint, self.spin)

        if what == "g2":
            return self.g2[iq_g, ik_g]
        if what == "g":
            return self.gvals[iq_g, ik_g]

        raise ValueError(f"Invalid {what=}")

    def get_gdf_at_qpt_kpt(self, qpoint, kpoint, what="g2") -> pd.DataFrame:
        """
        Build and return a dataframe with the |g(k,q)|^2 for the given (qpoint, kpoint) pair.

        Args:
            what="g2" for |g(k,q)|^2, "g" for g(k,q)
        """
        g2_slice = self.get_g_qpt_kpt(qpoint, kpoint, what)
        df = nparr_to_df(what, g2_slice, ["imode", "m_kq", "n_k"])

        # Shift band indices.
        df["m_kq"] += self.bstart_kq
        df["n_k"] += self.bstart_k

        return df

    def neq(self, other: Gqk, verbose: int) -> int:
        """
        Helper function to compare two GQK objects.
        """
        # This dimensions must agree in order to have a meaningful comparison.
        # so raise immediately if not equal.
        aname_list = ["spin", "nb_k", "nb_kq", "glob_nk", "glob_nq"]

        for aname in aname_list:
            val1, val2 = getattr(self, aname), getattr(other, aname)

            if isinstance(val1, (str, int, float)):
                eq = val1 == val2
            elif isinstance(val1, np.ndarray):
                eq = np.allclose(val1, val2)
            else:
                raise TypeError(f"Don't know how to handle comparison for type: {type(val1)}")

            if not eq:
                raise RuntimeError(f"Different values of {aname=}, {val1=}, {val2=}")

        ierr = 0
        kws = dict(verbose=verbose)  # , atol= rtol)

        # Compare v_nk or v_mn_k.
        if self.vk_cart_ibz is not None:
            if not _allclose("vk_cart_ibz", self.vk_cart_ibz, other.vk_cart_ibz, **kws):
                ierr += 1

        if self.vkmat_cart_ibz is not None:
            if not _allclose("vkmat_cart_ibz", self.vkmat_cart_ibz, other.vkmat_cart_ibz, **kws):
                ierr += 1

        # Compare g or g^2.
        if not _allclose("g2", self.g2, other.g2, **kws):
            ierr += 1

        if self.gvals is not None:
            if not _allclose("gvals", self.gvals, other.gvals, **kws):
                ierr += 1

        return ierr

    # @add_fig_kwargs
    # def plot_g2_hist(self, ax_list=None, **kwargs) -> Figure:

    #    natom = len(self.structure)
    #    nrows, ncols, gridspec_kw = natom, 3, None
    #    ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
    #                                           sharex=True, sharey=True, squeeze=False, gridspec_kw=gridspec_kw)
    #    ax_list = ax_list.ravel()

    #    # (glob_nq, glob_nk, natom3, m_kq, n_k)
    #    for imode, ax in zip(range(natom * 3), ax_list):
    #        data = self.g2[:,:,imode,:,:].flatten()
    #        ax.hist(data)

    #    return fig


class GstoreReader(BaseEphReader):
    """
    Reads data from file and constructs objects.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreReader
    """

    def __init__(self, filepath: PathLike):
        """
        Args:
            filepath: Path to the netcdf file.
        """
        super().__init__(filepath)

        # Read important dimensions.
        self.nsppol = self.read_dimvalue("number_of_spins")
        self.nkbz = self.read_dimvalue("gstore_nkbz")
        self.nkibz = self.read_dimvalue("gstore_nkibz")
        self.nqbz = self.read_dimvalue("gstore_nqbz")
        self.nqibz = self.read_dimvalue("gstore_nqibz")

        # Read important variables.
        self.completed = self.read_value("gstore_completed")
        self.done_spin_qbz = self.read_value("gstore_done_qbz_spin")
        self.with_vk = self.read_value("gstore_with_vk")
        self.qptopt = self.read_value("gstore_qptopt")
        self.kptopt = self.read_value("kptopt")
        self.kzone = self.read_string("gstore_kzone")
        self.qzone = self.read_string("gstore_qzone")
        self.kfilter = self.read_string("gstore_kfilter")
        self.gtype = self.read_string("gstore_gtype")

        # Note conversion Fortran --> C for the isym index.
        self.brange_k_spin = self.read_value("gstore_brange_k_spin")
        self.brange_k_spin[:, 0] -= 1
        self.brange_kq_spin = self.read_value("gstore_brange_kq_spin")
        self.brange_kq_spin[:, 0] -= 1

        self.erange_spin = self.read_value("gstore_erange_spin")
        # Total number of k/q points for each spin after filtering (if any)
        self.glob_spin_nq = self.read_value("gstore_glob_nq_spin")
        self.glob_nk_spin = self.read_value("gstore_glob_nk_spin")

        # K-points and q-points in the IBZ
        self.kibz = self.read_value("reduced_coordinates_of_kpoints")
        self.qibz = self.read_value("gstore_qibz")

        # K-points and q-points in the BZ
        self.kbz = self.read_value("gstore_kbz")
        self.qbz = self.read_value("gstore_qbz")
        self.ngqpt = self.read_value("gstore_ngqpt")

        # Mapping BZ --> IBZ. Note conversion Fortran --> C for the isym index.
        # nctkarr_t("gstore_kbz2ibz", "i", "six, gstore_nkbz"), &
        # nctkarr_t("gstore_qbz2ibz", "i", "six, gstore_nqbz"), &
        self.kbz2ibz = self.read_value("gstore_kbz2ibz")
        self.kbz2ibz[:, 0] -= 1

        self.qbz2ibz = self.read_value("gstore_qbz2ibz")
        self.qbz2ibz[:, 0] -= 1

        # Mapping q/k points in gqk --> BZ. Note conversion Fortran --> C for indexing.
        # nctkarr_t("gstore_qglob2bz", "i", "gstore_max_nq, number_of_spins"), &
        # nctkarr_t("gstore_kglob2bz", "i", "gstore_max_nk, number_of_spins") &
        self.qglob2bz = self.read_value("gstore_qglob2bz")
        self.qglob2bz -= 1
        self.kglob2bz = self.read_value("gstore_kglob2bz")
        self.kglob2bz -= 1

        self.use_lgk = self.read_value("gstore_use_lgk")
        self.use_lgq = self.read_value("gstore_use_lgq")

    def find_iq_glob_qpoint(self, qpoint, spin: int) -> tuple:
        """
        Find the internal index of the qpoint needed to access the gvals array.
        Return index of the q-point and qpoint as array.
        """
        qpoint = np.asarray(qpoint)
        for iq_g, iq_bz in enumerate(self.qglob2bz[spin]):
            if np.allclose(qpoint, self.qbz[iq_bz]):
                # print(f"Found {qpoint = } with index {iq_g = }")
                return iq_g, qpoint

        raise ValueError(f"Cannot find {qpoint=} in {self.path}")

    def find_ik_glob_kpoint(self, kpoint, spin: int) -> tuple:
        """
        Find the internal indices of the kpoint needed to access the gvals array.
        Return index of the k-point and kpoint as array.
        """
        kpoint = np.asarray(kpoint)
        for ik_g, ik_bz in enumerate(self.kglob2bz[spin]):
            if np.allclose(kpoint, self.kbz[ik_bz]):
                # print(f"Found {kpoint = } with index {ik_g = }")
                return ik_g, kpoint

        raise ValueError(f"Cannot find {kpoint=} in {self.path}")

    # TODO: This fix to read groups should be imported in pymatgen.
    @cached_property
    def path2group(self) -> dict:
        """Dictionary mapping path to group."""
        return self.rootgrp.groups


class GstoreRobot(Robot, RobotWithEbands):
    """
    This robot analyzes the results contained in multiple GSTORE.nc files.

    Usage example:

    .. code-block:: python

        robot = GstoreRobot.from_files([
            "t04o_GSTORE.nc",
            "t05o_GSTORE.nc",
            ])


    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GstoreRobot
    """

    EXT = "GSTORE"

    def neq(self, ref_basename: str | None = None, verbose: int = 0) -> int:
        """
        Compare all GSTORE.nc files stored in the GstoreRobot
        """
        # Find reference gstore. By default the first file in the robot is used.
        ref_gstore = self._get_ref_abifile_from_basename(ref_basename)

        exc_list = []
        ierr = 0
        for other_gstore in self.abifiles:
            if ref_gstore.filepath == other_gstore.filepath:
                continue
            print("Comparing: ", ref_gstore.basename, " with: ", other_gstore.basename)
            try:
                ierr += self._neq_two_gstores(ref_gstore, other_gstore, verbose)
                cprint("EQUAL", color="green")
            except Exception as exc:
                exc_list.append(str(exc))

        for exc in exc_list:
            cprint(exc, color="red")

        return ierr

    def _neq_two_gstores(self: GstoreRobot, gstore2: GstoreRobot, verbose: int) -> int:
        """
        Helper function to compare two GSTORE files.
        """
        # These quantities must be the same to have a meaningful comparison.
        aname_list = [
            "structure",
            "nsppol",
            "nkbz",
            "nkibz",
            "nqbz",
            "nqibz",
            "completed",
            "kzone",
            "qzone",
            "kfilter",
            "brange_k_spin",
            "brange_kq_spin",
            "erange_spin",
            "glob_spin_nq",
            "glob_nk_spin",
        ]

        for aname in aname_list:
            self._compare_attr_name(aname, self, gstore2)

        # Now compare the gkq objects for each spin.
        ierr = 0
        for spin in range(self.nsppol):
            gqk1, gqk2 = self.gqk_spin[spin], gstore2.gqk_spin[spin]
            ierr += gqk1.neq(gqk2, verbose)

        return ierr

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """
        # for fig in self.get_ebands_plotter().yield_figs(): yield fig

    def write_notebook(self, nbpath=None) -> str:
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend(
            [
                # nbv.new_markdown_cell("# This is a markdown cell"),
                nbv.new_code_cell("robot = abilab.GstoreRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
                # nbv.new_code_cell("ebands_plotter = robot.get_ebands_plotter()"),
            ]
        )

        # Mixins
        # nb.cells.extend(self.get_baserobot_code_cells())
        # nb.cells.extend(self.get_ebands_code_cells())

        return self._write_nb_nbpath(nb, nbpath)

    @add_fig_kwargs
    def compare_gvals_with_reconstruction(self, ax_list=None, tol: float = 1e-6, mag_tol: float = 1e-4,
                                           n_outliers: int = 1, use_hexbin: bool = False, bins: int = 50,
                                           verbose: int = 0, fontsize: int = 8,
                                           **kwargs) -> Figure:
        """
        This is a debugging tool.
        It reads gvals from two GSTORE files. In one case, the e-ph matrix elements
        have been computed in the full BZ both for k and q. In the other
        case, the gvals have been computed with k in the IBZ and q in the IBZ_k.
        and then reconstructed by symmetry at the end of the run.
        The goal of this method is to check that both calculations produce the same results.
        by printing some statistics and producing, for each spin, a 2x2 grid of parity plots for
        the e-ph matrix elements: |g| (top-left), phase of g in degrees (top-right), Re(g)
        (bottom-left), and Im(g) (bottom-right). A single complex matrix element can be wrong in
        different ways -- only in magnitude, only in phase, or in a way that's obvious in the
        Cartesian (Re, Im) view but not in the polar (|g|, phase) one or vice versa (e.g. a small
        |g| with a large phase error can still look fine in Re/Im if both components happen to
        be small) -- so looking at |g1 - g2| alone, as the original version of this method did,
        cannot distinguish these failure modes; all four views are needed together.

        Only the (k, q) points that were actually reconstructed by symmetry in the first file
        (gstore_glob_state_kqs == GSTORE_KQ_SYMMETRIZED) are used for the statistics and the plots.
        Including the directly-computed points (trivially expected to agree in both files) and
        flattening the whole (glob_nq, glob_nk, natom3, nb_kq, nb_k) array, as done previously,
        buries the actual test: e.g. (k, q) points that reconstruct exactly (diff ~1e-15) get
        averaged together with (k, q) points that are completely wrong (diff ~O(1)) into a single
        "mean |g1-g2|" scalar, and the parity plot is dominated by the many small-magnitude,
        off-diagonal-like matrix elements shared by both categories, so a ~30/70 exact/wrong split
        can look like an innocuous scatter around y = x.

        All four plots share the same matched/mismatched classification, based on the full complex
        difference (magnitude AND phase), so a point that lands on the y = x diagonal in the
        magnitude plot but is colored "mismatched" is a pure phase error, and vice versa.
        Matrix elements with |g| < mag_tol (in either file) are dropped from the phase plot only,
        since the phase of a near-zero complex number is numerically meaningless (Re/Im and |g|
        remain well-defined, and are informative, even when |g| is tiny, so they keep every
        point). Because phase is an angle, two equal phases can be reported ~360 degrees apart
        across the +-180 degree branch cut; three parallel y = x, y = x - 360, y = x + 360
        reference lines are drawn so such wrapped points still show up "on the diagonal" instead
        of looking like large, spurious mismatches.

        The worst matrix element(s) (largest |g1 - g2|, e.g. the "Max absolute difference"
        printed below) are marked with a black-edged star in ALL FOUR plots and their location
        (global q/k indices, perturbation, and bra/ket band) is printed, since with thousands of
        tiny (s=1, alpha=0.3) points a lone outlier is otherwise easy to miss by eye even though
        it is colored red.

        Each panel is also annotated with a least-squares linear fit (y = slope*x + intercept)
        of the pooled (matched + mismatched) points, together with the Pearson correlation
        coefficient and RMSE -- a quick, quantitative check on top of the visual parity (a
        perfect reconstruction gives slope=1, intercept=0, r=1). For the phase panel the fit is
        computed after realigning points across the +-180 degree branch cut (see ``wrap`` in
        :func:`_linreg_stats`), so a physically tiny phase difference near the cut doesn't look
        like a huge disagreement in the reported slope/RMSE.

        Args:
            tol: Tolerance on the per-(k, q) max |g1 - g2| used to classify a symmetrized point
                as an exact match. Drives the color-coding in both the magnitude and phase plots.
            mag_tol: Matrix elements with |g| below this threshold (in either file) are excluded
                from the phase parity plot only.
            n_outliers: Number of worst matrix elements (by |g1 - g2|) to mark and print per spin.
            use_hexbin: If True, render the (typically large) matched-point cloud in each panel
                as a density-colored hexbin plot instead of a scatter plot -- much faster and
                more legible when there are many thousands of points. Mismatched points are
                always drawn as an explicit scatter overlay (in both modes) since they are the
                ones worth inspecting individually.
            bins: hexbin grid resolution (only used when ``use_hexbin=True``).
            verbose: If > 0, print the state (gstore_glob_state_kqs) of every (spin, q, k) index.
        """
        if len(self.abifiles) != 2:
            raise ValueError(f"compare_gvals_with_reconstruction requires exactly 2 GSTORE files, got {len(self.abifiles)}")

        g1, g2 = self.abifiles[0], self.abifiles[1]

        if g1.nsppol != g2.nsppol:
            raise ValueError("The two GSTORE files must have the same nsppol.")

        # Check dimensions are the same before plotting
        for attr in ["glob_nk_spin", "glob_spin_nq"]:
            if not np.array_equal(getattr(g1.r, attr), getattr(g2.r, attr)):
                cprint(f"Warning: {attr} differs between {g1.basename} and {g2.basename}", "red")

        # Each spin gets its own 2x2 block of rows [2*spin, 2*spin+1]: top row is |g| (col 0)
        # and phase (col 1); bottom row is Re(g) (col 0) and Im(g) (col 1). Scales differ a lot
        # across the four (linear magnitude/Re/Im vs degrees in [-180, 180]) so axes aren't shared.
        ax_mat, fig, plt = get_axarray_fig_plt(
            ax_list, nrows=2 * g1.nsppol, ncols=2, sharex=False, sharey=False, squeeze=False, rescale_fig=True
        )

        # Read internal table with the state of (k, q, spin). Shape: (nsppol, glob_nq, glob_nk).
        state1_kqs = g1.r.read_value("gstore_glob_state_kqs")
        state2_kqs = g2.r.read_value("gstore_glob_state_kqs")

        if verbose:
            print("state1_kqs vs state2_kqs")
            for idx in np.ndindex(state1_kqs.shape):
                print(f"{idx}: {state1_kqs[idx]}    {state2_kqs[idx]}")

        for spin in range(g1.nsppol):
            gvals1 = g1.gqk_spin[spin].gvals
            gvals2 = g2.gqk_spin[spin].gvals

            if gvals1.shape != gvals2.shape:
                raise ValueError(f"Shape mismatch for spin {spin}: {gvals1.shape} vs {gvals2.shape}")

            sym_mask = state1_kqs[spin] == GSTORE_KQ_SYMMETRIZED  # (glob_nq, glob_nk)
            n_sym = int(sym_mask.sum())
            if n_sym == 0:
                cprint(f"Spin {spin}: no GSTORE_KQ_SYMMETRIZED (k, q) points found, skipping.", "yellow")
                continue

            gvals1_sym = gvals1[sym_mask]  # (n_sym, natom3, nb_kq, nb_k)
            gvals2_sym = gvals2[sym_mask]

            abs_diff = np.abs(gvals1_sym - gvals2_sym)
            max_abs_diff = np.max(abs_diff)
            mean_abs_diff = np.mean(abs_diff)

            # Classify each symmetrized (k, q) point on its own worst matrix element.
            # This is what actually reveals whether the reconstruction works or not.
            per_point_maxdiff = np.max(abs_diff.reshape(n_sym, -1), axis=1)
            is_ok = per_point_maxdiff < tol
            n_ok = int(np.sum(is_ok))

            print(f"Spin {spin}: {n_sym} symmetrized (k, q) points "
                  f"(out of {sym_mask.size} total).")
            print(f"  Exact matches (per-point max|g1-g2| < {tol:.1e}): "
                  f"{n_ok}/{n_sym} ({100 * n_ok / n_sym:.1f}%)")
            print(f"  Max absolute difference |g1 - g2|: {max_abs_diff:.4e}")
            print(f"  Mean absolute difference |g1 - g2|: {mean_abs_diff:.4e}")

            # Broadcast the per-(k, q) pass/fail flag to every matrix element of that point so
            # the split found above is visible in both parity plots.
            point_size = gvals1_sym[0].size
            matched = np.repeat(is_ok, point_size)
            x = np.abs(gvals1_sym).flatten()
            y = np.abs(gvals2_sym).flatten()

            # Locate and report the n_outliers matrix elements with the largest |g1 - g2| (flat
            # index into abs_diff/x/y/phase1/phase2, which all share abs_diff's element order).
            n_show = min(n_outliers, abs_diff.size)
            flat_diff = abs_diff.ravel()
            outlier_idx = np.argpartition(flat_diff, -n_show)[-n_show:]
            outlier_idx = outlier_idx[np.argsort(-flat_diff[outlier_idx])]
            sym_qk = np.argwhere(sym_mask)  # (n_sym, 2): row i --> (iq_g, ik_g) of gvals*_sym[i]

            for rank, flat_idx in enumerate(outlier_idx, start=1):
                i_sym, mu_o, m_o, n_o = np.unravel_index(flat_idx, abs_diff.shape)
                iq_g_o, ik_g_o = sym_qk[i_sym]
                print(f"  Outlier #{rank}: |g1-g2|={flat_diff[flat_idx]:.4e} at "
                      f"(iq_g={iq_g_o}, ik_g={ik_g_o}, mu={mu_o}, bra_band={m_o}, ket_band={n_o}) "
                      f"[|g1|={x[flat_idx]:.4e}, |g2|={y[flat_idx]:.4e}]")

            def mark_outliers(ax, valx, valy):
                """Star + rank-label the n_show worst matrix elements on ax (shared closure state)."""
                ax.scatter(valx[outlier_idx], valy[outlier_idx], s=100, marker="*", facecolors="yellow",
                           edgecolors="black", linewidths=0.8, zorder=5, label=f"worst {n_show} pt(s)")
                for rank, flat_idx in enumerate(outlier_idx, start=1):
                    ax.annotate(str(rank), (valx[flat_idx], valy[flat_idx]), textcoords="offset points",
                                xytext=(4, 4), fontsize=fontsize, fontweight="bold", zorder=6)

            # --- Magnitude parity plot --------------------------------------------------
            ax = ax_mat[2 * spin, 0]

            max_val = max(np.max(x), np.max(y)) if len(x) > 0 else 1.0

            if use_hexbin:
                # Density-colored hexbin for the (typically large) full point cloud -- much
                # faster to render/save than a per-point scatter. Mismatched points are still
                # drawn explicitly on top so they remain individually visible.
                h = ax.hexbin(x, y, gridsize=bins, cmap="Blues", mincnt=1, bins="log",
                              extent=(0, max_val, 0, max_val))
                fig.colorbar(h, ax=ax, shrink=0.85).set_label("point count (log)", fontsize=fontsize)
                ax.scatter(x[~matched], y[~matched], s=4, alpha=0.6, color="C3",
                           label=f"mismatched (k,q): {n_sym - n_ok}/{n_sym}", zorder=4)
            else:
                # Use rasterized=True to keep vector graphics small if many points
                ax.scatter(x[matched], y[matched], s=1, alpha=0.3, rasterized=True, color="C0",
                           label=f"matched (k,q): {n_ok}/{n_sym}")
                ax.scatter(x[~matched], y[~matched], s=1, alpha=0.3, rasterized=True, color="C3",
                           label=f"mismatched (k,q): {n_sym - n_ok}/{n_sym}")

            ax.plot([0, max_val], [0, max_val], 'k--', lw=0.5, label="y = x")
            mark_outliers(ax, x, y)
            _add_fit_annotation(ax, x, y, xlim=(0, max_val), fontsize=fontsize)

            ax.set_xlabel(f"$|g|$ from {g1.basename}", fontsize=fontsize)
            ax.set_ylabel(f"$|g|$ from {g2.basename}", fontsize=fontsize)
            ax.set_title(f"Spin {spin}: |g| parity, symmetrized (k,q) only", fontsize=fontsize)
            ax.legend(fontsize=fontsize, loc="lower right")
            ax.set_xlim(0, max_val)
            ax.set_ylim(0, max_val)
            ax.set_aspect("equal", adjustable="box")

            # --- Phase parity plot -------------------------------------------------------
            ax = ax_mat[2 * spin, 1]

            # Phase is meaningless noise when |g| ~ 0 in either file: drop those points here only.
            has_signal = (x > mag_tol) | (y > mag_tol)
            phase_matched = matched & has_signal
            phase_mismatched = (~matched) & has_signal
            n_signal = int(np.sum(has_signal))

            phase1 = np.degrees(np.angle(gvals1_sym)).flatten()
            phase2 = np.degrees(np.angle(gvals2_sym)).flatten()

            if use_hexbin:
                if n_signal > 0:
                    # hexbin/colorbar choke on a fully-empty array (can happen if mag_tol
                    # filters out every point), so only call it when there's something to bin.
                    h = ax.hexbin(phase1[has_signal], phase2[has_signal], gridsize=bins, cmap="Blues",
                                  mincnt=1, bins="log", extent=(-180, 180, -180, 180))
                    fig.colorbar(h, ax=ax, shrink=0.85).set_label("point count (log)", fontsize=fontsize)
                else:
                    ax.text(0.5, 0.5, f"no points with |g| > {mag_tol:.0e}", transform=ax.transAxes,
                            ha="center", va="center", fontsize=fontsize, color="0.5")
                ax.scatter(phase1[phase_mismatched], phase2[phase_mismatched], s=4, alpha=0.6,
                           color="C3", label=f"mismatched (k,q): {n_sym - n_ok}/{n_sym}", zorder=4)
            else:
                ax.scatter(phase1[phase_matched], phase2[phase_matched], s=1, alpha=0.3, rasterized=True,
                           color="C0", label=f"matched (k,q): {n_ok}/{n_sym}")
                ax.scatter(phase1[phase_mismatched], phase2[phase_mismatched], s=1, alpha=0.3, rasterized=True,
                           color="C3", label=f"mismatched (k,q): {n_sym - n_ok}/{n_sym}")

            # +-360 deg helper lines so points wrapped across the +-180 deg branch cut still
            # land "on the diagonal" instead of looking like spurious large mismatches.
            for offset, ref_label in ((0, "y = x"), (-360, None), (360, None)):
                ax.plot([-180, 180], [-180 + offset, 180 + offset], 'k--', lw=0.5, label=ref_label)
            mark_outliers(ax, phase1, phase2)
            _add_fit_annotation(ax, phase1[has_signal], phase2[has_signal], xlim=(-180, 180),
                                 wrap=True, fontsize=fontsize)

            ax.set_xlim(-180, 180)
            ax.set_ylim(-180, 180)
            ax.set_xlabel(f"phase(g) [deg] from {g1.basename}", fontsize=fontsize)
            ax.set_ylabel(f"phase(g) [deg] from {g2.basename}", fontsize=fontsize)
            ax.set_title(f"Spin {spin}: phase parity (|g| > {mag_tol:.0e}, {n_signal} pts)", fontsize=fontsize)
            ax.legend(fontsize=fontsize, loc="lower right")
            ax.set_aspect("equal", adjustable="box")

            # --- Re(g) parity plot -------------------------------------------------------
            ax = ax_mat[2 * spin + 1, 0]
            real1 = gvals1_sym.real.flatten()
            real2 = gvals2_sym.real.flatten()

            lim = max(np.max(np.abs(real1)), np.max(np.abs(real2))) if len(real1) > 0 else 1.0

            if use_hexbin:
                h = ax.hexbin(real1, real2, gridsize=bins, cmap="Blues", mincnt=1, bins="log",
                              extent=(-lim, lim, -lim, lim))
                fig.colorbar(h, ax=ax, shrink=0.85).set_label("point count (log)", fontsize=fontsize)
                ax.scatter(real1[~matched], real2[~matched], s=4, alpha=0.6, color="C3",
                           label=f"mismatched (k,q): {n_sym - n_ok}/{n_sym}", zorder=4)
            else:
                ax.scatter(real1[matched], real2[matched], s=1, alpha=0.3, rasterized=True, color="C0",
                           label=f"matched (k,q): {n_ok}/{n_sym}")
                ax.scatter(real1[~matched], real2[~matched], s=1, alpha=0.3, rasterized=True, color="C3",
                           label=f"mismatched (k,q): {n_sym - n_ok}/{n_sym}")

            ax.plot([-lim, lim], [-lim, lim], 'k--', lw=0.5, label="y = x")
            mark_outliers(ax, real1, real2)
            _add_fit_annotation(ax, real1, real2, xlim=(-lim, lim), fontsize=fontsize)

            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_xlabel(f"Re(g) from {g1.basename}", fontsize=fontsize)
            ax.set_ylabel(f"Re(g) from {g2.basename}", fontsize=fontsize)
            ax.set_title(f"Spin {spin}: Re(g) parity", fontsize=fontsize)
            ax.legend(fontsize=fontsize, loc="lower right")
            ax.set_aspect("equal", adjustable="box")

            # --- Im(g) parity plot -------------------------------------------------------
            ax = ax_mat[2 * spin + 1, 1]
            imag1 = gvals1_sym.imag.flatten()
            imag2 = gvals2_sym.imag.flatten()

            lim = max(np.max(np.abs(imag1)), np.max(np.abs(imag2))) if len(imag1) > 0 else 1.0

            if use_hexbin:
                h = ax.hexbin(imag1, imag2, gridsize=bins, cmap="Blues", mincnt=1, bins="log",
                              extent=(-lim, lim, -lim, lim))
                fig.colorbar(h, ax=ax, shrink=0.85).set_label("point count (log)", fontsize=fontsize)
                ax.scatter(imag1[~matched], imag2[~matched], s=4, alpha=0.6, color="C3",
                           label=f"mismatched (k,q): {n_sym - n_ok}/{n_sym}", zorder=4)
            else:
                ax.scatter(imag1[matched], imag2[matched], s=1, alpha=0.3, rasterized=True, color="C0",
                           label=f"matched (k,q): {n_ok}/{n_sym}")
                ax.scatter(imag1[~matched], imag2[~matched], s=1, alpha=0.3, rasterized=True, color="C3",
                           label=f"mismatched (k,q): {n_sym - n_ok}/{n_sym}")

            ax.plot([-lim, lim], [-lim, lim], 'k--', lw=0.5, label="y = x")
            mark_outliers(ax, imag1, imag2)
            _add_fit_annotation(ax, imag1, imag2, xlim=(-lim, lim), fontsize=fontsize)

            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_xlabel(f"Im(g) from {g1.basename}", fontsize=fontsize)
            ax.set_ylabel(f"Im(g) from {g2.basename}", fontsize=fontsize)
            ax.set_title(f"Spin {spin}: Im(g) parity", fontsize=fontsize)
            ax.legend(fontsize=fontsize, loc="lower right")
            ax.set_aspect("equal", adjustable="box")

        return fig

    def compare_gvals_gauge_invariant(self, tol: float = 1e-6, deg_tol: float = 1e-4, verbose: int = 0) -> None:
        """
        Gauge-invariant counterpart of ``compare_gvals_with_reconstruction``.

        Individual g_{mn} matrix elements are only defined up to an arbitrary unitary
        rotation within each degenerate subspace at k (the n/ket index) and at k+q
        (the m/bra index): the analytic reconstruction formula used by gstore_symmetrize
        assumes this rotation is the identity, which is exact for isolated non-degenerate
        bands but not for degenerate ones. Comparing raw matrix elements therefore flags
        many "mismatches" that are not bugs, just a different (equally valid) choice of
        basis within a degenerate block.

        This method groups the bra/ket bands into degenerate multiplets (using the IBZ
        eigenvalues) and, for each (bra-multiplet, ket-multiplet, perturbation) block,
        compares the singular values of the corresponding sub-matrix of g(k,q) instead of
        its raw elements. Singular values are invariant under G -> U G V^dagger for any
        unitary U, V, i.e. under independent gauge rotations of the bra and ket subspaces,
        while still being fully sensitive to genuine errors (including phase errors for
        non-degenerate, i.e. 1x1, blocks).

        Args:
            tol: Tolerance on the per-(k, q) max singular-value difference used to classify
                a symmetrized point as an exact match.
            deg_tol: Energy tolerance (eV) used to group bands into degenerate multiplets.
            verbose: If > 0, print per-point diagnostics.
        """
        if len(self.abifiles) != 2:
            raise ValueError(f"compare_gvals_gauge_invariant requires exactly 2 GSTORE files, got {len(self.abifiles)}")

        g1, g2 = self.abifiles[0], self.abifiles[1]

        if g1.nsppol != g2.nsppol:
            raise ValueError("The two GSTORE files must have the same nsppol.")

        def find_multiplets(energies: np.ndarray, etol: float) -> list:
            """Group band indices (local, 0-based) into degenerate multiplets."""
            order = np.argsort(energies)
            groups = [[int(order[0])]]
            for i in range(1, len(order)):
                if abs(energies[order[i]] - energies[order[i - 1]]) < etol:
                    groups[-1].append(int(order[i]))
                else:
                    groups.append([int(order[i])])
            return groups

        for spin in range(g1.nsppol):
            gqk1, gqk2 = g1.gqk_spin[spin], g2.gqk_spin[spin]
            gvals1, gvals2 = gqk1.gvals, gqk2.gvals

            if gvals1.shape != gvals2.shape:
                raise ValueError(f"Shape mismatch for spin {spin}: {gvals1.shape} vs {gvals2.shape}")

            state1_kqs = g1.r.read_value("gstore_glob_state_kqs")
            sym_mask = state1_kqs[spin] == GSTORE_KQ_SYMMETRIZED  # (glob_nq, glob_nk)
            if not np.any(sym_mask):
                cprint(f"Spin {spin}: no GSTORE_KQ_SYMMETRIZED (k, q) points found, skipping.", "yellow")
                continue

            r = g1.r
            eigens_ibz = r.read_value("eigenvalues") * 27.211386245988  # Ha --> eV, indexed like r.kibz.
            bstart_k, nb_k = gqk1.bstart_k, gqk1.nb_k
            bstart_kq, nb_kq = gqk1.bstart_kq, gqk1.nb_kq

            kbz, kbz2ibz, kglob2bz = r.kbz, r.kbz2ibz, r.kglob2bz[spin]
            qbz, qglob2bz = r.qbz, r.qglob2bz[spin]

            n_multiplet_cache, m_multiplet_cache = {}, {}

            def get_multiplets(cache: dict, ik_ibz: int, bstart: int, nb: int) -> list:
                if ik_ibz not in cache:
                    e = eigens_ibz[spin, ik_ibz, bstart:bstart + nb]
                    cache[ik_ibz] = find_multiplets(e, deg_tol)
                return cache[ik_ibz]

            def find_kbz_index(kpt: np.ndarray) -> int:
                diffs = kbz - kpt[None, :]
                diffs -= np.round(diffs)
                idx = np.where(np.all(np.abs(diffs) < 1e-6, axis=1))[0]
                if len(idx) == 0:
                    raise ValueError(f"Cannot locate k-point {kpt} in the kbz mesh")
                return int(idx[0])

            natom3 = gvals1.shape[2]
            per_point_maxdiff = []

            for ik_g in range(gqk1.glob_nk):
                ik_bz = kglob2bz[ik_g]
                ik_ibz = kbz2ibz[ik_bz, 0]
                kpt = kbz[ik_bz]
                n_mult = get_multiplets(n_multiplet_cache, ik_ibz, bstart_k, nb_k)

                for iq_g in range(gqk1.glob_nq):
                    if not sym_mask[iq_g, ik_g]:
                        continue

                    qpt = qbz[qglob2bz[iq_g]]
                    ikq_bz = find_kbz_index(kpt + qpt)
                    ikq_ibz = kbz2ibz[ikq_bz, 0]
                    m_mult = get_multiplets(m_multiplet_cache, ikq_ibz, bstart_kq, nb_kq)

                    worst = 0.0
                    for mu in range(natom3):
                        G1, G2 = gvals1[iq_g, ik_g, mu], gvals2[iq_g, ik_g, mu]
                        for mb in m_mult:
                            for nb_ in n_mult:
                                block1 = G1[np.ix_(mb, nb_)]
                                block2 = G2[np.ix_(mb, nb_)]
                                sv1 = np.linalg.svd(block1, compute_uv=False)
                                sv2 = np.linalg.svd(block2, compute_uv=False)
                                worst = max(worst, float(np.max(np.abs(sv1 - sv2))))

                    per_point_maxdiff.append(worst)
                    if verbose and worst >= tol:
                        print(f"  spin={spin} ik_g={ik_g} iq_g={iq_g}: max singular-value diff = {worst:.4e}")

            per_point_maxdiff = np.array(per_point_maxdiff)
            n_tot = len(per_point_maxdiff)
            n_ok = int(np.sum(per_point_maxdiff < tol))

            print(f"Spin {spin}: {n_tot} symmetrized (k, q) points (out of {sym_mask.size} total).")
            print(f"  Gauge-invariant exact matches (per-point max singular-value diff < {tol:.1e}): "
                  f"{n_ok}/{n_tot} ({100 * n_ok / n_tot:.1f}%)")
            print(f"  Max singular-value difference: {per_point_maxdiff.max():.4e}")
            print(f"  Mean singular-value difference: {per_point_maxdiff.mean():.4e}")
