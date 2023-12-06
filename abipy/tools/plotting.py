# coding: utf-8
"""
Utilities for generating matplotlib plots.

.. note::

    Avoid importing matplotlib or plotly in the module namespace otherwise startup is very slow.
"""
from __future__ import annotations

import os
import time
import itertools
import functools
import numpy as np
import pandas as pd

from collections import namedtuple, OrderedDict
from typing import Any, Callable, Iterator
from monty.string import list_strings
from pymatgen.util.plotting import add_fig_kwargs
from abipy.tools import duck
from abipy.tools.iotools import dataframe_from_filepath
from abipy.tools.typing import Figure, Axes, VectorLike
from abipy.tools.numtools import data_from_cplx_mode

import matplotlib.collections as mcoll
from plotly.tools import mpl_to_plotly

__all__ = [
    "set_axlims",
    "add_fig_kwargs",
    "get_ax_fig_plt",
    "get_axarray_fig_plt",
    "get_ax3d_fig_plt",
    "plot_array",
    "ArrayPlotter",
    "data_from_cplx_mode",
    "Marker",
    "plot_unit_cell",
    "GenericDataFilePlotter",
    "GenericDataFilesPlotter",
    "add_plotly_fig_kwargs",
    "get_fig_plotly",
    "get_figs_plotly",
]


# https://matplotlib.org/gallery/lines_bars_and_markers/linestyles.html
linestyles = OrderedDict(
    [('solid',               (0, ())),
     ('loosely_dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely_dotted',      (0, (1, 1))),
     #
     ('loosely_dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely_dashed',      (0, (5, 1))),
     #
     ('loosely_dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely_dashdotted',  (0, (3, 1, 1, 1))),
     #
     ('loosely_dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely_dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]
)


class FilesPlotter:
    """
    Use matplotlib to plot multiple png files on a grid.

    Example:

        FilesPlotter(["file1.png", file2.png"]).plot()
    """
    def __init__(self, filepaths: list[str]):
        self.filepaths = list_strings(filepaths)

    @add_fig_kwargs
    def plot(self, **kwargs) -> Figure:
        """Loop through the PNG files and display them in subplots."""
        # Build grid of plots.
        num_plots, ncols, nrows = len(self.filepaths), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        for i, (filepath, ax) in enumerate(zip(self.filepaths, ax_list)):
            ax.axis('off')
            ax.imshow(plt.imread(filepath))

        return fig


@functools.cache
def get_color_symbol(style: str="VESTA") -> dict:
    """
    Dictionary mapping chemical symbols to RGB color.

    Args:
        style: "VESTA" or "Jmol".
    """
    from monty.serialization import loadfn
    from pymatgen import vis
    colors = loadfn(os.path.join(os.path.dirname(vis.__file__), "ElementColorSchemes.yaml"))
    if style not in colors:
        raise KeyError(f"Invalid {style=}. Should be in {colors.keys()}")
    color_symbol = {el: [j / 256.001 for j in colors[style][el]] for el in colors[style]}
    return color_symbol


###################
# Matplotlib tools
###################

def get_ax_fig_plt(ax=None, **kwargs):
    """
    Helper function used in plot functions supporting an optional Axes argument.
    If ax is None, we build the `matplotlib` figure and create the Axes else
    we return the current active figure.

    Args:
        ax (Axes, optional): Axes object. Defaults to None.
        kwargs: keyword arguments are passed to plt.figure if ax is not None.

      Returns:
        ax: :class:`Axes` object
        figure: matplotlib figure
        plt: matplotlib pyplot module.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure(**kwargs)
        ax = fig.gca()
    else:
        fig = plt.gcf()

    return ax, fig, plt


def get_ax3d_fig_plt(ax=None, **kwargs):
    """
    Helper function used in plot functions supporting an optional Axes3D
    argument. If ax is None, we build the `matplotlib` figure and create the
    Axes3D else we return the current active figure.

    Args:
        ax (Axes3D, optional): Axes3D object. Defaults to None.
        kwargs: keyword arguments are passed to plt.figure if ax is not None.

    Returns:
        tuple[Axes3D, Figure]: matplotlib Axes3D and corresponding figure objects
    """
    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure(**kwargs)
        ax = fig.add_subplot(projection="3d")
    else:
        fig = plt.gcf()

    return ax, fig, plt


def get_axarray_fig_plt(
    ax_array, nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None, **fig_kw
):
    """
    Helper function used in plot functions that accept an optional array of Axes
    as argument. If ax_array is None, we build the `matplotlib` figure and
    create the array of Axes by calling plt.subplots else we return the
    current active figure.

    Returns:
        ax: Array of Axes objects
        figure: matplotlib figure
        plt: matplotlib pyplot module.
    """
    import matplotlib.pyplot as plt

    if ax_array is None:
        fig, ax_array = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex=sharex,
            sharey=sharey,
            squeeze=squeeze,
            subplot_kw=subplot_kw,
            gridspec_kw=gridspec_kw,
            **fig_kw,
        )
    else:
        fig = plt.gcf()
        ax_array = np.reshape(np.array(ax_array), (nrows, ncols))
        if squeeze:
            if ax_array.size == 1:
                ax_array = ax_array[0]
            elif any(s == 1 for s in ax_array.shape):
                ax_array = ax_array.ravel()

    return ax_array, fig, plt


def is_mpl_figure(obj: Any) -> bool:
    """Return True if obj is a matplotlib Figure."""
    from matplotlib import pyplot as plt
    return isinstance(obj, plt.Figure)


def ax_append_title(ax, title, loc="center", fontsize=None) -> str:
    """Add title to previous ax.title. Return new title."""
    prev_title = ax.get_title(loc=loc)
    new_title = prev_title + title
    ax.set_title(new_title, loc=loc, fontsize=fontsize)
    return new_title


def ax_share(xy_string: str, *ax_list) -> None:
    """
    Share x- or y-axis of two or more subplots after they are created.

    Args:
        xy_string: "x" to share x-axis, "xy" for both
        ax_list: List of axes to share.

    Example:

        ax_share("y", ax0, ax1)
        ax_share("xy", *(ax0, ax1, ax2))
    """
    if "x" in xy_string:
        for ix, ax in enumerate(ax_list):
            others = [a for a in ax_list if a != ax]
            ax.get_shared_x_axes().join(*others)

    if "y" in xy_string:
        for ix, ax in enumerate(ax_list):
            others = [a for a in ax_list if a != ax]
            ax.get_shared_y_axes().join(*others)


def set_axlims(ax, lims: tuple, axname: str) -> tuple:
    """
    Set the data limits for the axis ax.

    Args:
        lims: tuple(2) for (left, right), tuple(1) or scalar for left only.
        axname: "x" for x-axis, "y" for y-axis.

    Return: (left, right)
    """
    left, right = None, None
    if lims is None: return left, right

    len_lims = None
    try:
        len_lims = len(lims)
    except TypeError:
        # Assume Scalar
        left = float(lims)

    if len_lims is not None:
        if len(lims) == 2:
            left, right = lims[0], lims[1]
        elif len(lims) == 1:
            left = lims[0]

    set_lim = getattr(ax, {"x": "set_xlim", "y": "set_ylim"}[axname])
    if left != right:
        set_lim(left, right)

    return left, right


def set_ax_xylabels(ax, xlabel: str, ylabel: str, exchange_xy: bool = False) -> None:
    """
    Set the x- and the y-label of axis ax, exchanging x and y if exchange_xy.
    """
    if exchange_xy: xlabel, ylabel = ylabel, xlabel
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def set_logscale(ax_or_axlist, xy_log) -> None:
    """
    Activate logscale

    Args:
        ax_or_axlist: Axes or list of axes.
        xy_log: None or empty string for linear scale. "x" for log scale on x-axis.
            "xy" for log scale on x- and y-axis. "x:semilog" for semilog scale on x-axis.
    """
    if not xy_log: return

    # Parse xy_log string.
    xy, log_type = xy_log, "log"
    if ":" in xy_log:
        xy, log_type = xy_log.split(":")

    ax_list = [ax_or_axlist] if not duck.is_listlike(ax_or_axlist) else ax_or_axlist

    for ix, ax in enumerate(ax_list):
        if "x" in xy:
            ax.set_xscale(log_type)
        if "y" in xy:
            ax.set_yscale(log_type)


def set_ticks_fontsize(ax_or_axlist, fontsize: int, xy_string="xy", **kwargs) -> None:
    """
    Set tick properties for one axis or a list of axis.

    Args:
        ax_or_axlist: Axes or list of axes.
        xy_string: "x" to share x-axis, "xy" for both.
    """
    ax_list = [ax_or_axlist] if not duck.is_listlike(ax_or_axlist) else ax_or_axlist

    for ix, ax in enumerate(ax_list):
        if "x" in xy_string:
            ax.tick_params(axis='x', labelsize=fontsize, **kwargs)

        if "y" in xy_string:
            ax.tick_params(axis='y', labelsize=fontsize, **kwargs)


def set_grid_legend(ax_or_axlist, fontsize: int,
                    xlabel=None, ylabel=None, grid=True, legend=True, direction=None, title=None, legend_loc="best") -> None:
    """
    Activate grid and legend for one axis or a list of axis.

    Args:
        grid: True to activate the grid.
        legend: True to activate the legend.
        direction: Use "x" ("y") if to add xlabel (ylabel) only to the last ax.
        title: Title string
    """
    if duck.is_listlike(ax_or_axlist):
        for ix, ax in enumerate(ax_or_axlist):
            ax.grid(grid)
            if legend: ax.legend(loc=legend_loc, fontsize=fontsize, shadow=True)
            if xlabel:
                doit = direction is None or (direction == "y" and ix == len(ax_or_axlist) -1)
                if doit: ax.set_xlabel(xlabel)
            if ylabel:
                doit = direction is None or (direction == "x" and ix == len(ax_or_axlist) -1)
                if doit: ax.set_ylabel(ylabel)
            if title: ax.set_title(title, fontsize=fontsize)
    else:
        ax = ax_or_axlist
        ax.grid(grid)
        if legend: ax.legend(loc=legend_loc, fontsize=fontsize, shadow=True)
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
        if title: ax.set_title(title, fontsize=fontsize)


def set_visible(ax, boolean: bool, *args) -> None:
    """
    Hide/Show the artists of axis ax listed in args.
    ax can be a single axis, a list or axis or a numpy arrays.
    """
    if duck.is_listlike(ax):
        if isinstance(ax, np.ndarray):
            for _ in ax.ravel():
                set_visible(_, *args)
        else:
            for _ in ax:
                set_visible(_, *args)
        return

    if "legend" in args and ax.legend():
        ax.legend().set_visible(boolean)
    if "title" in args and ax.title:
        ax.title.set_visible(boolean)
    if "xlabel" in args and ax.xaxis.label:
        ax.xaxis.label.set_visible(boolean)
    if "ylabel" in args and ax.yaxis.label:
        ax.yaxis.label.set_visible(boolean)
    if "xticklabels" in args:
        for label in ax.get_xticklabels():
            label.set_visible(boolean)
    if "yticklabels" in args:
        for label in ax.get_yticklabels():
            label.set_visible(boolean)


def rotate_ticklabels(ax, rotation: float, axname: str ="x") -> None:
    """Rotate the ticklables of axis ``ax``"""
    if "x" in axname:
        for tick in ax.get_xticklabels():
            tick.set_rotation(rotation)
    if "y" in axname:
        for tick in ax.get_yticklabels():
            tick.set_rotation(rotation)


def hspan_ax_line(ax, line, abs_conv, hatch, alpha=0.2, with_label=True) -> None:
    """
    Add hspan to ax showing the convergence region of width `abs_conv`.
    Use same color as line. Return immediately if abs_conv is None or x-values are strings.
    """
    if abs_conv is None: return
    xs = line.get_xdata()
    ys = line.get_ydata()
    if duck.is_string(xs[0]): return

    color = line.get_color()
    span_style = dict(alpha=0.2, color=color, hatch=hatch)

    x_max = xs[-1]
    x_inds = np.where(xs == x_max)[0]
    # This to support the case in which we have multiple ys for the same x_max
    for i, ix in enumerate(x_inds):
        y_xmax = ys[ix]
        ax.axhspan(y_xmax - abs_conv, y_xmax + abs_conv,
                   label=r"$|y-y(x_{max})| \leq %s$" % abs_conv if (with_label and i == 0) else None,
                   **span_style)


@add_fig_kwargs
def plot_xy_with_hue(data: pd.DataFrame, x: str, y: str, hue: str,
                     decimals=None, ax=None, xlims=None, ylims=None, fontsize=8, **kwargs) -> Figure:
    """
    Plot y = f(x) relation for different values of `hue`.
    Useful for convergence tests done wrt two parameters.

    Args:
        data: |pandas-DataFrame| containing columns `x`, `y`, and `hue`.
        x: Name of the column used as x-value
        y: Name of the column(s) used as y-value
        hue: Variable that define subsets of the data, which will be drawn on separate lines
        decimals: Number of decimal places to round `hue` columns. Ignore if None
        ax: |matplotlib-Axes| or None if a new figure should be created.
        xlims, ylims: Set the data limits for the x(y)-axis. Accept tuple e.g. `(left, right)`
            or scalar e.g. `left`. If left (right) is None, default values are used
        fontsize: Legend fontsize.
        kwargs: Keyword arguments are passed to ax.plot method.

    Returns: |matplotlib-Figure|
    """
    if isinstance(y, (list, tuple)):
        # Recursive call for each ax in ax_list.
        num_plots, ncols, nrows = len(y), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)

        ax_list = ax_list.ravel()
        if num_plots % ncols != 0: ax_list[-1].axis('off')

        for ykey, ax in zip(y, ax_list):
            plot_xy_with_hue(data, x, str(ykey), hue, decimals=decimals, ax=ax,
                             xlims=xlims, ylims=ylims, fontsize=fontsize, show=False, **kwargs)
        return fig

    # Check here because pandas error messages are a bit criptic.
    miss = [k for k in (x, y, hue) if k not in data]
    if miss:
        raise ValueError("Cannot find `%s` in dataframe.\nAvailable keys are: %s" % (str(miss), str(data.keys())))

    # Truncate values in hue column so that we can group.
    if decimals is not None:
        data = data.round({hue: decimals})

    ax, fig, plt = get_ax_fig_plt(ax=ax)
    for key, grp in data.groupby(by=hue):
        # Sort xs and rearrange ys
        xy = np.array(sorted(zip(grp[x], grp[y]), key=lambda t: t[0]))
        xvals, yvals = xy[:, 0], xy[:, 1]

        label = "%s" % (str(key))
        if not kwargs:
            ax.plot(xvals, yvals, 'o-', label=label)
        else:
            ax.plot(xvals, yvals, label=label, **kwargs)

    ax.grid(True)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    set_axlims(ax, xlims, "x")
    set_axlims(ax, ylims, "y")
    ax.legend(loc="best", fontsize=fontsize, shadow=True)

    return fig


def linear_fit_ax(ax, xs, ys, fontsize, with_label=True, with_ideal_line=False, **kwargs) -> tuple[float]:
    """
    Calculate a linear least-squares regression for two sets of measurements.
    kwargs are passed to ax.plot.
    """
    from scipy.stats import linregress
    fit = linregress(xs, ys)
    label = r"Linear fit $\alpha={:.2f}$, $r^2$={:.2f}".format(fit.slope, fit.rvalue**2)
    if "color" not in kwargs:
        kwargs["color"] = "r"

    ax.plot(xs, fit.slope*xs + fit.intercept, label=label if with_label else None, **kwargs)
    if with_ideal_line:
        # Plot y = x line
        ax.plot([xs[0], xs[-1]], [ys[0], ys[-1]], color='k', linestyle='-',
                linewidth=1, label='Ideal' if with_label else None)
    return fit


@add_fig_kwargs
def plot_array(array, color_map=None, cplx_mode="abs", **kwargs) -> Figure:
    """
    Use imshow for plotting 2D or 1D arrays. Return: |matplotlib-Figure|

    Example::

        plot_array(np.random.rand(10,10))

    See <http://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python>

    Args:
        array: Array-like object (1D or 2D).
        color_map: color map.
        cplx_mode:
            Flag defining how to handle complex arrays. Possible values in ("re", "im", "abs", "angle")
            "re" for the real part, "im" for the imaginary part.
            "abs" means that the absolute value of the complex number is shown.
            "angle" will display the phase of the complex number in radians.
    """
    # Handle vectors
    array = np.atleast_2d(array)
    array = data_from_cplx_mode(cplx_mode, array)

    import matplotlib as mpl
    from matplotlib import pyplot as plt
    if color_map is None:
        # make a color map of fixed colors
        color_map = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                                 ['blue', 'black', 'red'], 256)

    img = plt.imshow(array, interpolation='nearest', cmap=color_map, origin='lower')

    # Make a color bar
    plt.colorbar(img, cmap=color_map)

    # Set grid
    plt.grid(True, color='white')

    fig = plt.gcf()
    return fig


class ConvergenceAnalyzer:
    """
    This object allows one to plot the convergence of an arbitrary list
    of quantities as a function of the same x.
    """

    # Colors for the different convergence criteria.
    color_ilevel = ["red", "blue", "green"]

    # matplotlib option to fill convergence window.
    HATCH = "/"

    @classmethod
    def from_xy_label_vals(cls, xlabel, xs, ylabel, yvalues, tols) -> ConvergenceAnalyzer:
        """
        Simplified interface to analyze a single list of values.
        """
        yvals_dict = {ylabel: yvalues}
        ytols_dict = {ylabel: tols}
        return cls(xlabel, xs, yvals_dict, ytols_dict)

    @classmethod
    def from_file(cls, filepath: str, xkey: str, ytols_dict: dict, **kwargs) -> ConvergenceAnalyzer:
        """
        High-level constructor to build the object from a file containing data
        that can be converted to pandas DataFrame. kwargs are passed to the pandas IO routines.

        Args:
            filepath: Filename.
            xkey: name of the x-variable.
            ytols_dict: dict mapping the name of the y-variable to absolute tolerance(s).
        """
        df = dataframe_from_filepath(filepath, **kwargs)
        return cls.from_dataframe(df, xkey, ytols_dict)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, xkey: str, ytols_dict: dict) -> ConvergenceAnalyzer:
        """
        Build the object from a pandas dataframe.

        Args:
            df: DataFrame
            xkey: name of the x-variable.
            ytols_dict: dict mapping the name of the y-variable to tolerance(s).
        """
        df = df.sort_values(xkey)
        xs = df[xkey].values
        yvals_dict = {k: df[k].values for k in ytols_dict}
        return cls(xkey, xs, yvals_dict, ytols_dict)

    def __init__(self, xkey: str, xs: VectorLike, yvals_dict: dict[str, VectorLike], ytols_dict: dict):
        """
        Args:
            xkey:
            xs:
            yvals_dict:
            ytols_dict: dict mapping the name of the y-variable to absolute tolerance(s).

        Example::

            plotter = ConvergencePlotter("ecut", ecut_value, yvals_dict, ytols_dict)
            plotter.plot()
        """
        # Convert to numpy arrays and store data in self.
        self.xkey = self.xlabel = xkey
        self.xs = np.array(xs)
        if not np.all(self.xs[:-1] <= self.xs[1:]):
            raise ValueError("xs values should be in ascending order")

        self.yvals_dict = {k: np.array(v) for k, v in yvals_dict.items()}
        self.ykey2label = {k: k for k in yvals_dict}

        if len(self.yvals_dict) > len(self.color_ilevel):
            raise ValueError(f"Not programmed for more than {len(self.color_ilevel)} convergence levels")

        # Handle ytols_dict.
        self.ytols_dict = {}
        for ykey, ytols in ytols_dict.items():
            if not duck.is_listlike(ytols): ytols = [ytols]

            if any(yt <= 0 for yt in ytols):
                raise ValueError(f"tolerances cannot be negative: {ytols}")

            # Sort input tolerances just to be on the safe side.
            self.ytols_dict[ykey] = np.sort(np.array(ytols))[::-1]

        # Compute the first index in xs that gives value within the convergence window.
        # -1 or None indicates that convergence has not been achieved.
        self.ykey_ixs = {}
        self.ykey_best_xs = {}

        for ykey, ys in self.yvals_dict.items():
            if len(ys) != len(xs):
                raise ValueError(f"len(ys) != len(xs): {len(ys)} and {len(xs)}")

            tol_levels = self.ytols_dict[ykey]

            # Init values assuming no convergence achieved.
            self.ykey_ixs[ykey] = [-1] * len(tol_levels)
            self.ykey_best_xs[ykey] = [None] * len(tol_levels)

            # For each y-tolerance.
            num_x = len(self.xs)
            y_xmax = ys[-1]
            for il, ytol in enumerate(tol_levels):
                for _, xx in enumerate(self.xs[::-1]):
                    ix = -_ + num_x - 1
                    if abs(y_xmax - ys[ix]) > ytol:
                        self.ykey_ixs[ykey][il] = ix + 1
                        break

                ix = self.ykey_ixs[ykey][il]
                if ix != -1:
                    # If converged, use linear interpolation to get a better estimate of the
                    # converged xx. This is useful especially if the xs grid is too coarse.
                    best_xx = self.xs[ix]
                    if ix - 1 >= 0:
                        x0, y0 = xs[ix-1], ys[ix-1]
                        x1, y1 = xs[ix], ys[ix]
                        alpha = (y1 - y0) / (x1 - x0)
                        # y(x) = alpha * (x - x0) + y0
                        #print("best_xx 1", best_xx)
                        if (y0 - y_xmax) >= 0: best_xx = x0 + ( ytol + y_xmax - y0) / alpha
                        if (y0 - y_xmax) <  0: best_xx = x0 + (-ytol + y_xmax - y0) / alpha
                        #print("best_xx 2", best_xx)

                    self.ykey_best_xs[ykey][il] = best_xx

        # Here we change the x-y labels for the plots using an hard-coded mapping
        # in order to add additional info on units and normalization.
        auto_key_label = dict(
            ecut=r"$E_{cut}$ (Ha)",
            energy_per_atom=r"$E/N_{at}$ (eV)",
            pressure="P (GPa)",
        )

        for key, label in auto_key_label.items():
            self.set_label(key, label, ignore_exc=True)

    def set_label(self, key: str, label: str, ignore_exc=False) -> None:
        """
        Set the label for `key` to be used in the plot.
        Dont't raise exception if `ignore_exc` is True.
        """
        if key in self.ykey2label:
            self.ykey2label[key] = label
        elif key == self.xkey:
            self.xlabel = label
        else:
            if not ignore_exc:
                raise ValueError(
                    f"key:`{key}` should be either in {list(self.ykey2label.keys())} or {self.xname}")

    def get_ylabel(self, ykey: str) -> str:
        """Return the ylabel to be used for `ykey` in the plot."""
        return self.ykey2label[ykey]

    def ytol_ix_xx(self, ykey) -> Iterator[tuple]:
        """
        Iterate over (ytols, ixs, and xs) for the given ``ykey`.
        """
        return zip(self.ytols_dict[ykey], self.ykey_ixs[ykey] ,self.ykey_best_xs[ykey])

    def get_dataframe_ykey(self, ykey: str) -> pd.DataFrame:
        """Return dataframe with convergence params for `ykey`."""
        rows = []
        for ytol, ix, xx in self.ytol_ix_xx(ykey):
            rows.append(dict(ytol=ytol, ix=ix, xx=xx, ykey=ykey))
            #rows.append(dict(ytol=ytol, ix=ix, xx=xx), xx_best=xx_best, ykey=ykey)
        return pd.DataFrame(rows)

    def to_string(self, verbose=0) -> str:
        """
        String representation with verbosity level `verbose`.
        """
        lines = []; app = lines.append
        app(f"Number of points for x-axis: {len(self.xs)}")
        for ykey in self.yvals_dict:
            app("ykey: %s" % ykey)
            df = self.get_dataframe_ykey(ykey)
            app(str(df))

        return "\n".join(lines)

    def __str__(self) -> str:
        return self.to_string()

    def _decorate_ax(self, ax, ykey, ys, yscale) -> None:
        """
        Decorate axis ax by adding patches showing the convergence window
        and vertical lines where convergence is achieved.

        Args:
            ax: matplotlib axes.
            ykey: y-name
            ys: y-values
            yscale: "linear" or "log"
        """
        # Precompute y-limits of the converge window for each tolerance.
        y_xmax = ys[-1]
        ytols = self.ytols_dict[ykey]
        ntols = len(ytols)
        ylims = np.empty((ntols, 2))
        ylims_log = np.empty((ntols, 2))

        for il, ytol in enumerate(ytols):
            # Absolute tolerance.
            y0, y1 = y_xmax - ytol, y_xmax + ytol
            y1_log = ytol
            ylims[il] = [y0, y1]
            ylims_log[il] = [0, y1_log]

        # Loop again as ylimits are known.
        for il, ytol in enumerate(ytols):
            label = r"$|y-y_\infty| \leq %s$" % ytol
            span_style = dict(alpha=0.2, color=self.color_ilevel[il], zorder=abs(ytol), hatch=self.HATCH)

            y0, y1 = ylims[il]
            y0_log, y1_log = ylims_log[il]

            if il == ntols - 1:
                if yscale == "linear":
                    ax.axhspan(y0, y1, label=label, **span_style)
                elif yscale == "log":
                    ax.axhspan(y0_log, y1_log, label=label, **span_style)
                else:
                    raise ValueError(f"Invalid yscale: {yscale}")
            else:
                # Use limits of the next window to avoid overlapping patches.
                if yscale == "linear":
                    ax.axhspan(y0, ylims[il+1,0], label=label, **span_style)
                    ax.axhspan(ylims[il+1,1], y1, **span_style)
                elif yscale == "log":
                    ax.axhspan(y0_log, ylims_log[il+1,0], label=label, **span_style)
                    ax.axhspan(ylims_log[il+1,1], y1_log, **span_style)
                else:
                    raise ValueError(f"Invalid yscale: {yscale}")

            # Add vertical line to show best_xx.
            best_xx = self.ykey_best_xs[ykey][il]
            line_style = dict(lw=1, color=self.color_ilevel[il], ls=":")
            if best_xx is not None:
                ax.axvline(best_xx, **line_style)

    @add_fig_kwargs
    def plot(self, ax_mat=None, fontsize=8, **kwargs) -> Figure:
        """
        Plot convergence profile. A new grid is build if `ax_mat` is None:
        """
        nrows, ncols = len(self.yvals_dict), 2

        ax_mat, fig, plt = get_axarray_fig_plt(ax_mat, nrows=nrows, ncols=ncols,
                                               sharex=False, sharey=False, squeeze=False)

        # TODO
        #for icol in range(ncols):
        #    ax_share("x", ax_mat[0,icol], ax_mat[1,icol])

        for irow, ((ykey, ys), ax_row) in enumerate(zip(self.yvals_dict.items(), ax_mat)):
            # Plot y(x)
            ax1, ax2 = ax_row
            ax1.plot(self.xs, ys, marker="o", color="k")
            ax1.set_ylabel(self.get_ylabel(ykey))
            self._decorate_ax(ax1, ykey, ys, "linear")

            # Plot |y(x) - y_xmax| on log scale.
            abs_diffs = np.abs(ys - ys[-1])
            ax2.plot(self.xs, abs_diffs, marker="o", color="k")
            ax2.set_yscale("log")
            self._decorate_ax(ax2, ykey, ys, "log")
            ax2.set_xlim(self.xs[0] - 1, self.xs[-2] + 1)

            title = ""
            for i, (ytol, ix, xx) in enumerate(self.ytol_ix_xx(ykey)):
                pre_str = "" if i == 0 else ", "
                ytol_string = str(ytol)
                #print("ytol_string:", ytol_string, "pre_str:", pre_str, "ytol_string:", ytol_string, "xx:", xx)
                if xx is not None:
                    s = "x: %.1f for $\Delta$: %s" % (xx, ytol_string)
                else:
                    s = "x: ?? for $\Delta$: %s" % (ytol_string)
                title += pre_str + s

            ax2.set_title(title, fontsize=fontsize)
            ax2.set_ylabel(r"$|y-y(x_{max})|$", fontsize=fontsize)

            set_grid_legend(ax_row, fontsize,
                            xlabel=self.xlabel if irow == (nrows - 1) else None,
                            grid=False, legend=True)

        fig.tight_layout()
        return fig


class ArrayPlotter:

    def __init__(self, *labels_and_arrays):
        """
        Args:
            labels_and_arrays: list [("label1", arr1), ("label2", arr2)]
        """
        self._arr_dict = {}
        for label, array in labels_and_arrays:
            self.add_array(label, array)

    def __len__(self) -> int:
        return len(self._arr_dict)

    def __iter__(self):
        return self._arr_dict.__iter__()

    def keys(self):
        return self._arr_dict.keys()

    def items(self):
        return self._arr_dict.items()

    def add_array(self, label: str, array) -> None:
        """Add array with the given name."""
        if label in self._arr_dict:
            raise ValueError("%s is already in %s" % (label, list(self._arr_dict.keys())))

        self._arr_dict[label] = array

    def add_arrays(self, labels: list, arr_list: list) -> None:
        """
        Add a list of arrays

        Args:
            labels: List of labels.
            arr_list: List of arrays.
        """
        assert len(labels) == len(arr_list)
        for label, arr in zip(labels, arr_list):
            self.add_array(label, arr)

    @add_fig_kwargs
    def plot(self, cplx_mode="abs", colormap="jet", fontsize=8, **kwargs) -> Figure:
        """
        Args:
            cplx_mode: "abs" for absolute value, "re", "im", "angle"
            colormap: matplotlib colormap.
            fontsize: legend and label fontsize.
        """
        # Build grid of plots.
        num_plots, ncols, nrows = len(self), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + (num_plots % ncols)

        import matplotlib.pyplot as plt
        fig, ax_mat = plt.subplots(nrows=nrows, ncols=ncols, sharex=False, sharey=False, squeeze=False)
        # Don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_mat[-1, -1].axis("off")

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        from matplotlib.ticker import MultipleLocator

        for ax, (label, arr) in zip(ax_mat.flat, self.items()):
            data = data_from_cplx_mode(cplx_mode, arr)
            # Use origin to place the [0, 0] index of the array in the lower left corner of the axes.
            img = ax.matshow(data, interpolation='nearest', cmap=colormap, origin='lower', aspect="auto")
            ax.set_title("(%s) %s" % (cplx_mode, label), fontsize=fontsize)

            # Make a color bar for this ax
            # Create divider for existing axes instance
            # http://stackoverflow.com/questions/18266642/multiple-imshow-subplots-each-with-colorbar
            divider3 = make_axes_locatable(ax)
            # Append axes to the right of ax, with 10% width of ax
            cax3 = divider3.append_axes("right", size="10%", pad=0.05)
            # Create colorbar in the appended axes
            # Tick locations can be set with the kwarg `ticks`
            # and the format of the ticklabels with kwarg `format`
            cbar3 = plt.colorbar(img, cax=cax3, ticks=MultipleLocator(0.2), format="%.2f")
            # Remove xticks from ax
            ax.xaxis.set_visible(False)
            # Manually set ticklocations
            #ax.set_yticks([0.0, 2.5, 3.14, 4.0, 5.2, 7.0])

            # Set grid
            ax.grid(True, color='white')

        fig.tight_layout()
        return fig


#TODO use object and introduce c for color, client code should be able to customize it.
# Rename it to ScatterData
class Marker(namedtuple("Marker", "x y s")):
    """
    Stores the position and the size of the marker.
    A marker is a list of tuple(x, y, s) where x, and y are the position
    in the graph and s is the size of the marker.
    Used for plotting purpose e.g. QP data, energy derivatives...

    Example::

        x, y, s = [1, 2, 3], [4, 5, 6], [0.1, 0.2, -0.3]
        marker = Marker(x, y, s)
        marker.extend((x, y, s))

    """
    def __new__(cls, *xys):
        """Extends the base class adding consistency check."""
        if not xys:
            xys = ([], [], [])
            return super().__new__(cls, *xys)

        if len(xys) != 3:
            raise TypeError("Expecting 3 entries in xys got %d" % len(xys))

        x = np.asarray(xys[0])
        y = np.asarray(xys[1])
        s = np.asarray(xys[2])
        xys = (x, y, s)

        for s in xys[-1]:
            if np.iscomplex(s):
                raise ValueError("Found ambiguous complex entry %s" % str(s))

        return super().__new__(cls, *xys)

    def __bool__(self):
        return bool(len(self.s))

    __nonzero__ = __bool__

    def extend(self, xys):
        """
        Extend the marker values.
        """
        if len(xys) != 3:
            raise TypeError("Expecting 3 entries in xys got %d" % len(xys))

        self.x.extend(xys[0])
        self.y.extend(xys[1])
        self.s.extend(xys[2])

        lens = np.array((len(self.x), len(self.y), len(self.s)))
        if np.any(lens != lens[0]):
            raise TypeError("x, y, s vectors should have same lengths but got %s" % str(lens))

    def posneg_marker(self) -> tuple[Marker, Marker]:
        """
        Split data into two sets: the first one contains all the points with positive size.
        The first set contains all the points with negative size.
        """
        pos_x, pos_y, pos_s = [], [], []
        neg_x, neg_y, neg_s = [], [], []

        for x, y, s in zip(self.x, self.y, self.s):
            if s >= 0.0:
                pos_x.append(x)
                pos_y.append(y)
                pos_s.append(s)
            else:
                neg_x.append(x)
                neg_y.append(y)
                neg_s.append(s)

        return self.__class__(pos_x, pos_y, pos_s), self.__class__(neg_x, neg_y, neg_s)


class Exposer:
    """
    Base class for Exposer objects.

    Example:

        kws = dict(show=False)
        with Exposer.as_exposer("panel") as e:
            e(obj.plot1(**plot_kws))
            e(obj.plot2(**plot_kws))
    """

    @classmethod
    def as_exposer(cls, exposer, **kwargs) -> Exposer:
        """
        Return an instance of Exposer, usually from a string with then name.

        Args:
            exposer: "mpl" for MplExposer, "panel" for PanelExposer.
        """
        if isinstance(exposer, cls): return exposer

        # Assume string.
        exposer_cls = dict(
            mpl=MplExposer,
            panel=PanelExposer,
        )[exposer]
        return exposer_cls(**kwargs)

    def add_obj_with_yield_figs(self, obj: Any) -> None:
        """
        Add an object implementing a `yield_figs` method to the Exposer.
        """
        if not hasattr(obj, "yield_figs"):
            raise TypeError(f"object of type {type(obj)} does not implement `yield_figs` method")

        for fig in obj.yield_figs():
            self.add_fig(fig)

    def __call__(self, obj: Any):
        """
        Add an object to the Exposer
        Support mpl figure, list of figures or generator yielding figures.
        """
        import types
        if isinstance(obj, (types.GeneratorType, list, tuple)):
            for fig in obj:
                self.add_fig(fig)
        else:
            self.add_fig(obj)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement. """
        if exc_type is not None: return
        self.expose()


class MplExposer(Exposer): # pragma: no cover
    """
    Context manager used to produce several matplotlib figures and show
    all of them at once so that users do not have to close the window
    to visualize the next one.

    Example:

        plot_args = dict(show=False)
        with MplExposer() as e:
            e(obj.plot1(**plot_args))
            e(obj.plot2(**plot_args))
    """

    def __init__(self, slide_mode=False, slide_timeout=None, verbose=1, **kwargs):
        """
        Args:
            slide_mode: If True, iterate over figures. Default: Expose all figures at once.
            slide_timeout: Close figure after slide-timeout seconds. Block if None.
            verbose: verbosity level
        """
        self.figures = []
        self.slide_mode = bool(slide_mode)
        self.timeout_ms = slide_timeout
        self.verbose = verbose
        if self.timeout_ms is not None:
            self.timeout_ms = int(self.timeout_ms * 1000)
            assert self.timeout_ms >= 0

        if self.verbose:
            if self.slide_mode:
                print("\nSliding matplotlib figures with slide timeout: %s [s]" % slide_timeout)
            else:
                print("\nLoading all matplotlib figures before showing them. It may take some time...")

        self.start_time = time.time()

    def add_fig(self, fig: Figure) -> None:
        """
        Add a matplotlib figure.
        """
        if fig is None: return

        if not self.slide_mode:
            self.figures.append(fig)
        else:
            import matplotlib.pyplot as plt
            if self.timeout_ms is not None:
                # Creating a timer object
                # timer calls plt.close after interval milliseconds to close the window.
                timer = fig.canvas.new_timer(interval=self.timeout_ms)
                timer.add_callback(plt.close, fig)
                timer.start()

            plt.show()
            fig.clear()

    def expose(self) -> None:
        """
        Show all figures. Clear figures if needed.
        """
        if not self.slide_mode:
            print("All figures in memory, elapsed time: %.3f s" % (time.time() - self.start_time))
            import matplotlib.pyplot as plt
            plt.show()
            for fig in self.figures:
                fig.clear()


class PanelExposer(Exposer):  # pragma: no cover
    """
    Context manager used to produce several matplotlib/plotly figures
    and show all of them inside the web browser using a panel template.

    Example:

        with PanelExposer() as e:
            e(obj.plot1(show=False))
            e(obj.plot2(show=False))
    """
    def __init__(self, title=None, dpi=92, verbose=1, **kwargs):
        """
        Args:
            title: String to be show in the header.
            verbose: verbosity level
        """
        self.title = title
        self.figures = []
        self.verbose = verbose
        self.dpi = int(dpi)

        if self.verbose:
            print("\nLoading all figures before showing them. It may take some time...")

        self.start_time = time.time()

    def add_fig(self, fig: Figure) -> None:
        """Add a matplotlib figure."""
        if fig is None: return
        self.figures.append(fig)

    def expose(self):
        """Show all figures. Clear figures if needed."""
        import panel as pn
        pn.config.sizing_mode = 'stretch_width'
        from abipy.panels.core import get_template_cls_from_name
        cls = get_template_cls_from_name("FastGridTemplate")

        template = cls(
            title=self.title if self.title is not None else self.__class__.__name__,
            header_background="#ff8c00 ", # Dark orange
        )
        #pn.config.sizing_mode = 'stretch_width'
        from abipy.panels.core import mpl, ply
        for i, fig in enumerate(self.figures):
            row, col = divmod(i, 2)
            if is_plotly_figure(fig):
                p = ply(fig, with_divider=False)
            elif is_mpl_figure(fig):
                p = mpl(fig, with_divider=False, dpi=self.dpi)
            else:
                raise TypeError(f"Don't know how to handle type: `{type(fig)}`")

            if hasattr(template.main, "append"):
                template.main.append(p)
            else:
                # Assume .main area acts like a GridSpec
                row_slice = slice(3 * row, 3 * (row + 1))
                if col == 0: template.main[row_slice, :6] = p
                if col == 1: template.main[row_slice, 6:] = p

        return template.show()


def plot_unit_cell(lattice, ax=None, **kwargs) -> tuple[Figure, Axes]:
    """
    Adds the unit cell of the lattice to a matplotlib Axes3D

    Args:
        lattice: Lattice object
        ax: matplotlib :class:`Axes3D` or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to black
            and linewidth to 3.

    Returns:
        matplotlib figure and ax
    """
    ax, fig, plt = get_ax3d_fig_plt(ax)

    if "color" not in kwargs: kwargs["color"] = "k"
    if "linewidth" not in kwargs: kwargs["linewidth"] = 3

    v = 8 * [None]
    v[0] = lattice.get_cartesian_coords([0.0, 0.0, 0.0])
    v[1] = lattice.get_cartesian_coords([1.0, 0.0, 0.0])
    v[2] = lattice.get_cartesian_coords([1.0, 1.0, 0.0])
    v[3] = lattice.get_cartesian_coords([0.0, 1.0, 0.0])
    v[4] = lattice.get_cartesian_coords([0.0, 1.0, 1.0])
    v[5] = lattice.get_cartesian_coords([1.0, 1.0, 1.0])
    v[6] = lattice.get_cartesian_coords([1.0, 0.0, 1.0])
    v[7] = lattice.get_cartesian_coords([0.0, 0.0, 1.0])

    for i, j in ((0, 1), (1, 2), (2, 3), (0, 3), (3, 4), (4, 5), (5, 6),
                 (6, 7), (7, 4), (0, 7), (1, 6), (2, 5), (3, 4)):
        ax.plot(*zip(v[i], v[j]), **kwargs)

    # Plot cartesian frame
    ax_add_cartesian_frame(ax)

    return fig, ax


def ax_add_cartesian_frame(ax, start=(0, 0, 0)) -> Axes:
    """
    Add cartesian frame to 3d axis at point `start`.
    """
    # https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import proj3d
    arrow_opts = {"color": "k"}
    arrow_opts.update(dict(lw=1, arrowstyle="-|>",))

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            FancyArrowPatch.draw(self, renderer)

    start = np.array(start)
    for end in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        end = start + np.array(end)
        xs, ys, zs = list(zip(start, end))
        p = Arrow3D(xs, ys, zs,
                   connectionstyle='arc3', mutation_scale=20,
                   alpha=0.8, **arrow_opts)
        ax.add_artist(p)

    return ax


def plot_structure(structure,
                   ax=None, to_unit_cell=False, alpha=0.7,
                   style="points+labels", color_scheme="VESTA", **kwargs) -> Figure:
    """
    Plot structure with matplotlib (minimalistic version).

    Args:
        structure: |Structure| object
        ax: matplotlib :class:`Axes3D` or None if a new figure should be created.
        alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
        to_unit_cell: True if sites should be wrapped into the first unit cell.
        style: "points+labels" to show atoms sites with labels.
        color_scheme: color scheme for atom types. Allowed values in ("Jmol", "VESTA")

    Returns: |matplotlib-Figure|
    """
    fig, ax = plot_unit_cell(structure.lattice, ax=ax, linewidth=1)

    from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
    from pymatgen.vis.structure_vtk import EL_COLORS
    xyzs, colors = np.empty((len(structure), 4)), []

    for i, site in enumerate(structure):
        symbol = site.specie.symbol
        color = tuple(i / 255 for i in EL_COLORS[color_scheme][symbol])
        radius = CovalentRadius.radius[symbol]
        if to_unit_cell and hasattr(site, "to_unit_cell"): site = site.to_unit_cell()
        # Use cartesian coordinates.
        x, y, z = site.coords
        xyzs[i] = (x, y, z, radius)
        colors.append(color)
        if "labels" in style:
            ax.text(x, y, z, symbol)

    # The definition of sizes is not optimal because matplotlib uses points
    # wherease we would like something that depends on the radius (5000 seems to give reasonable plots)
    # For possibile approaches, see
    # https://stackoverflow.com/questions/9081553/python-scatter-plot-size-and-style-of-the-marker/24567352#24567352
    # https://gist.github.com/syrte/592a062c562cd2a98a83
    if "points" in style:
        x, y, z, s = xyzs.T.copy()
        s = 5000 * s ** 2
        ax.scatter(x, y, zs=z, s=s, c=colors, alpha=alpha)  #facecolors="white", #edgecolors="blue"

    ax.set_title(structure.composition.formula)
    ax.set_axis_off()

    return fig


def _generic_parser_fh(fh) -> dict:
    """
    Parse file with data in tabular format. Supports multi datasets a la gnuplot.
    Mainly used for files without any schema, not even CSV

    Args:
        fh: File object

    Returns:
        dict title --> numpy array
        where title is taken from the first (non-empty) line preceding the dataset
    """
    arr_list = [None]
    data = []
    head_list = []
    count = -1
    last_header = None
    for l in fh:
        l = l.strip()
        if not l or l.startswith("#"):
            count = -1
            last_header = l
            if arr_list[-1] is not None: arr_list.append(None)
            continue

        count += 1
        if count == 0: head_list.append(last_header)
        if arr_list[-1] is None: arr_list[-1] = []
        data = arr_list[-1]
        data.append(list(map(float, l.split())))

    if len(head_list) != len(arr_list):
        raise RuntimeError("len(head_list) != len(arr_list), %d != %d" % (len(head_list), len(arr_list)))

    od = {}
    for key, data in zip(head_list, arr_list):
        key = " ".join(key.split())
        if key in od:
            print("Header %s already in dictionary. Using new key %s" % (key, 2 * key))
            key = 2 * key
        od[key] = np.array(data).T.copy()

    return od


class GenericDataFilePlotter:
    """
    Extract data from a generic text file with results in tabular format and plot data with matplotlib.
    Multiple datasets are supported.
    No attempt is made to handle metadata (e.g. column name)
    Mainly used to handle text files written without any schema.
    """
    def __init__(self, filepath: str):
        with open(filepath, "rt") as fh:
            self.od = _generic_parser_fh(fh)

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation with verbosity level `verbose`."""
        lines = []
        for key, arr in self.od.items():
            lines.append("key: `%s` --> array shape: %s" % (key, str(arr.shape)))
        return "\n".join(lines)

    @add_fig_kwargs
    def plot(self, use_index=False, fontsize=8, **kwargs) -> Figure:
        """
        Plot all arrays. Use multiple axes if datasets.

        Args:
            use_index: By default, the x-values are taken from the first column.
                If use_index is False, the x-values are the row index.
            fontsize: fontsize for title.
            kwargs: options passed to ``ax.plot``.

        Return: |matplotlib-figure|
        """
        # build grid of plots.
        num_plots, ncols, nrows = len(self.od), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # Don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        for ax, (key, arr) in zip(ax_list, self.od.items()):
            ax.set_title(key, fontsize=fontsize)
            ax.grid(True)
            xs = arr[0] if not use_index else list(range(len(arr[0])))
            for ys in arr[1:] if not use_index else arr:
                ax.plot(xs, ys)

        return fig


class GenericDataFilesPlotter:

    @classmethod
    def from_files(cls, filepaths: list[str]) -> GenericDataFilesPlotter:
        """
        Build object from a list of `filenames`.
        """
        new = cls()
        for filepath in filepaths:
            new.add_file(filepath)
        return new

    def __init__(self):
        self.odlist = []
        self.filepaths = []

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        lines = []
        app = lines.append
        for od, filepath in zip(self.odlist, self.filepaths):
            app("File: %s" % filepath)
            for key, arr in od.items():
                lines.append("\tkey: `%s` --> array shape: %s" % (key, str(arr.shape)))

        return "\n".join(lines)

    def add_file(self, filepath: str) -> None:
        """Add data from `filepath`"""
        with open(filepath, "rt") as fh:
            self.odlist.append(_generic_parser_fh(fh))
            self.filepaths.append(filepath)

    @add_fig_kwargs
    def plot(self, use_index=False, fontsize=8, colormap="viridis", **kwargs) -> Figure:
        """
        Plot all arrays. Use multiple axes if datasets.

        Args:
            use_index: By default, the x-values are taken from the first column.
                If use_index is False, the x-values are the row index.
            fontsize: fontsize for title.
            colormap: matplotlib color map.
            kwargs: options passed to ``ax.plot``.

        Return: |matplotlib-figure|
        """
        if not self.odlist: return None

        # Compute intersection of all keys.
        # Here we loose the initial ordering in the dict but oh well!
        klist = [list(d.keys()) for d in self.odlist]
        keys = set(klist[0]).intersection(*klist)
        if not keys:
            print("Warning: cannot find common keys in files. Check input data")
            return None

        # Build grid of plots.
        num_plots, ncols, nrows = len(keys), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots // ncols) + (num_plots % ncols)

        ax_list, fig, plt = get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                                sharex=False, sharey=False, squeeze=False)
        ax_list = ax_list.ravel()

        # Don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: ax_list[-1].axis("off")

        cmap = plt.get_cmap(colormap)
        line_cycle = itertools.cycle(["-", ":", "--", "-.",])

        # One ax for key, each ax may show multiple arrays
        # so we need different line styles that are consistent with input data.
        # Figure may be crowded but it's difficult to do better without metadata
        # so I'm not gonna spend time to implement more complicated logic.
        for ax, key in zip(ax_list, keys):
            ax.set_title(key, fontsize=fontsize)
            ax.grid(True)
            for iod, (od, filepath) in enumerate(zip(self.odlist, self.filepaths)):
                if key not in od: continue
                arr = od[key]
                color = cmap(iod / len(self.odlist))
                xvals = arr[0] if not use_index else list(range(len(arr[0])))
                arr_list = arr[1:] if not use_index else arr
                for iarr, (ys, linestyle) in enumerate(zip(arr_list, line_cycle)):
                    ax.plot(xvals, ys, color=color, linestyle=linestyle,
                            label=os.path.relpath(filepath) if iarr == 0 else None)

            ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig


##########################
# Plotly helper functions
##########################

_LATEX_GREEK_TO_UNICODE = dict(
    alpha="α",
    beta="β",
    gamma="ɣ",
    delta="δ",
    epsilon="ε",
    zeta="ζ",
    eta="η",
    theta="θ",
    iota="ι",
    kappa="κ",
    #lambda="λ",
    mu="μ",
    nu="ν",
    xi="ξ",
    omicron="ο",
    pi="π",
    rho="ρ",
    sigma="σ",
    tau="τ",
    upsilon="υ",
    phi="φ",
    chi="χ",
    psi="ψ",
    omega="ω",
    # Capital case:
    Alpha="Α",
    Beta="Β",
    Gamma="Γ",
    Delta="Δ",
    Epsilon="Ε",
    Zeta="Ζ",
    Eta="Η",
    Theta="Θ",
    Iota="Ι",
    Kappa="Κ",
    Lambda="Λ",
    Mu="Μ",
    Nu="Ν",
    Xi="Ξ",
    Omicron="Ο",
    Po="Π",
    Rho="Ρ",
    Sigma="Σ",
    Tau="Τ",
    Upsilon="Υ",
    Phi="Φ",
    Chi="Χ",
    Psi="Ψ",
    Omega="Ω",
)

_LATEX_GREEK_TO_UNICODE["lambda"] = "λ"


def latex_greek_2unicode(latex: str) -> str:
    """
    Convert a single greek letter in latex notation into unicode
    """
    s = latex.replace("$", "").replace("\\", "").strip()
    return _LATEX_GREEK_TO_UNICODE[s]


def is_plotly_figure(obj: Any) -> bool:
    """Return True if obj is a plotly Figure."""
    import plotly.graph_objs as go
    return isinstance(obj, go.Figure)
    #return isinstance(obj, (go.Figure, go.FigureWidget))


class PlotlyRowColDesc:
    """
    This object specifies the position of a plotly subplot inside a grid.

    rcd: PlotlyRowColDesc object used when fig is not None to specify the (row, col) of the subplot in the grid.
    """

    @classmethod
    def from_object(cls, obj: Any) -> PlotlyRowColDesc:
        """
        Build an instance for a generic object.
        If oject is None, a simple descriptor corresponding to a (1,1) grid is returned.
        """
        if obj is None: return cls(0, 0, 1, 1)
        if isinstance(obj, cls): return obj

        # Assume list with 4 integers
        try:
            return cls(*obj)
        except Exception as exc:
            raise TypeError(f"Dont know how to convert `{type(obj)}` into `{cls}`")

    def __init__(self, py_row: int, py_col: int, nrows: int, ncols: int):
        """
        Args:
            py_row, py_col: python index of the subplot in the grid (starts from 0)
            nrows, ncols: Number of rows/cols in the grid.
        """
        self.py_row, self.py_col = (py_row, py_col)
        self.nrows, self.ncols = (nrows, ncols)
        self.iax = 1 + self.py_col + self.py_row * self.ncols
        # Note that plotly col and row start from 1.
        if nrows == 1 and ncols == 1:
            self.ply_row, self.ply_col = (None, None)
        else:
            self.ply_row, self.ply_col = (self.py_row + 1, self.py_col + 1)

    def __str__(self) -> str:
        lines = []
        app = lines.append
        app("py_rowcol: (%d, %d) in grid: (%d, %d)" % (self.py_row, self.py_col, self.nrows, self.ncols))
        app("plotly_rowcol: (%s, %s)" % (self.ply_row, self.ply_col))

        return "\n".join(lines)

    #@lazy_property
    #def rowcol_dict(self):
    #    if self.nrows == 1 and self.ncols == 1: return {}
    #    return dict(row=self.ply_row, col=self.ply_col)


def get_figs_plotly(nrows=1, ncols=1, subplot_titles=(), sharex=False, sharey=False, **fig_kw):
    """
    Helper function used in plot functions that build the `plotly` figure by calling plotly.subplots.

    Returns:
        figure: plotly graph_objects figure
        go: plotly graph_objects module.
    """
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go

    fig = make_subplots(rows=int(nrows), cols=int(ncols), subplot_titles=subplot_titles, shared_xaxes=sharex,
                        shared_yaxes=sharey, **fig_kw)

    return fig, go


def get_fig_plotly(fig=None, **fig_kw):
    """
    Helper function used in plot functions that build the `plotly` figure by calling
    plotly.graph_objects.Figure if fig is None else return fig

    Returns:
        figure: plotly graph_objects figure
        go: plotly graph_objects module.
    """
    import plotly.graph_objects as go

    if fig is None:
        fig = go.Figure(**fig_kw)
        #fig = go.FigureWidget(**fig_kw)

    return fig, go


def plotly_set_lims(fig, lims, axname, iax=None) -> tuple:
    """
    Set the data limits for the axis ax.

    Args:
        fig: Plotly Figure.
        lims: tuple(2) for (left, right), if tuple(1) or scalar for left only, none is set.
        axname: "x" for x-axis, "y" for y-axis.
        iax: An int, use iax=n to decorate the nth axis when the fig has subplots.

    Return: (left, right)
    """
    left, right = None, None
    if lims is None: return (left, right)

    # iax = kwargs.pop("iax", 1)
    # xaxis = 'xaxis%u' % iax
    #fig.layout[xaxis].title.text = "Wave Vector"

    axis = dict(x=fig.layout.xaxis, y=fig.layout.yaxis)[axname]

    len_lims = None
    try:
        len_lims = len(lims)
    except TypeError:
        # Assume Scalar
        left = float(lims)

    if len_lims is not None:
        if len(lims) == 2:
            left, right = lims[0], lims[1]
        elif len(lims) == 1:
            left = lims[0]

    ax_range = axis.range
    if ax_range is None and (left is None or right is None):
        return None, None

    #if left is not None: ax_range[0] = left
    #if right is not None: ax_range[1] = right

    # Example: fig.update_layout(yaxis_range=[-4,4])
    k = dict(x="xaxis", y="yaxis")[axname]
    if iax:
        k= k + str(iax)
    fig.layout[k].range = [left, right]

    return left, right


_PLOTLY_DEFAULT_SHOW = [True]


def set_plotly_default_show(true_or_false: bool) -> None:
    """
    Set the default value of show in the add_plotly_fig_kwargs decorator.
    Useful for instance when generating the sphinx gallery of plotly plots.
    """
    _PLOTLY_DEFAULT_SHOW[0] = true_or_false


def add_plotly_fig_kwargs(func: Callable) -> Callable:
    """
    Decorator that adds keyword arguments for functions returning plotly figures.
    The function should return either a plotly figure or None to signal some
    sort of error/unexpected event.
    See doc string below for the list of supported options.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        # pop the kwds used by the decorator.
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", _PLOTLY_DEFAULT_SHOW[0])
        hovermode = kwargs.pop("hovermode", False)
        savefig = kwargs.pop("savefig", None)
        write_json = kwargs.pop("write_json", None)
        config = kwargs.pop("config", None)
        renderer = kwargs.pop("renderer", None)
        chart_studio = kwargs.pop("chart_studio", False)
        template = kwargs.pop("template", None)

        # Allow users to specify the renderer via shell env.
        if renderer is not None and os.getenv("PLOTLY_RENDERER", default=None) is not None:
            renderer = None

        # Call func and return immediately if None is returned.
        fig = func(*args, **kwargs)
        if fig is None:
            return fig

        # Operate on plotly figure.
        if title is not None:
            fig.update_layout(title_text=title, title_x=0.5)

        if template is not None:
            fig.update_layout(template=template)

        if savefig:
            # https://plotly.github.io/plotly.py-docs/generated/plotly.io.write_image.html
            if savefig.endswith("html"):
                from plotly.offline import plot as show_plotly
                show_plotly(fig, include_mathjax="cdn", filename=savefig, auto_open=False)

            else:
                try:
                    import kaleido
                except ImportError:
                    kaleido = False

                if kaleido is None:
                    raise ValueError(
                        "kaleido package required to save static ploty images\n"
                        "please install it using:\npip install kaleido"
                    )

                fig.write_image(savefig, engine="kaleido", scale=5, width=750, height=750)
                #fig.write_image(savefig)

        if write_json:
            import plotly.io as pio
            pio.write_json(fig, write_json)

        fig.layout.hovermode = hovermode

        if show: # and _PLOTLY_DEFAULT_SHOW:
            my_config = dict(
                responsive=True,
                #showEditInChartStudio=True,
                showLink=True,
                plotlyServerURL="https://chart-studio.plotly.com",
            )

            if config is not None:
                my_config.update(config)

            #add_template_buttons(fig)

            fig.show(renderer=renderer, config=my_config)

        if chart_studio:
            push_to_chart_studio(fig)

        return fig

    # Add docstring to the decorated method.
    s = (
            "\n\n"
            + """\
                Keyword arguments controlling the display of the figure:
                ================  ====================================================================
                kwargs            Meaning
                ================  ====================================================================
                title             Title of the plot (Default: None).
                show              True to show the figure (default: True).
                hovermode         True to show the hover info (default: False)
                savefig           "abc.png" , "abc.jpeg" or "abc.webp" to save the figure to a file.
                write_json        Write plotly figure to `write_json` JSON file.
                                  Inside jupyter-lab, one can right-click the `write_json` file from
                                  the file menu and open with "Plotly Editor".
                                  Make some changes to the figure, then use the file menu to save
                                  the customized plotly plot.
                                  Requires `jupyter labextension install jupyterlab-chart-editor`.
                                  See https://github.com/plotly/jupyterlab-chart-editor
                renderer          (str or None (default None)) –
                                  A string containing the names of one or more registered renderers
                                  (separated by ‘+’ characters) or None. If None, then the default
                                  renderers specified in plotly.io.renderers.default are used.
                                  See https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html
                config (dict)     A dict of parameters to configure the figure. The defaults are set in plotly.js.
                chart_studio      True to push figure to chart_studio server. Requires authenticatios.
                                  Default: False.
                template          Plotly template. See https://plotly.com/python/templates/
                                  ["plotly", "plotly_white", "plotly_dark", "ggplot2",
                                   "seaborn", "simple_white", "none"]
                                  Default is None that is the default template is used.
                ================  ====================================================================
        """
    )

    if wrapper.__doc__ is not None:
        # Add s at the end of the docstring.
        wrapper.__doc__ += "\n" + s
    else:
        # Use s
        wrapper.__doc__ = s

    return wrapper


def plotlyfigs_to_browser(figs, filename=None, browser=None):
    """
    Save a list of plotly figures in an HTML file and open it the browser.
    Useful to display multiple figures generated by different AbiPy methods
    without having to construct a plotly subplot grid.

    Args:
        figs: List of plotly figures.
        filename: File name to save in. Use temporary filename if filename is None.
        browser: Open webpage in ``browser``. Use $BROWSER if None.

    Example:

        fig1 = plotter.combiplotly(renderer="browser", title="foo", show=False)
        fig2 = plotter.combiplotly(renderer="browser", title="bar", show=False)
        from abipy.tools.plotting import plotlyfigs_to_browser
        plotlyfigs_to_browser([fig1, fig2])

    Return: path to HTML file.
    """
    if filename is None:
        import tempfile
        fd, filename = tempfile.mkstemp(text=True, suffix=".html")

    if not isinstance(figs, (list, tuple)): figs = [figs]

    # Based on https://stackoverflow.com/questions/46821554/multiple-plotly-plots-on-1-page-without-subplot
    with open(filename, "wt") as fp:
        for i, fig in enumerate(figs):
            first = True if i == 0 else False
            fig.write_html(fp, include_plotlyjs=first, include_mathjax="cdn" if first else False)

    import webbrowser
    print("Opening HTML file:", filename)
    webbrowser.get(browser).open_new_tab("file://" + filename)

    return filename


def plotly_klabels(labels: list, allow_dupes=False) -> list:
    """
    This helper function polish a list of k-points labels before calling plotly by:

        - Checking if we have two equivalent consecutive labels (only the first one is shown and the second one is set to "")
        - Replacing particular Latex tokens with unicode as plotly support for Latex is far from optimal.

    Return: New list labels, same length as input labels.
    """
    new_labels = labels.copy()

    if not allow_dupes:
        # Don't show label if previous k-point is the same.
        for il in range(1, len(new_labels)):
            if new_labels[il] == new_labels[il - 1]: new_labels[il] = ""

    replace = {
        r"$\Gamma$": "Γ",
    }

    for il in range(len(new_labels)):
        if new_labels[il] in replace:
            new_labels[il] = replace[new_labels[il]]

    return new_labels


def plotly_set_xylabels(fig, xlabel, ylabel, exchange_xy):
    """
    Set the x- and the y-label of axis ax, exchanging x and y if exchange_xy
    """
    if exchange_xy: xlabel, ylabel = ylabel, xlabel
    fig.layout.xaxis.title.text = xlabel
    fig.layout.yaxis.title.text = ylabel


_PLOTLY_AUTHEHTICATED = False


def plotly_chartstudio_authenticate():
    """
    Authenticate the user on the chart studio portal by reading `PLOTLY_USERNAME` and `PLOTLY_API_KEY`
    from the pymatgen configuration file located in $HOME/.pmgrc.yaml.

        PLOTLY_USERNAME: johndoe
        PLOTLY_API_KEY: XXXXXXXXXXXXXXXXXXXX

    """
    global _PLOTLY_AUTHEHTICATED
    if _PLOTLY_AUTHEHTICATED: return

    try:
        from pymatgen.core import SETTINGS
        #from pymatgen.settings import SETTINGS
    except ImportError:
        from pymatgen import SETTINGS

    example = """
Add it to $HOME/.pmgrc.yaml using the follow syntax:

PLOTLY_USERNAME: john_doe
PLOTLY_API_KEY: secret  # to get your api_key go to profile > settings > regenerate key

"""

    username = SETTINGS.get("PLOTLY_USERNAME")
    if username is None:
        raise RuntimeError(f"Cannot find PLOTLY_USERNAME in pymatgen settings.\n{example}")

    api_key = SETTINGS.get("PLOTLY_API_KEY")
    if api_key is None:
        raise RuntimeError(f"Cannot find PLOTLY_API_KEY in pymatgen settings.\n{example}")

    import chart_studio
    # https://towardsdatascience.com/how-to-create-a-plotly-visualization-and-embed-it-on-websites-517c1a78568b
    chart_studio.tools.set_credentials_file(username=username, api_key=api_key)
    _PLOTLY_AUTHEHTICATED = True


def push_to_chart_studio(figs) -> None:
    """
    Push a plotly figure or a list of figures to the chart studio cloud.
    """
    plotly_chartstudio_authenticate()
    import chart_studio.plotly as py
    if not isinstance(figs, (list, tuple)): figs = [figs]
    for fig in figs:
        py.plot(fig, auto_open=True)


####################################################
# This code is shamelessy taken from Adam's package
####################################################


def go_points(points, size=4, color="black", labels=None, **kwargs):

    #textposition = 'top right',
    #textfont = dict(color='#E58606'),
    mode = "markers" if labels is None else "markers+text"
    #text = labels

    if labels is not None:
        labels = plotly_klabels(labels, allow_dupes=True)

    import plotly.graph_objects as go
    return go.Scatter3d(
        x=[v[0] for v in points],
        y=[v[1] for v in points],
        z=[v[2] for v in points],
        marker=dict(size=size, color=color),
        mode=mode,
        text=labels,
        **kwargs
    )


def _add_if_not_in(d, key, value):
    if key not in d:
        d[key] = value


def go_line(v1, v2, color="black", width=2, mode="lines", **kwargs):

    _add_if_not_in(kwargs, "line_color", "black")
    _add_if_not_in(kwargs, "line_width", 2)

    import plotly.graph_objects as go
    return go.Scatter3d(
        mode=mode,
        x=[v1[0], v2[0]],
        y=[v1[1], v2[1]],
        z=[v1[2], v2[2]],
        #line=dict(color=color, width=width),
        **kwargs
    )


def go_lines(V, name=None, color="black", width=2, **kwargs):
    import plotly.graph_objects as go
    gen = ((v1, v2) for (v1, v2) in V)
    v1, v2 = next(gen)
    out = [
        go_line(v1, v2, width=width, color=color, name=name, legendgroup=name, **kwargs)
    ]
    out.extend(
        go_line(
            v1,
            v2,
            width=width,
            color=color,
            showlegend=False,
            legendgroup=name,
            **kwargs
        )
        for (v1, v2) in gen
    )
    return out


def vectors(lattice, name=None, color="black", width=4, **kwargs):
    gen = zip(lattice, ["a", "b", "c"])
    v, label = next(gen)

    out = [
        go_line(
            [0, 0, 0],
            v,
            text=["", label],
            width=width,
            color=color,
            name=name,
            legendgroup=name,
            mode="lines+text",
            **kwargs
        )
    ]
    out.extend(
        go_line(
            [0, 0, 0],
            v,
            text=["", label],
            width=width,
            color=color,
            showlegend=False,
            legendgroup=name,
            mode="lines+text",
            **kwargs
        )
        for (v, label) in gen
    )
    return out


def get_vectors(lattice_mat, name=None, color="black", width=2, **kwargs):
    return go_lines([[[0, 0, 0], v] for v in lattice_mat], **kwargs)


def get_box(lattice_mat, **kwargs):
    a, b, c = lattice_mat
    segments = [
        [[0, 0, 0], a],
        [[0, 0, 0], b],
        [[0, 0, 0], c],
        [a, a + b],
        [a, a + c],
        [b, b + a],
        [b, b + c],
        [c, c + a],
        [c, c + b],
        [a + b, a + b + c],
        [a + c, a + b + c],
        [b + c, a + b + c],
    ]
    return go_lines(segments, **kwargs)


def plot_fcc_conv():

    fcc_conv = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    fcc_vectors = vectors(
        fcc_conv, name="conv lattice vectors", color="darkblue", width=6
    )
    fcc_box = get_box(fcc_conv, name="conv lattice")

    atoms = go_points(
        [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        size=10,
        color="orange",
        name="atoms",
        legendgroup="atoms",
    )

    import plotly.graph_objects as go
    fig = go.Figure(data=[*fcc_box, *fcc_vectors, atoms])
    return fig


def plot_fcc_prim():
    fcc_prim = np.array([[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]])

    fcc_prim_vectors = vectors(
        fcc_prim, name="prim lattice vectors", color="green", width=6
    )
    fcc_prim_box = get_box(fcc_prim, name="prim lattice", color="green")

    atoms = go_points(
        [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        size=10,
        color="orange",
        name="atoms",
        legendgroup="atoms",
    )

    fcc_conv = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    fcc_conv_box = get_box(fcc_conv, name="conv lattice")

    import plotly.graph_objects as go
    fig = go.Figure(data=[*fcc_prim_box, *fcc_prim_vectors, *fcc_conv_box, atoms])

    return fig


def plot_fcc_100():

    # fcc_100_cell = np.array([[0, 0.5, -0.5], [0, 0.5, 0.5], [1.0, 0.0, 0]])
    fcc_100_cell = np.array([[0.5, -0.5, 0], [0.5, 0.5, 0], [0.0, 0, 1.0]])

    fcc_100_vectors = vectors(
        fcc_100_cell, name="100 lattice vectors", color="red", width=6
    )
    fcc_100_box = get_box(fcc_100_cell, name="100 lattice", color="red")

    fig = plot_fcc_conv()
    fig.add_traces([*fcc_100_box, *fcc_100_vectors])

    return fig


def plot_fcc_110():
    fcc_110_cell = np.array([[0, 0.0, 1.0], [0.5, -0.5, 0], [0.5, 0.5, 0.0]])

    fcc_110_vectors = vectors(
        fcc_110_cell, name="reduced lattice vectors", color="red", width=6
    )
    fcc_110_box = get_box(fcc_110_cell, name="reduced lattice", color="red")

    fig = plot_fcc_conv()
    fig.add_traces([*fcc_110_box, *fcc_110_vectors])
    return fig


def plot_fcc_111():
    fcc_111_cell = np.array([[0.5, 0, -0.5], [0, 0.5, -0.5], [1, 1, 1]])

    fcc_111_vectors = vectors(
        fcc_111_cell, name="reduced lattice vectors", color="red", width=6
    )
    fcc_111_box = get_box(fcc_111_cell, name="reduced lattice", color="red")

    fig = plot_fcc_conv()
    fig.add_traces([*fcc_111_box, *fcc_111_vectors])
    return fig


def plotly_structure(structure, ax=None, to_unit_cell=False, alpha=0.7,
                     style="points+labels", color_scheme="VESTA", **kwargs):
    """
    Plot structure with plotly (minimalistic version).

    Args:
        structure: |Structure| object
        ax: matplotlib :class:`Axes3D` or None if a new figure should be created.
        alpha: The alpha blending value, between 0 (transparent) and 1 (opaque)
        to_unit_cell: True if sites should be wrapped into the first unit cell.
        style: "points+labels" to show atoms sites with labels.
        color_scheme: color scheme for atom types. Allowed values in ("Jmol", "VESTA")

    Returns: |matplotlib-Figure|
    """
    #fig, ax = plot_unit_cell(structure.lattice, ax=ax, linewidth=1)

    box = get_box(structure.lattice.matrix) #, **kwargs):

    from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
    from pymatgen.vis.structure_vtk import EL_COLORS

    #symb2data = {}
    #for symbol in structure.symbol_set:
    #    symb2data[symbol] = d = {}
    #    d["color"] = color = tuple(i / 255 for i in EL_COLORS[color_scheme][symbol])
    #    d["radius"] = CovalentRadius.radius[symbol]
    #    inds = structure.indices_from_symbol(symbol)
    #    sites = [structure[i] for i in inds]
    #    d["xyz"] = []
    #    for site in sites:
    #       if to_unit_cell and hasattr(site, "to_unit_cell"): site = site.to_unit_cell()
    #       Use cartesian coordinates.
    #       x, y, z = site.coords
    #       d["xyz"].append((x, y ,z)

    xyz, sizes, colors = np.empty((len(structure), 3)), [], []
    for i, site in enumerate(structure):
        symbol = site.specie.symbol
        color = tuple(i / 255 for i in EL_COLORS[color_scheme][symbol])
        radius = CovalentRadius.radius[symbol]
        if to_unit_cell and hasattr(site, "to_unit_cell"): site = site.to_unit_cell()
        # Use cartesian coordinates.
        x, y, z = site.coords
        xyz[i] = (x, y, z) # , radius)
        sizes.append(radius)
        colors.append(color)
        #if "labels" in style:
        #    ax.text(x, y, z, symbol)

    atoms = go_points(
        #[[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        xyz,
        size=10,
        color="orange",
        name="atoms",
        legendgroup="atoms",
    )

    #marker = [dict(size=size, color=color) for (size, color) in zip(sizes, colors)]

    #atoms = go.Scatter3d(
    #    x=[v[0] for v in xyz],
    #    y=[v[1] for v in xyz],
    #    z=[v[2] for v in xyz],
    #    #marker=dict(size=size, color=color),
    #    marker=marker,
    #    mode="markers",
    #    #**kwargs
    #)

    # The definition of sizes is not optimal because matplotlib uses points
    # whereas we would like something that depends on the radius (5000 seems to give reasonable plots)
    # For possibile approaches, see
    # https://stackoverflow.com/questions/9081553/python-scatter-plot-size-and-style-of-the-marker/24567352#24567352
    # https://gist.github.com/syrte/592a062c562cd2a98a83
    #if "points" in style:
    #    x, y, z, s = xyzs.T.copy()
    #    s = 5000 * s ** 2
    #    ax.scatter(x, y, zs=z, s=s, c=colors, alpha=alpha)  #facecolors="white", #edgecolors="blue"

    #ax.set_title(structure.composition.formula)
    #ax.set_axis_off()

    #fig = go.Figure(data=[*box, *vectors, atoms])
    import plotly.graph_objects as go
    fig = go.Figure(data=[*box, atoms])
    return fig


# This is the matplotlib API to plot the BZ.

def plotly_wigner_seitz(lattice, fig=None, **kwargs):
    """
    Adds the skeleton of the Wigner-Seitz cell of the lattice to a plotly figure.

    Args:
        lattice: Lattice object
        fig: plotly figure or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to black
            and linewidth to 1.

    Returns: Plotly figure
    """
    #ax, fig, plt = get_ax3d_fig_plt(ax)
    fig, go = get_fig_plotly(fig=fig) #, **fig_kw)

    if "line_color" not in kwargs:
        kwargs["line_color"] = "black"
    if "line_width" not in kwargs:
        kwargs["line_width"] = 1

    bz = lattice.get_wigner_seitz_cell()
    #ax, fig, plt = get_ax3d_fig_plt(ax)

    for iface in range(len(bz)):  # pylint: disable=C0200
        for line in itertools.combinations(bz[iface], 2):
            for jface in range(len(bz)):
                if (iface < jface
                    and any(np.all(line[0] == x) for x in bz[jface])
                    and any(np.all(line[1] == x) for x in bz[jface])):
                    #ax.plot(*zip(line[0], line[1]), **kwargs)
                    fig.add_trace(go_line(line[0], line[1], showlegend=False, **kwargs))

    return fig


def plotly_lattice_vectors(lattice, fig=None, **kwargs):
    """
    Adds the basis vectors of the lattice provided to a plotly figure.

    Args:
        lattice: Lattice object
        fig: plotly figure or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to green
            and linewidth to 3.

    Returns: plotly figure
    """
    fig, go = get_fig_plotly(fig=fig)

    if "line_color" not in kwargs:
        kwargs["line_color"] = "green"
    if "line_width" not in kwargs:
        kwargs["line_width"] = 3
    if "showlegend" not in kwargs:
        kwargs["showlegend"] = False

    vertex1 = lattice.get_cartesian_coords([0.0, 0.0, 0.0])
    vertex2 = lattice.get_cartesian_coords([1.0, 0.0, 0.0])
    fig.add_trace(go_line(vertex1, vertex2, name="a", **kwargs))
    vertex2 = lattice.get_cartesian_coords([0.0, 1.0, 0.0])
    fig.add_trace(go_line(vertex1, vertex2, name="b", **kwargs))
    vertex2 = lattice.get_cartesian_coords([0.0, 0.0, 1.0])
    fig.add_trace(go_line(vertex1, vertex2, name="c", **kwargs))

    return fig


def plotly_path(line, lattice=None, coords_are_cartesian=False, fig=None, **kwargs):
    """
    Adds a line passing through the coordinates listed in 'line' to a plotly figure.

    Args:
        line: list of coordinates.
        lattice: Lattice object used to convert from reciprocal to cartesian coordinates
        coords_are_cartesian: Set to True if you are providing
            coordinates in cartesian coordinates. Defaults to False.
            Requires lattice if False.
        fig: plotly figure or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to red
            and linewidth to 3.

    Returns: plotly figure
    """
    fig, go = get_fig_plotly(fig=fig)

    if "line_color" not in kwargs:
        kwargs["line_color"] = "red"
    if "line_width" not in kwargs:
        kwargs["line_width"] = 3

    for k in range(1, len(line)):
        vertex1 = line[k - 1]
        vertex2 = line[k]
        if not coords_are_cartesian:
            if lattice is None:
                raise ValueError("coords_are_cartesian False requires the lattice")
            vertex1 = lattice.get_cartesian_coords(vertex1)
            vertex2 = lattice.get_cartesian_coords(vertex2)

        fig.add_trace(go_line(vertex1, vertex2, showlegend=False, **kwargs))

    return fig


#def plotly_labels(labels, lattice=None, coords_are_cartesian=False, ax=None, **kwargs):
#    """
#    Adds labels to a matplotlib Axes
#
#    Args:
#        labels: dict containing the label as a key and the coordinates as value.
#        lattice: Lattice object used to convert from reciprocal to cartesian coordinates
#        coords_are_cartesian: Set to True if you are providing.
#            coordinates in cartesian coordinates. Defaults to False.
#            Requires lattice if False.
#        ax: matplotlib :class:`Axes` or None if a new figure should be created.
#        kwargs: kwargs passed to the matplotlib function 'text'. Color defaults to blue
#            and size to 25.
#
#    Returns:
#        matplotlib figure and matplotlib ax
#    """
#    ax, fig, plt = get_ax3d_fig_plt(ax)
#
#    if "color" not in kwargs:
#        kwargs["color"] = "b"
#    if "size" not in kwargs:
#        kwargs["size"] = 25
#
#    for k, coords in labels.items():
#        label = k
#        if k.startswith("\\") or k.find("_") != -1:
#            label = "$" + k + "$"
#        off = 0.01
#        if coords_are_cartesian:
#            coords = np.array(coords)
#        else:
#            if lattice is None:
#                raise ValueError("coords_are_cartesian False requires the lattice")
#            coords = lattice.get_cartesian_coords(coords)
#        ax.text(*(coords + off), s=label, **kwargs)
#
#    return fig, ax


def plotly_points(points, lattice=None, coords_are_cartesian=False, fold=False, labels=None, fig=None, **kwargs):
    """
    Adds points to a plotly figure.

    Args:
        points: list of coordinates
        lattice: Lattice object used to convert from reciprocal to cartesian coordinates
        coords_are_cartesian: Set to True if you are providing
            coordinates in cartesian coordinates. Defaults to False.
            Requires lattice if False.
        fold: whether the points should be folded inside the first Brillouin Zone.
            Defaults to False. Requires lattice if True.
        fig: plotly figure or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'scatter'. Color defaults to blue

    Returns: plotly figure
    """
    fig, go = get_fig_plotly(fig=fig) #, **fig_kw)

    if "marker_color" not in kwargs:
        kwargs["marker_color"] = "blue"

    if (not coords_are_cartesian or fold) and lattice is None:
        raise ValueError("coords_are_cartesian False or fold True require the lattice")

    from pymatgen.electronic_structure.plotter import fold_point
    vecs = []
    for p in points:

        if fold:
            p = fold_point(p, lattice, coords_are_cartesian=coords_are_cartesian)

        elif not coords_are_cartesian:
            p = lattice.get_cartesian_coords(p)

        vecs.append(p)

    kws = dict(textposition="top right", showlegend=False) #, textfont=dict(color='#E58606'))
    kws.update(kwargs)
    fig.add_trace(go_points(vecs, labels=labels, **kws))

    return fig


@add_plotly_fig_kwargs
def plotly_brillouin_zone_from_kpath(kpath, fig=None, **kwargs):
    """
    Gives the plot (as a matplotlib object) of the symmetry line path in
        the Brillouin Zone.

    Args:
        kpath (HighSymmKpath): a HighSymmKPath object
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        **kwargs: provided by add_fig_kwargs decorator

    Returns: plotly figure.
    """
    lines = [[kpath.kpath["kpoints"][k] for k in p] for p in kpath.kpath["path"]]
    return plotly_brillouin_zone(
        bz_lattice=kpath.prim_rec,
        lines=lines,
        fig=fig,
        labels=kpath.kpath["kpoints"],
        show=False,
        **kwargs,
    )


@add_plotly_fig_kwargs
def plotly_brillouin_zone(
    bz_lattice,
    lines=None,
    labels=None,
    kpoints=None,
    fold=False,
    coords_are_cartesian=False,
    fig=None,
    **kwargs,
):
    """
    Plots a 3D representation of the Brillouin zone of the structure.
    Can add to the plot paths, labels and kpoints

    Args:
        bz_lattice: Lattice object of the Brillouin zone
        lines: list of lists of coordinates. Each list represent a different path
        labels: dict containing the label as a key and the coordinates as value.
        kpoints: list of coordinates
        fold: whether the points should be folded inside the first Brillouin Zone.
            Defaults to False. Requires lattice if True.
        coords_are_cartesian: Set to True if you are providing
            coordinates in cartesian coordinates. Defaults to False.
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        kwargs: provided by add_fig_kwargs decorator

    Returns: plotly figure
    """

    fig = plotly_lattice_vectors(bz_lattice, fig=fig)
    plotly_wigner_seitz(bz_lattice, fig=fig)
    if lines is not None:
        for line in lines:
            plotly_path(line, bz_lattice, coords_are_cartesian=coords_are_cartesian, fig=fig)

    if labels is not None:
        # TODO
        #plotly_labels(labels, bz_lattice, coords_are_cartesian=coords_are_cartesian, ax=ax)
        plotly_points(
            labels.values(),
            lattice=bz_lattice,
            coords_are_cartesian=coords_are_cartesian,
            fold=False,
            labels=list(labels.keys()),
            fig=fig,
        )

    if kpoints is not None:
        plotly_points(
            kpoints,
            lattice=bz_lattice,
            coords_are_cartesian=coords_are_cartesian,
            fold=fold,
            fig=fig,
        )

    return fig


def add_colorscale_dropwdowns(fig):
    """
    Add dropdown widgets to change/reverse the colorscale.
    Based on: https://plotly.com/python/dropdowns/#update-several-data-attributes
    """
    button_layer_1_height = 1.30

    # Create list of buttons
    # A single button has the form:
    #
    #    dict(
    #        args=["colorscale", "Viridis"],
    #        label="Viridis",
    #        method="restyle"
    #    ),

    colorscales = ["Viridis", "Cividis", "Blues", "Greens"]

    colorscale_buttons = []
    for cscale in colorscales:
        colorscale_buttons.append(dict(
                args=["colorscale", cscale],
                label=cscale,
                method="restyle",
        ))

    fig.update_layout(
        updatemenus=[
            dict(
                buttons=colorscale_buttons,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.1,
                xanchor="left",
                y=button_layer_1_height,
                yanchor="top"
            ),
            dict(
                buttons=list([
                    dict(
                        args=["reversescale", False],
                        label="False",
                        method="restyle"
                    ),
                    dict(
                        args=["reversescale", True],
                        label="True",
                        method="restyle"
                    )
                ]),
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.37,
                xanchor="left",
                y=button_layer_1_height,
                yanchor="top"
            ),
        ]
    )

    y = button_layer_1_height - 0.02

    fig.update_layout(
        annotations=[
            dict(text="colorscale", x=0, xref="paper", y=y, yref="paper",
                 align="left", showarrow=False),
            dict(text="Reverse<br>Colorscale", x=0.25, xref="paper", y=y,
                 yref="paper", showarrow=False),
    ])

    return fig

def mpl_to_ply(fig, latex=False):
    # Nasty workaround for plotly latex rendering in legend/breaking exception
    def parse_latex(label):
        # Remove latex symobols
        new_label = label.replace("$", "")
        new_label = new_label.replace("\\", "") if not latex else new_label
        new_label = new_label.replace("{", "") if not latex else new_label
        new_label = new_label.replace("}", "") if not latex else new_label
        # plotly latex needs an extra \ for parsing python strings
        # new_label = new_label.replace(" ", "\\ ") if latex else new_label
        # Wrap the label in dollar signs for LaTeX, if needed unless empty``
        new_label = f"${new_label}$" if latex and len(new_label) > 0 else new_label

        return new_label

    for ax in fig.get_axes():
        # TODO improve below logic to add new scatter plots?
        # Loop backwards through the collections to avoid modifying the list as we iterate
        for coll in ax.collections[::-1]:
            if isinstance(coll, mcoll.PathCollection):
                # Use the remove() method to remove the scatter plot collection from the axes
                coll.remove()

        # Process the axis title, x-label, and y-label
        for label in [ax.get_title(), ax.get_xlabel(), ax.get_ylabel()]:
            # Few differences in how mpl and ply parse/encode symbols
            new_label = parse_latex(label)
            # Set the new label
            if label == ax.get_title():
                ax.set_title(new_label)
            elif label == ax.get_xlabel():
                ax.set_xlabel(new_label)
            elif label == ax.get_ylabel():
                ax.set_ylabel(new_label)

        # Check if the axis has a legend
        if ax.get_legend():
            legend = ax.get_legend()
            # Get the legend's text entries
            for text in legend.get_texts():
                label = text.get_text()
                # Remove any existing dollar signs
                new_label = parse_latex(label)
                # Set the new label
                text.set_text(new_label)

    # Convert to plotly figure
    plotly_fig = mpl_to_plotly(fig)

    plotly_fig.update_layout(template  = "plotly_white", title = {
                                "xanchor": "center",
                                "yanchor": "top",
                                "x": 0.5,
                                "font": {
                                    "size": 14
                                },
                            })

    # Iterate over the axes in the figure to retrieve the custom line attributes
    for ax in fig.get_axes():
        if hasattr(ax, '_custom_rc_lines'):
            for rc, color in ax._custom_rc_lines:
                # Add vertical lines to the Plotly figure
                plotly_fig.add_vline(
                    x=rc,
                    line_width=2,
                    line_dash="dash",
                    line_color=color
                )

    # # Loop through each trace and update the hover labels to remove $
    for trace in plotly_fig.data:
        # Retrieve the current label and remove any $ signs
        new_label = trace.name.replace("$", "")

        # Update the trace's name (which is used for the legend label)
        trace.name = new_label

    return plotly_fig
