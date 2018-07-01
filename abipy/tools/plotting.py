# coding: utf-8
"""
Utilities for generating matplotlib plots.

.. note::

    Avoid importing matplotlib in the module namespace otherwise startup is very slow.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import time
import itertools
import numpy as np

from collections import OrderedDict, namedtuple
from monty.string import list_strings
from monty.functools import lazy_property
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig_plt, get_ax3d_fig_plt, get_axarray_fig_plt


__all__ = [
    "set_axlims",
    "get_ax_fig_plt",
    "get_ax3d_fig_plt",
    "plot_array",
    "ArrayPlotter",
    "data_from_cplx_mode",
    "Marker",
    "plot_unit_cell",
    "GenericDataFilePlotter",
    "GenericDataFilesPlotter",
]


# https://matplotlib.org/gallery/lines_bars_and_markers/linestyles.html
linestyles = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]
)


def ax_append_title(ax, title, loc="center", fontsize=None):
    """Add title to previous ax.title. Return new title."""
    prev_title = ax.get_title(loc=loc)
    new_title = prev_title + title
    ax.set_title(new_title, loc=loc, fontsize=fontsize)
    return new_title


#def set_grid(fig, boolean):
#    if hasattr(fig, "axes"):
#        for ax in fig.axes:
#            if ax.grid: ax.grid.set_visible(boolean)
#    else:
#            if ax.grid: ax.grid.set_visible(boolean)


def set_axlims(ax, lims, axname):
    """
    Set the data limits for the axis ax.

    Args:
        lims: tuple(2) for (left, right), tuple(1) or scalar for left only.
        axname: "x" for x-axis, "y" for y-axis.

    Return: (left, right)
    """
    left, right = None, None
    if lims is None: return (left, right)

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
    set_lim(left, right)

    return left, right


def set_visible(ax, boolean, *args):
    """
    Hide/Show the artists of axis ax listed in args.
    """
    if "legend" in args and ax.legend():
        #handles, labels = ax.get_legend_handles_labels()
        #if handles:
        ax.legend().set_visible(boolean)
    if "title" in args and ax.title:
        ax.title.set_visible(boolean)
    if "xlabel" in args and ax.xaxis.label:
        ax.xaxis.label.set_visible(boolean)
    if "ylabel" in args and ax.yaxis.label:
        ax.yaxis.label.set_visible(boolean)


def rotate_ticklabels(ax, rotation, axname="x"):
    """Rotate the ticklables of axis ``ax``"""
    if "x" in axname:
        for tick in ax.get_xticklabels():
            tick.set_rotation(rotation)
    if "y" in axname:
        for tick in ax.get_yticklabels():
            tick.set_rotation(rotation)


def data_from_cplx_mode(cplx_mode, arr):
    """
    Extract the data from the numpy array ``arr`` depending on the values of ``cplx_mode``.

    Args:
        cplx_mode: Possible values in ("re", "im", "abs", "angle")
            "re" for the real part,
            "im" for the imaginary part.
            "abs" means that the absolute value of the complex number is shown.
            "angle" will display the phase of the complex number in radians.
    """
    if cplx_mode == "re": return arr.real
    if cplx_mode == "im": return arr.imag
    if cplx_mode == "abs": return np.abs(arr)
    if cplx_mode == "angle": return np.angle(arr, deg=False)
    raise ValueError("Unsupported mode `%s`" % str(cplx_mode))


@add_fig_kwargs
def plot_xy_with_hue(data, x, y, hue, decimals=None, ax=None,
                     xlims=None, ylims=None, fontsize=12, **kwargs):
    """
    Plot y = f(x) relation for different values of `hue`.
    Useful for convergence tests done wrt to two parameters.

    Args:
        data: |pandas-DataFrame| containing columns `x`, `y`, and `hue`.
        x: Name of the column used as x-value
        y: Name of the column used as y-value
        hue: Variable that define subsets of the data, which will be drawn on separate lines
        decimals: Number of decimal places to round `hue` columns. Ignore if None
        ax: |matplotlib-Axes| or None if a new figure should be created.
        xlims ylims: Set the data limits for the x(y)-axis. Accept tuple e.g. `(left, right)`
            or scalar e.g. `left`. If left (right) is None, default values are used
        fontsize: Legend fontsize.
        kwargs: Keywork arguments are passed to ax.plot method.

    Returns: |matplotlib-Figure|
    """
    # Check here because pandas error messages are a bit criptic.
    miss = [k for k in (x, y, hue) if k not in data]
    if miss:
        raise ValueError("Cannot find `%s` in dataframe.\nAvailable keys are: %s" % (str(miss), str(data.keys())))

    # Truncate values in hue column so that we can group.
    if decimals is not None:
        data = data.round({hue: decimals})

    ax, fig, plt = get_ax_fig_plt(ax=ax)
    for key, grp in data.groupby(hue):
        # Sort xs and rearrange ys
        xy = np.array(sorted(zip(grp[x], grp[y]), key=lambda t: t[0]))
        xvals, yvals = xy[:, 0], xy[:, 1]

        label = "{} = {}".format(hue, key)
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


@add_fig_kwargs
def plot_array(array, color_map=None, cplx_mode="abs", **kwargs):
    """
    Use imshow for plotting 2D or 1D arrays.

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

    Returns: |matplotlib-Figure|
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


class ArrayPlotter(object):

    def __init__(self, *labels_and_arrays):
        """
        Args:
            labels_and_arrays: List [("label1", arr1), ("label2", arr2")]
        """
        self._arr_dict = OrderedDict()
        for label, array in labels_and_arrays:
            self.add_array(label, array)

    def __len__(self):
        return len(self._arr_dict)

    def __iter__(self):
        return self._arr_dict.__iter__()

    def keys(self):
        return self._arr_dict.keys()

    def items(self):
        return self._arr_dict.items()

    def add_array(self, label, array):
        """Add array with the given name."""
        if label in self._arr_dict:
            raise ValueError("%s is already in %s" % (label, list(self._arr_dict.keys())))

        self._arr_dict[label] = array

    def add_arrays(self, labels, arr_list):
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
    def plot(self, cplx_mode="abs", colormap="jet", fontsize=8, **kwargs):
        """
        Args:
            cplx_mode: "abs" for absolute value, "re", "im", "angle"
            colormap: matplotlib colormap.
            fontsize: legend and label fontsize.

        Returns: |matplotlib-Figure|
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
            return super(cls, Marker).__new__(cls, *xys)

        if len(xys) != 3:
            raise TypeError("Expecting 3 entries in xys got %d" % len(xys))

        x = np.asarray(xys[0])
        y = np.asarray(xys[1])
        s = np.asarray(xys[2])
        xys = (x, y, s)

        for s in xys[-1]:
            if np.iscomplex(s):
                raise ValueError("Found ambiguous complex entry %s" % str(s))

        return super(cls, Marker).__new__(cls, *xys)

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

    def posneg_marker(self):
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

        return self.__class__(pos_x, pos_y, pos_s), Marker(neg_x, neg_y, neg_s)


class MplExpose(object): # pragma: no cover
    """
    Example:

        with MplExpose() as e:
            e(obj.plot1(show=False))
            e(obj.plot2(show=False))
    """
    def __init__(self, slide_mode=False, slide_timeout=None, verbose=1):
        """
        Args:
            slide_mode: If true, iterate over figures. Default: Expose all figures at once.
            slide_timeout: Close figure after slide-timeout seconds Block if None.
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

    def __call__(self, obj):
        """
        Add an object to MplExpose. Support mpl figure, list of figures or
        generator yielding figures.
        """
        import types
        if isinstance(obj, (types.GeneratorType, list, tuple)):
            for fig in obj:
                self.add_fig(fig)
        else:
            self.add_fig(obj)

    def add_fig(self, fig):
        """Add a matplotlib figure."""
        if fig is None: return

        if not self.slide_mode:
            self.figures.append(fig)
        else:
            #print("Printing and closing", fig)
            import matplotlib.pyplot as plt
            if self.timeout_ms is not None:
                # Creating a timer object
                # timer calls plt.close after interval milliseconds to close the window.
                timer = fig.canvas.new_timer(interval=self.timeout_ms)
                timer.add_callback(plt.close, fig)
                timer.start()

            plt.show()
            fig.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement. """
        self.expose()

    def expose(self):
        """Show all figures. Clear figures if needed."""
        if not self.slide_mode:
            print("All figures in memory, elapsed time: %.3f s" % (time.time() - self.start_time))
            import matplotlib.pyplot as plt
            plt.show()
            for fig in self.figures:
                fig.clear()


def plot_unit_cell(lattice, ax=None, **kwargs):
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

    if "color" not in kwargs:
        kwargs["color"] = "k"
    if "linewidth" not in kwargs:
        kwargs["linewidth"] = 3

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

    return fig, ax


def plot_structure(structure, ax=None, to_unit_cell=False, alpha=0.7,
                   style="points+labels", color_scheme="VESTA", **kwargs):
    """
    Plot structure with matplotlib (minimalistic version)

    Args:
        structure: Structure object
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
        if to_unit_cell and hasattr(site, "to_unit_cell"): site = site.to_unit_cell
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
        s = 5000 * s **2
        ax.scatter(x, y, zs=z, s=s, c=colors, alpha=alpha)  #facecolors="white", #edgecolors="blue"

    ax.set_title(structure.composition.formula)
    ax.set_axis_off()

    return fig


def _generic_parser_fh(fh):
    """
    Parse file with data in tabular format. Supports multi datasets a la gnuplot.
    Mainly used for files without any schema, not even CSV

    Args:
        fh: File object

    Returns:
        OrderedDict title --> numpy array
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

    od = OrderedDict()
    for key, data in zip(head_list, arr_list):
        key = " ".join(key.split())
        if key in od:
            print("Header %s already in dictionary. Using new key %s" % (key, 2 * key))
            key = 2 * key
        od[key] = np.array(data).T.copy()

    return od


class GenericDataFilePlotter(object):
    """
    Extract data from a generic text file with results
    in tabular format and plot data with matplotlib.
    Multiple datasets are supported.
    No attempt is made to handle metadata (e.g. column name)
    Mainly used to handle text files written without any schema.
    """
    def __init__(self, filepath):
        with open(filepath, "rt") as fh:
            self.od = _generic_parser_fh(fh)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        lines = []
        for key, arr in self.od.items():
            lines.append("key: `%s` --> array shape: %s" % (key, str(arr.shape)))
        return "\n".join(lines)

    @add_fig_kwargs
    def plot(self, use_index=False, fontsize=8, **kwargs):
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


class GenericDataFilesPlotter(object):

    @classmethod
    def from_files(cls, filepaths):
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

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        lines = []
        app = lines.append
        for od, filepath in zip(self.odlist, self.filepaths):
            app("File: %s" % filepath)
            for key, arr in od.items():
                lines.append("\tkey: `%s` --> array shape: %s" % (key, str(arr.shape)))

        return "\n".join(lines)

    def add_file(self, filepath):
        """Add data from `filepath`"""
        with open(filepath, "rt") as fh:
            self.odlist.append(_generic_parser_fh(fh))
            self.filepaths.append(filepath)

    @add_fig_kwargs
    def plot(self, use_index=False, fontsize=8, colormap="viridis", **kwargs):
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
