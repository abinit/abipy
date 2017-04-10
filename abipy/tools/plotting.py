# coding: utf-8
"""
Utilities for generating matplotlib plots.

.. note:

    Avoid importing matplotlib in the module namespace otherwise startup is very slow.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import collections

from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig_plt


__all__ = [
    "set_axlims",
    "plot_array",
    "ArrayPlotter",
    "data_from_cplx_mode",
    "Marker",
]


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
    if len(lims) == 2:
        left, right = lims[0], lims[1]
    elif len(lims) == 1:
        left = lims[0]
    else:
        left = float(lims)

    set_lim = getattr(ax, {"x": "set_xlim", "y": "set_ylim"}[axname])
    set_lim(left, right)
    return left, right


def data_from_cplx_mode(cplx_mode, arr):
    """
    Extract the data from the numpy array `arr` depending on the values of `cplx_mode`.
    cplx_mode in ("re", "im", "abs", "angle")
    "re" for the real part, "im" for the imaginary part.
    "abs" means that the absolute value of the complex number is shown.
    "angle" will display the phase of the complex number in radians.
    """
    if cplx_mode == "re": return arr.real
    if cplx_mode == "im": return arr.imag
    if cplx_mode == "abs": return np.abs(arr)
    if cplx_mode == "angle": return np.angle(arr, deg=False)
    raise ValueError("Unsupported mode %s" % cplx_mode)


@add_fig_kwargs
def plot_array(array, color_map=None, cplx_mode="abs", **kwargs):
    """
    Use imshow for plotting 2D or 1D arrays.

    Example::

        plot_array(np.random.rand(10,10))

    See http://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python

    Args:
        array: Array-like object (1D or 2D).
        color_map: color map.
        cplx_mode:
            Flag defining how to handle complex arrays. Possible values are in
            ("re", "im", "abs", "angle")
            "re" for the real part, "im" for the imaginary part.
            "abs" means that the absolute value of the complex number is shown.
            "angle" will display the phase of the complex number in radians.

    Returns:
        `matplotlib` figure.
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

    # make a color bar
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
        self._arr_dict = collections.OrderedDict()

        for label, array in labels_and_arrays:
            self.add_array(label, array)

    def __len__(self):
        return len(self._arr_dict)

    def __iter__(self):
        return self._arr_dict.__iter__()

    def items(self):
        return self._arr_dict.items()

    def add_array(self, label, array):
        """Add array with the given name."""
        if label in self._arr_dict:
            raise ValueError("%s is already in %s" % (label, self._arr_dict.keys()))

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
    def plot(self, cplx_mode="abs", color_map="jet", **kwargs):
        """
        Args:
            cplx_mode: "abs" for absolute value, "re", "im", "angle"
            color_map: matplotlib colormap
        """
        #ax, fig, plt = get_ax_fig_plt(None)
        #plt.axis("off")

        # Build grid of plots.
        num_plots, ncols, nrows = len(self), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + (num_plots % ncols)

        import matplotlib.pyplot as plt
        fig, axmat = plt.subplots(nrows=nrows, ncols=ncols, sharex=False, sharey=False, squeeze=False)
        # don't show the last ax if num_plots is odd.
        if num_plots % ncols != 0: axmat[-1, -1].axis("off")

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        from matplotlib.ticker import MultipleLocator

        for ax, (label, arr) in zip(axmat.flat, self.items()):
            data = data_from_cplx_mode(cplx_mode, arr)
            # use origin to place the [0,0] index of the array in the lower left corner of the axes.
            #img = ax.imshow(data, interpolation='nearest', cmap=color_map, origin='lower', aspect="auto")
            img = ax.matshow(data, interpolation='nearest', cmap=color_map, origin='lower', aspect="auto")
            #img = ax.matshow(data, cmap=color_map)

            ax.set_title(label + " (%s)" % cplx_mode)

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


class Marker(collections.namedtuple("Marker", "x y s")):
    """Stores the position and the size of the marker."""

    def __new__(cls, *xys):
        """Extends the base class adding consistency check."""
        assert len(xys) == 3
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

    def extend(self, *xys):
        """Extend the marker values."""
        assert len(xys) == 3

        self.x.extend(xys[0])
        self.y.extend(xys[1])
        self.s.extend(xys[2])

        lens = map(len, self.x, self.y, self.s)
        assert np.all(lens == lens[0])

    def posneg_marker(self):
        """
        Split data into two sets: the first one contains all the points with positive size.
        the first set contains all the points with negative size.
        """
        pos_x, pos_y, pos_s = [], [], []
        neg_x, neg_y, neg_s = [], [], []

        for (x, y, s) in zip(self.x, self.y, self.s):
            if s >= 0.0:
                pos_x.append(x)
                pos_y.append(y)
                pos_s.append(s)
            else:
                neg_x.append(x)
                neg_y.append(y)
                neg_s.append(s)

        return Marker(pos_x, pos_y, pos_s), Marker(neg_x, neg_y, neg_s)


if __name__ == "__main__":
    plot_array(np.random.rand(10, 10))
