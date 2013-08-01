"""
Utilities for generating matplotlib plots.

.. note:
    Avoid importing these tools in the top-level namespace of the module
    since they rely on matplotlib whose startup is very slow.
"""
from __future__ import division, print_function

import numpy as np
import collections
import matplotlib as mpl
import matplotlib.cm as cmap

from matplotlib import pyplot as plt

__all__ = [
    "plot_array",
    "ArrayPlotter"
]


def plot_array(array, color_map=None, cplx_mode="abs", *args, **kwargs):
    """
    Use imshow for plotting 2D or 1D arrays.

    Example::
        plot_array(np.random.rand(10,10))

    See http://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python

    Args:
        array:
            Array-like object (1D or 2D).
        color_map:
            color map.
        cplx_mode:
            Flag defining how to handle complex arrays. Possible values are "abs", "angle"
            "abs" means that the absolute value of the complex number is shown.
            "angle" will display the phase of the complex number in radians.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

    Returns:
        `matplotlib` figure.
    """
    title = kwargs.pop("title", None)
    show = kwargs.pop("show", True)
    savefig = kwargs.pop("savefig", None)

    # Handle vectors
    array = np.atleast_2d(array)

    # Handle complex arrays.
    if np.iscomplexobj(array):
        if cplx_mode == "abs":
            array = np.abs(array)
        elif cplx_mode == "angle":
            array = np.angle(array, deg=False)
        else:
            raise ValueError("Unsupported mode %s" % cplx_mode)

    if color_map is None:
        # make a color map of fixed colors
        color_map = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                                 ['blue', 'black', 'red'], 256)

    img = plt.imshow(array, interpolation='nearest', cmap=color_map, origin='lower')

    # make a color bar
    plt.colorbar(img, cmap=color_map)

    # Set grid
    plt.grid(True, color='white')

    if title:
        plt.title(title)

    if show:
        plt.show()

    fig = plt.gcf()
    if savefig is not None:
        fig.savefig(savefig)

    return fig


class ArrayPlotter(object):

    def __init__(self):
        self._arr_dict = collections.OrderedDict()

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
            labels:
                List of labels.
            arr_list:
                List of arrays.
        """
        assert len(labels) == len(arr_list)
        for (label, arr) in zip(labels, arr_list):
            self.add_array(label, arr)

    def plot(self, color_map="jet", **kwargs):
        """
        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        # Grid parameters
        num_plots, ncols, nrows = len(self), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        import matplotlib.pyplot as plt

        fig = plt.figure()

        for n, (label, arr) in enumerate(self.items()):
            fig.add_subplot(nrows, ncols, n)
            img = plt.imshow(np.abs(arr), interpolation='nearest', cmap=color_map, origin='lower')
            # make a color bar
            plt.colorbar(img, cmap=color_map)
            plt.title(label)
            # Set grid
            plt.grid(True, color='white')

        if title:
            fig.suptitle(title)

        #plt.tight_layout()

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig


class Marker(collections.namedtuple("Marker", "x y s")):
    """Stores the position and the size of the marker."""

    def __new__(cls, *xys):
        """Extends the base class adding consistency check."""
        assert len(xys) == 3

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
