# coding: utf-8
"""Customiced version fo PrettyTable with methods to plot the data with matplotlib."""
from __future__ import division, print_function, absolute_import

from collections import namedtuple

import prettytable 
import numpy as np


class PrettyTable(prettytable.PrettyTable):
    """
    Extends ``prettytable.PrettyTable`` adding methods for numerical analysis and 
    methods to produce ``matplotlib`` plots.
    """
    def quadfit(self, xname=None, yname=None):
        """
        Quadratic fit. Return namedtuple with the parameters of the fix a*x**2 + b*x + c,
        the position of the minimum in x0 and the value of the minimum in y0
        """
        xvals, yvals = self.get_float_cols(xname, yname)

        a, b, c = np.polyfit(xvals, yvals, 2)
        x0 = -b/(2*a)
        y0 = a*x0**2 + b*x0 + c

        return namedtuple("quadfit_results", "a, b, c, x0, y0")(a=a, b=b, c=c, x0=x0, y0=y0)

    def plot(self, xname=None, yname=None, **kwargs):
        """
        Plot yname as function of xname. xname and yname can be omitted if only two columns are present.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            `matplotlib` figure
        """
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        # Extract data from the table.
        xvals, yvals = self.get_float_cols(xname, yname)

        ax.plot(xvals, yvals, **kwargs)

        if title is not None: fig.suptitle(title)
        if show: plt.show()
        if savefig is not None: fig.savefig(savefig)

        return fig

    def get_float_cols(self, xname, yname):
        """
        Extract values in columns xname, yname. Return numpy array of floats.
        if xname and yname are None, the first two columns are returned.
        """
        if xname is None and yname is None:
            if len(self.field_names) != 2:
                raise ValueError("xname and yname must be specified if more than 2 columns are present")
            xname, yname = self.field_names
                                                                                                         
        ix = self.field_names.index(xname)
        xvals = np.array([float(row[ix]) for row in self._rows])
        iy = self.field_names.index(yname)
        yvals = np.array([float(row[iy]) for row in self._rows])

        return xvals, yvals


