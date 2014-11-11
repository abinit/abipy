"""Customiced version fo PrettyTable with methods to plot the data with matplotlib."""
from __future__ import division, print_function, absolute_import

import prettytable 


class PrettyTable(prettytable.PrettyTable):
    """Extends ``prettytable.PrettyTable`` adding ``matplotlib`` plots."""

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
        if xname is None and yname is None:
            if len(self.field_names) != 2:
                raise ValueError("xname and yname must be specified if more than 2 columns are present")
            xname, yname = self.field_names

        ix = self.field_names.index(xname)
        xvals = [row[ix] for row in self._rows]
        iy = self.field_names.index(yname)
        yvals = [row[iy] for row in self._rows]

        ax.plot(xvals, yvals, **kwargs)

        if title is not None: fig.suptitle(title)
        if show: plt.show()
        if savefig is not None: fig.savefig(savefig)

        return fig
