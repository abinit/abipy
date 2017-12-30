# coding: utf-8
"""Widgets for Pandas Dataframes."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pandas as pd
import ipywidgets as ipw
import abipy.display.utils as ut

from functools import wraps


@wraps(pd.DataFrame.plot)
def plot(data, **kwargs):

    def plot_dataframe(x, y, kind, sharex, sharey, subplots, grid, legend,
                      logx, logy, loglog, colorbar, sort_columns): # pragma: no cover

        x, y = ut.widget2py(x, y)
        sharex, colorbar = ut.str2bool_or_none(sharex, colorbar)
        data.plot(x=x, y=y, kind=kind, subplots=subplots, sharex=None, sharey=sharey,
                layout=None, figsize=None, use_index=True, title=None, grid=grid, legend=legend, style=None,
                logx=logx, logy=logy, loglog=loglog, xticks=None, yticks=None, xlim=None, ylim=None,
                rot=None, fontsize=None, colormap=colorbar, table=False, yerr=None, xerr=None, secondary_y=False,
                sort_columns=sort_columns, **kwargs)
                # There's a typo in the documentation (colorbar/colormap!)
        #import matplotlib.pyplot as plt
        #return plt.gcf()

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                plot_dataframe,
                x=allcols,
                y=allcols,
                sharex=["None", "True", "False"],
                sharey=False,
                kind=["line", "bar", "barh", "hist", "box", "kde", "density", "area", "pie", "scatter", "hexbin"],
                subplots=False,
                grid=True,
                legend=True,
                logx=False,
                logy=False,
                loglog=False,
                colorbar=["None", "True", "False"],
                sort_columns=False,
            )