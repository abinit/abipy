# coding: utf-8
r"""
Widgets for Pandas Dataframes based on seaborn API.

API reference

Distribution plots
jointplot(x, y[, data, kind, stat_func, ...])	Draw a plot of two variables with bivariate and univariate graphs.
pairplot(data[, hue, hue_order, palette, ...])	Plot pairwise relationships in a dataset.
distplot(a[, bins, hist, kde, rug, fit, ...])	Flexibly plot a univariate distribution of observations.
kdeplot(data[, data2, shade, vertical, ...])	Fit and plot a univariate or bivariate kernel density estimate.
rugplot(a[, height, axis, ax])	Plot datapoints in an array as sticks on an axis.

Regression plots
lmplot(x, y, data[, hue, col, row, palette, ...])	Plot data and regression model fits across a FacetGrid.
regplot(x, y[, data, x_estimator, x_bins, ...])	Plot data and a linear regression model fit.
residplot(x, y[, data, lowess, x_partial, ...])	Plot the residuals of a linear regression.
coefplot(formula, data[, groupby, ...])	Plot the coefficients from a linear model.

# Categorical plots
factorplot([x, y, hue, data, row, col, ...])	Draw a categorical plot onto a FacetGrid.
boxplot([x, y, hue, data, order, hue_order, ...])	Draw a box plot to show distributions with respect to categories.
violinplot([x, y, hue, data, order, ...])	Draw a combination of boxplot and kernel density estimate.
stripplot([x, y, hue, data, order, ...])	Draw a scatterplot where one variable is categorical.
swarmplot([x, y, hue, data, order, ...])	Draw a categorical scatterplot with non-overlapping points.
pointplot([x, y, hue, data, order, ...])	Show point estimates and confidence intervals using scatter plot glyphs.
barplot([x, y, hue, data, order, hue_order, ...])	Show point estimates and confidence intervals as rectangular bars.
countplot([x, y, hue, data, order, ...])	Show the counts of observations in each categorical bin using bars.

Matrix plots
heatmap(data[, vmin, vmax, cmap, center, ...])	Plot rectangular data as a color-encoded matrix.
clustermap(data[, pivot_kws, method, ...])	Plot a hierarchically clustered heatmap of a pandas DataFrame

Timeseries plots
tsplot(data[, time, unit, condition, value, ...])	Plot one or more timeseries with flexible representation of uncertainty.

Miscellaneous plots
palplot(pal[, size])	Plot the values in a color palette as a horizontal array.

Axis grids
FacetGrid(data[, row, col, hue, col_wrap, ...])	Subplot grid for plotting conditional relationships.
PairGrid(data[, hue, hue_order, palette, ...])	Subplot grid for plotting pairwise relationships in a dataset.
JointGrid(x, y[, data, size, ratio, space, ...])	Grid for drawing a bivariate plot with marginal univariate plots.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import ipywidgets as ipw
import seaborn as sns
import abipy.display.utils as ut

from functools import wraps
from collections import OrderedDict
from IPython.display import display, clear_output

__all__ = [
    #"api_selector",
    # Distribution plots
    "jointplot",
    "pairplot",
    #"distplot",
    #"kdeplot",
    #"rugplot",
    # Regression plots
    "lmplot",
    #"regplot",
    #"residplot",
    #"coefplot",
    # Categorical plots
    "factorplot",
    "boxplot",
    "violinplot",
    "stripplot",
    "swarmplot",
    "pointplot",
    "barplot",
    "countplot",
    # Matrix plots
    #"heatmap",
    #"clustermap",
    # Timeseries plots
    #"tsplot",
    # Miscellaneous plots
    #"palplot",
]


#def api_selector(data, funcname="countplot"):
#    """
#    A widgets with ToogleButtons that allow the user to select and display
#    the widget associated to the different seaborn functions.
#    """
#    this_module = sys.modules[__name__]
#    name2wfunc = OrderedDict()
#    for a in __all__:
#        if a == "api_selector": continue
#        func = this_module.__dict__.get(a)
#        if not callable(func): continue
#        name2wfunc[func.__name__] = func
#
#    w1 = ipw.ToggleButtons(description='seaborn API', options=list(name2wfunc.keys()))
#    w1.value = funcname
#    w2 = name2wfunc[funcname](data)
#    box = ipw.VBox(children=[w1, w2])
#
#    def on_value_change(change):
#        #print(change)
#        box.close()
#        clear_output()
#        api_selector(data, funcname=change["new"])
#    w1.observe(on_value_change, names='value')
#
#    return display(box)

#def _help_button(funcname):
#    btn = ipw.Button(description="Help", button_style='info', tooltip='Click me', icon='check')
#    btn.value = "help"
#    def func(btn):
#        open_seaborn_doc(funcname)
#    btn.on_click(func)
#    return btn
#
#
#def open_seaborn_doc(funcname):
#    # http://seaborn.pydata.org/generated/seaborn.pairplot.html#seaborn.pairplot"
#    url = "http://seaborn.pydata.org/generated/seaborn.%s.html#seaborn.%s" % (funcname, funcname)
#    import webbrowser
#    try:
#        return webbrowser.open(url)
#    except webbrowser.Error as exc:
#        # Warn the user and ignore the exception.
#        import warnings
#        warnings.warn(str(exc) + "\nwhile trying to open url: %s" % url)


@wraps(sns.jointplot)
def jointplot(data, joint_kws=None, marginal_kws=None, annot_kws=None, **kwargs):

    def sns_joinplot(x, y, kind, color):  # pragma: no cover
        x, y, color = ut.widget2py(x, y, color)
        # TODO: stat_func
        return sns.jointplot(x, y, data=data, kind=kind, # stat_func=<function pearsonr>,
                            color=color, size=6, ratio=5, space=0.2, dropna=True, xlim=None, ylim=None,
                            joint_kws=joint_kws, marginal_kws=marginal_kws, annot_kws=annot_kws, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_joinplot,
                x=allcols,
                y=allcols,
                kind=["scatter", "reg", "resid", "kde", "hex"],
                color=ut.colors_dropdow(),
                #button=_help_button("joinplot"),
            )


@wraps(sns.pairplot)
def pairplot(data, plot_kws=None, diag_kws=None, grid_kws=None):
    # TODO: Write widget with multiple checkboxes to implement lists.

    def sns_pairplot(x_vars, y_vars, hue, kind, diag_kind): # pragma: no cover
        x_vars, y_vars, hue = ut.widget2py(x_vars, y_vars, hue)
        return sns.pairplot(data, hue=hue, hue_order=None, palette=None, vars=None, x_vars=x_vars, y_vars=y_vars,
                     kind=kind, diag_kind=diag_kind, markers=None, size=2.5, aspect=1, dropna=True,
                     plot_kws=plot_kws, diag_kws=diag_kws, grid_kws=grid_kws)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_pairplot,
                x_vars=allcols,
                y_vars=allcols,
                hue=allcols,
                kind=["scatter", "ref"],
                diag_kind=["hist", "kde"],
            )


"""
@wraps(sns.distplot)
def distplot(data, fit=None, hist_kws=None, kde_kws=None, rug_kws=None, fit_kws=None):

    def sns_distplot(hist, kde, rug, color, vertical, norm_hist):
        color = ut.widget2py(color)
        ax, fig, _ = ut.get_ax_fig_plt()
        return sns.distplot(a, bins=None, hist=hist, kde=kde, rug=rug, fit=fit,
                     hist_kws=hist_kws, kde_kws=kde_kws, rug_kws=rug_kws, fit_kws=fit_kws,
                     color=clor, vertical=vertical, norm_hist=norm_hist, axlabel=None, label=None, ax=ax)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_distplot,
                hist=True,
                kde=True,
                rug=False,
                color=ut.colors_dropdow(),
                vertical=False,
                norm_hist=False,
            )


@wraps(sns.kdeplot)
def kdeplot(data, **kwargs):

    def sns_kdeplot()
        color = ut.widget2py(color)
        ax, fig, _ = ut.get_ax_fig_plt()

        return sns.kdeplot(data, data2=None, shade=False, vertical=False, kernel='gau', bw='scott',
                    gridsize=100, cut=3, clip=None, legend=True, cumulative=False, shade_lowest=True, ax=ax, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_kdeplot,
                color=ut.colors_dropdow(),
            )
"""

####################
# Regression plots #
####################

@wraps(sns.lmplot)
def lmplot(data, scatter_kws=None, line_kws=None):

    def sns_lmplot(x, y, hue, col, row, legend, size):  # pragma: no cover
        x, y, hue, col, row = ut.widget2py(x, y, hue, col, row)

        return sns.lmplot(x, y, data, hue=hue, col=col, row=row, palette=None, col_wrap=None,
                   size=size, aspect=1, markers='o', sharex=True, sharey=True, hue_order=None,
                   col_order=None, row_order=None, legend=legend, legend_out=True,
                   x_estimator=None, x_bins=None, x_ci='ci', scatter=True, fit_reg=True,
                   ci=95, n_boot=1000, units=None, order=1, logistic=False, lowess=False, robust=False,
                   logx=False, x_partial=None, y_partial=None, truncate=False, x_jitter=None, y_jitter=None,
                   scatter_kws=scatter_kws, line_kws=line_kws)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_lmplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                col=allcols,
                row=allcols,
                legend=True,
                size=ut.size_slider(default=5),
            )


#@wraps(sns.interactplot)
#def interactplot(data, contour_kws=None, scatter_kws=None, **kwargs):
#
#    def sns_interactplot(x1, x2, y, filled, colorbar, logistic):
#        ax, fig, _ = ut.get_ax_fig_plt()
#        return sns.interactplot(x1, x2, y, data=data, filled=filled, cmap='RdBu_r', colorbar=colorbar,
#                         levels=30, logistic=logistic, contour_kws=contour_kws, scatter_kws=scatter_kws,
#                         ax=ax, **kwargs)
#
#    allcols = list(data.keys())
#    return ipw.interact_manual(
#                sns_interactplot,
#                x1=allcols,
#                x2=allcols,
#                y=allcols,
#                filled=False,
#                colorbar=True,
#                logistic=False,
#            )


#####################
# Categorical plots #
#####################

@wraps(sns.factorplot)
def factorplot(data, facet_kws=None, **kwargs):

    def sns_factorplot(x, y, hue, color, kind, size, legend):  # pragma: no cover
        x, y, hue, color = ut.widget2py(x, y, hue, color)
        return sns.factorplot(x=x, y=y, hue=hue, data=data, row=None, col=None, col_wrap=None, # estimator=<function mean>,
                       ci=95, n_boot=1000, units=None, order=None, hue_order=None, row_order=None, col_order=None,
                       kind=kind, size=size, aspect=1, orient=None, color=color, palette=None,
                       legend=legend, legend_out=True, sharex=True, sharey=True, margin_titles=False,
                       facet_kws=facet_kws, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_factorplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                color=ut.colors_dropdow(),
                kind=["point", "bar", "count", "box", "violin", "strip"],
                size=ut.size_slider(default=4),
                legend=True,
            )


@wraps(sns.boxplot)
def boxplot(data, **kwargs):

    def sns_boxplot(x, y, hue, orient, color, saturation, notch): # pragma: no cover
        x, y, hue, orient, color = ut.widget2py(x, y, hue, orient, color)
        ax, fig, _ = ut.get_ax_fig_plt()
        return sns.boxplot(x=x, y=y, hue=hue, data=data, order=None, hue_order=None, orient=orient,
                          color=color, palette=None, saturation=saturation, width=0.8, fliersize=5, linewidth=None,
                          whis=1.5, notch=notch, ax=ax, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_boxplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                orient=["None", "v", "h"],
                color=ut.colors_dropdow(),
                saturation=ut.saturation_slider(default=0.75),
                notch=False,
            )


@wraps(sns.violinplot)
def violinplot(data, **kwargs):

    def sns_violinplot(x, y, hue, bw, scale, inner, split, orient, color, saturation): # pragma: no cover
        x, y, hue, inner, orient, color = ut.widget2py(x, y, hue, inner, orient, color)
        ax, fig, _ = ut.get_ax_fig_plt()

        sns.violinplot(x=x, y=y, hue=hue, data=data, order=None, hue_order=None,
                       bw=bw, cut=2, scale=scale, scale_hue=True,
                       gridsize=100, width=0.8, inner=inner, split=split, orient=orient,
                       linewidth=None, color=color, palette=None, saturation=saturation, ax=ax, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_violinplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                bw=["scott", "silverman", "float"],
                scale=["area", "count", "width"],
                inner=["box", "quartile", "point", "stick", "None"],
                split=False,
                orient=["None", "v", "h"],
                color=ut.colors_dropdow(),
                saturation=ut.saturation_slider(default=0.75),
            )


@wraps(sns.stripplot)
def stripplot(data, **kwargs):

    def sns_stripplot(x, y, hue, split, orient, color, size, linewidth): # pragma: no cover
        x, y, hue, orient, color = ut.widget2py(x, y, hue, orient, color)
        ax, fig, _ = ut.get_ax_fig_plt()
        return sns.stripplot(x=x, y=y, hue=hue, data=data, order=None, hue_order=None, jitter=False,
                            split=split, orient=orient, color=color, palette=None, size=size, edgecolor='gray',
                            linewidth=linewidth, ax=ax, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_stripplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                split=False,
                orient=["None", "v", "h"],
                color=ut.colors_dropdow(),
                size=ut.size_slider(default=5),
                linewidth=ut.linewidth_slider(default=0),
            )


@wraps(sns.swarmplot)
def swarmplot(data, **kwargs):

    def sns_swarmplot(x, y, hue, split, orient, color, size, linewidth): # pragma: no cover
        x, y, hue, orient, color = ut.widget2py(x, y, hue, orient, color)
        ax, fig, _ = ut.get_ax_fig_plt()
        return sns.swarmplot(x=x, y=y, hue=hue, data=data, order=None, hue_order=None,
                            split=split, orient=orient, color=color, palette=None, size=size,
                            edgecolor='gray', linewidth=linewidth, ax=ax, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_swarmplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                split=False,
                orient=["None", "v", "h"],
                color=ut.colors_dropdow(),
                size=ut.size_slider(default=5),
                linewidth=ut.linewidth_slider(default=0),
            )


@wraps(sns.pointplot)
def pointplot(data, **kwargs):

    def sns_pointplot(x, y, hue, split, join, orient, color, linewidth): # pragma: no cover
        x, y, hue, orient, color = ut.widget2py(x, y, hue, orient, color)
        ax, fig, _ = ut.get_ax_fig_plt()
        return sns.pointplot(x=x, y=y, hue=hue, data=data, order=None, hue_order=None, # estimator=<function mean>,
                            ci=95, n_boot=1000, units=None, markers='o', linestyles='-', dodge=False, join=join, scale=1,
                            orient=orient, color=color, palette=None, ax=ax, errwidth=None, capsize=None, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_pointplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                split=False,
                join=True,
                orient=["None", "v", "h"],
                color=ut.colors_dropdow(),
                linewidth=ut.linewidth_slider(default=0),
            )


@wraps(sns.barplot)
def barplot(data, **kwargs):

    def sns_barplot(x, y, hue, orient, color, saturation): # pragma: no cover
        x, y, hue, orient, color = ut.widget2py(x, y, hue, orient, color)
        ax, fig, _ = ut.get_ax_fig_plt()
        return sns.barplot(x=x, y=y, hue=hue, data=data, order=None, hue_order=None, # estimator=<function mean>,
                           ci=95, n_boot=1000, units=None, orient=orient, color=color, palette=None,
                           saturation=saturation, errcolor='.26', ax=ax, **kwargs) # errwidth=None, capsize=None, # New args added in ??

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_barplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                orient=["None", "v", "h"],
                color=ut.colors_dropdow(),
                saturation=ut.saturation_slider(default=0.75),
            )


@wraps(sns.countplot)
def countplot(data, **kwargs):

    def sns_countplot(x, y, hue, color, saturation): # pragma: no cover
        x, y, hue, color = ut.widget2py(x, y, hue, color)
        ax, fig, _ = ut.get_ax_fig_plt()
        return sns.countplot(x=x, y=y, hue=hue, data=data, order=None, hue_order=None, orient=None,
                             color=color, palette=None, saturation=saturation, ax=ax, **kwargs)

    allcols = ["None"] + list(data.keys())
    return ipw.interact_manual(
                sns_countplot,
                x=allcols,
                y=allcols,
                hue=allcols,
                color=ut.colors_dropdow(),
                saturation=ut.saturation_slider(default=0.75),
            )

################
# Matrix plots #
################

@wraps(sns.heatmap)
def heatmap(data, annot_kws=None, cbar_kws=None, **kwargs):

    def sns_heatmap(): # pragma: no cover
        ax, fig, _ = ut.get_ax_fig_plt()
        return sns.heatmap(data, vmin=None, vmax=None, cmap=None, center=None, robust=False, annot=None,
                           fmt='.2g', annot_kws=annot_kws, linewidths=0, linecolor='white', cbar=True,
                           cbar_kws=cbar_kws, cbar_ax=None, square=False, ax=ax,
                           xticklabels=True, yticklabels=True, mask=None, **kwargs)

    return ipw.interact_manual(
                sns_heatmap,
            )


@wraps(sns.clustermap)
def clustermap(data, pivot_kws=None, cbar_kws=None, **kwargs):

    def sns_clustermap(): # pragma: no cover
        return sns.clustermap(data, pivot_kws=pivot_kws, method='average', metric='euclidean',
                              z_score=None, standard_scale=None, figsize=None, cbar_kws=cbar_kws,
                              row_cluster=True, col_cluster=True, row_linkage=None, col_linkage=None,
                              row_colors=None, col_colors=None, mask=None, **kwargs)

    return ipw.interact_manual(
                sns_clustermap,
            )
