# coding: utf-8
"""Tools to build ipython widgets."""
from __future__ import print_function, division, unicode_literals, absolute_import

import ipywidgets as ipw

from collections import OrderedDict
from abipy.tools.plotting import get_ax_fig_plt


#def add_docstrings(*tuples):
#    """
#    This decorator adds to the docstring the documentation for functions.
#    When writing high-level API, it's quite common to call thirdy-party functions
#    with a restricted set of arguments while optional keyword arguments are
#    collected in an optional dictionary.
#
#    The first item of the tuple contains the function (python object) wrapped by the code.
#    The second item is list of strings with the name of the actual arguments passed to function.
#    """
#    from functools import wraps
#    def wrapper(func):
#        @wraps(func)
#        def wrapped_func(*args, **kwargs):
#            return func(*args, **kwargs)
#
#        # Add docstrings for the functions that will be called by func.
#        lines = []
#        app = lines.append
#        for t in tuples:
#            fname = t[0].__name__
#            # List of strings or string.
#            if isinstance(t[1], (list, tuple)):
#                fargs = ",".join("`%s`" % a for a in t[1])
#            else:
#                fargs = "`%s`" % t[1]
#            app("\n%s are passed to function :func:`%s` in module :mod:`%s`" % (fargs, fname, t[0].__module__))
#            app("Docstring of `%s`:" % fname)
#            app(t[0].__doc__)
#        s = "\n".join(lines)
#
#        if wrapped_func.__doc__ is not None:
#            # Add s at the end of the docstring.
#            wrapped_func.__doc__ += "\n" + s
#        else:
#            # Use s
#            wrapped_func.__doc__ = s
#        return wrapped_func
#    return wrapper


def widget2py(*args):
    l = [None if a == "None" else a for a in args]
    return l[0] if len(l) == 1 else l


def str2bool_or_none(*args):
    d = {"None": None, "True": True, "False": False}
    l = [d[a] for a in args]
    return l[0] if len(l) == 1 else l


# Taken from matplotlib.markers.MarkerStyle (replaced dict with OrderedDict).
_mpl_markers = OrderedDict([
    ('.', 'point'),
    (',', 'pixel'),
    ('o', 'circle'),
    ('v', 'triangle_down'),
    ('^', 'triangle_up'),
    ('<', 'triangle_left'),
    ('>', 'triangle_right'),
    ('1', 'tri_down'),
    ('2', 'tri_up'),
    ('3', 'tri_left'),
    ('4', 'tri_right'),
    ('8', 'octagon'),
    ('s', 'square'),
    ('p', 'pentagon'),
    ('*', 'star'),
    ('h', 'hexagon1'),
    ('H', 'hexagon2'),
    ('+', 'plus'),
    ('x', 'x'),
    ('D', 'diamond'),
    ('d', 'thin_diamond'),
    ('|', 'vline'),
    ('_', 'hline'),
    #(TICKLEFT: 'tickleft',
    #(TICKRIGHT: 'tickright',
    #(TICKUP: 'tickup',
    #(TICKDOWN: 'tickdown',
    #(CARETLEFT: 'caretleft',
    #(CARETRIGHT: 'caretright',
    #(CARETUP: 'caretup',
    #(CARETDOWN: 'caretdown',
    ("None", 'nothing'),
    (None, 'nothing'),
    (' ', 'nothing'),
    ('', 'nothing'),
])


def markers_dropdown(default="o"):
    return ipw.Dropdown(
        options={name: key for key, name in _mpl_markers.items()},
        value=default,
        description='marker',
    )


_mpl_colors = OrderedDict([
    ("None", "None"),
    ("blue", "b"),
    ("green", "g"),
    ("red", "r"),
    ("cyan", "c"),
    ("magenta", "m"),
    ("yellow", "y"),
    ("black", "k"),
    ("white", "w"),
])


def colors_dropdow(default="None"):
    return ipw.Dropdown(
        options=_mpl_colors,
        value=default,
        description='color',
    )


def linewidth_slider(default=1, orientation="horizontal"):
    return ipw.FloatSlider(
        value=default,
        min=0,
        max=10,
        step=0.5,
        description='linewidth',
        orientation=orientation,
        readout_format='.1f'
    )


def size_slider(default=5, orientation="horizontal"):
    return ipw.FloatSlider(
        value=default,
        min=0,
        max=20,
        step=0.5,
        description='size',
        orientation=orientation,
        readout_format='.1f'
    )


def saturation_slider(default=0.75, orientation="horizontal"):
    return ipw.FloatSlider(
        value=default,
        min=0,
        max=1,
        step=0.05,
        description='saturation',
        orientation=orientation,
        readout_format='.1f'
    )

# Have colormaps separated into categories:
# http://matplotlib.org/examples/color/colormaps_reference.html
_mpl_categ_cmaps = OrderedDict([
        #('Perceptually Uniform Sequential',
    ('Uniform',        ['viridis', 'inferno', 'plasma', 'magma']),
    ('Sequential',     ['Blues', 'BuGn', 'BuPu',
                        'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
                        'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
                        'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']),
    ('Sequential(2)',  ['afmhot', 'autumn', 'bone', 'cool',
                        'copper', 'gist_heat', 'gray', 'hot',
                        'pink', 'spring', 'summer', 'winter']),
    ('Diverging',      ['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                        'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                        'seismic']),
    ('Qualitative',    ['Accent', 'Dark2', 'Paired', 'Pastel1',
                        'Pastel2', 'Set1', 'Set2', 'Set3']),
    ('Miscellaneous',  ['gist_earth', 'terrain', 'ocean', 'gist_stern',
                        'brg', 'CMRmap', 'cubehelix',
                        'gnuplot', 'gnuplot2', 'gist_ncar',
                        'nipy_spectral', 'jet', 'rainbow',
                        'gist_rainbow', 'hsv', 'flag', 'prism'])
])

# flat list.
_mpl_cmaps = [cm for sublist in _mpl_categ_cmaps.values() for cm in sublist]

def colormap_widget(default=None):
    options = _mpl_cmaps
    value = options[0]
    if default is not None:
        value = default
        if default not in _mpl_cmaps: options[:].insert(0, value)
    return ipw.Dropdown(options=options, value=value, description='colormap')


#def colormap_widget():
#    from IPython.display import display, clear_output
#    w_type = ipw.Dropdown(options=list(_mpl_categ_cmaps.keys()), description='colormap category')
#    w_cmap = ipw.Dropdown(options=_mpl_categ_cmaps["Uniform"], description='colormap name')
#
#    def on_value_change(change):
#        print(change)
#        print(w_cmap.value)
#        w_cmap.options = _mpl_categ_cmaps[w_type.value]
#        print(w_cmap.value)
#
#    w_type.observe(on_value_change, names='value')
#    box = ipw.HBox(children=[w_type, w_cmap])
#    return display(box)
