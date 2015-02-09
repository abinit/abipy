"""Tools for ipython notebooks."""
from __future__ import print_function, division, unicode_literals


def mpld3_enable_notebook():
    """Change the default plugins, enable ipython notebook mode and return mpld3 module."""
    import mpld3
    from mpld3 import plugins as plugs
    plugs.DEFAULT_PLUGINS = [plugs.Reset(), plugs.Zoom(), plugs.BoxZoom(), plugs.MousePosition()]
    mpld3.enable_notebook()
    return mpld3
