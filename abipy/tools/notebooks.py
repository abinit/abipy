"""Tools for ipython notebooks."""
from __future__ import print_function, division, unicode_literals, absolute_import


def mpld3_enable_notebook():
    """Change the default plugins, enable ipython notebook mode and return mpld3 module."""
    import mpld3
    from mpld3 import plugins as plugs
    plugs.DEFAULT_PLUGINS = [plugs.Reset(), plugs.Zoom(), plugs.BoxZoom(), plugs.MousePosition()]
    mpld3.enable_notebook()
    return mpld3


def print_source_in_module(function, module):
    """
    For use inside an IPython notebook: given a module and a function, print the source code.

    Based on:
        http://stackoverflow.com/questions/20665118/how-to-show-source-code-of-a-package-function-in-ipython-notebook
    """
    from inspect import getmembers, isfunction, getsource
    from pygments import highlight
    from pygments.lexers import PythonLexer
    from pygments.formatters import HtmlFormatter
    from IPython.core.display import HTML

    internal_module = __import__(module)

    internal_functions = dict(getmembers(internal_module, isfunction))

    return HTML(highlight(getsource(internal_functions[function]), PythonLexer(), HtmlFormatter(full=True)))


def print_source(function):
    """For use inside an IPython notebook: given a function, print the source code."""
    from inspect import getsource
    from pygments import highlight
    from pygments.lexers import PythonLexer
    from pygments.formatters import HtmlFormatter
    from IPython.core.display import HTML

    return HTML(highlight(getsource(function), PythonLexer(), HtmlFormatter(full=True)))
