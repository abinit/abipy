"""Tools for ipython notebooks."""
from __future__ import print_function, division, unicode_literals, absolute_import


def print_source_in_module(function, module):  # pragma: no cover
    """
    For use inside an jupyter_ notebook: given a module and a function, print the source code.

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


def print_source(function, **kwargs):  # pragma: no cover
    """
    For use inside a jupyter_ notebook: given a function, print the source code.

    Args:
        **kwargs: Passed to HtmlFormatter

    Return:
        HTML string.
    """
    from inspect import getsource
    from pygments import highlight
    from pygments.lexers import PythonLexer
    from pygments.formatters import HtmlFormatter
    from IPython.core.display import HTML

    if "full" not in kwargs: kwargs["full"] = True
    return HTML(highlight(getsource(function), PythonLexer(), HtmlFormatter(**kwargs)))


def print_doc(function, **kwargs):  # pragma: no cover
    """
    For use inside a jupyter_ notebook: given a function, print the docstring.

    Args:
        **kwargs: Passed to HtmlFormatter

    Return:
        HTML string.
    """
    from inspect import getsource
    from pygments import highlight
    from pygments.lexers import PythonLexer
    from pygments.formatters import HtmlFormatter
    from IPython.core.display import HTML

    # Extract source code up to end of docstring.
    lines, count = [], 0
    for l in getsource(function).splitlines():
        lines.append(l)
        if l.lstrip().startswith('"""'): count += 1
        if count == 2: break

    if "full" not in kwargs: kwargs["full"] = True
    return HTML(highlight("\n".join(lines), PythonLexer(), HtmlFormatter(**kwargs)))


def ipw_listdir(top=".", recurse=True, widget_type="dropdown"):   # pragma: no cover
    """
    Return an ipython widget listing all the files located within the directory ``top``
    that can be inspected with :ref:`abiopen.py`. The user can select the file in the widget
    and print info on the corresponding file inside the notebook.

    Args:
        top: Initial directory.
        recurse: False to ignore directories within ``top``.
        widget_type: Specify the widget to create. Possible values in:
            ["tooglebuttons", "dropdown", "radiobuttons"]
    """
    from abipy import abilab
    from IPython.display import display, clear_output
    import ipywidgets as ipw

    # Select the widget class from widget_type
    d = dict(
        tooglebuttons=ipw.ToggleButtons,
        dropdown=ipw.Dropdown,
        radiobuttons=ipw.RadioButtons,
    )
    try:
        widget_class = d[widget_type]
    except KeyError:
        raise KeyError("Invalid `widget_type`: %s, Choose among: %s" % (widget_type, str(list(d.keys()))))

    def on_value_change(change):
        """Callback"""
        clear_output()
        path = change["new"]
        #print(change)
        with abilab.abiopen(path) as abifile:
            print(abifile)
            #display(abifile)

    # Get dict: dirname --> list_of_files supported by abiopen.
    dir2files = abilab.dir2abifiles(top, recurse=recurse)
    children = []
    for dirname, files in dir2files.items():
        w = widget_class(options=files, description="%s:" % dirname)
        # TODO: Should register the callback of "selected" but I didn't find the event type!
        w.observe(on_value_change, names='value', type="change")
        children.append(w)
    box = ipw.VBox(children=children)

    return display(box)
