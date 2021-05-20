"""Utilities for pandas dataframes"""

import sys


def print_dataframe(frame, title=None, precision=6, sortby=None, file=sys.stdout, display=None):
    """
    Print entire pandas DataFrame.

    Args:
        frame: pandas DataFrame.
        title: Optional string to print as initial title.
        precision: Floating point output precision (number of significant digits).
            This is only a suggestion [default: 6] [currently: 6]
        sortby: string name or list of names which refer to the axis items to be sorted (dataframe is not changed)
        file: a file-like object (stream); defaults to the current sys.stdout.
            If file == "string", a temporary stream is created and a string is returned.
        display: Use ipython rich display protocol by invoking _repr_`display_ and returning the result.
            Use e.g. display="html" to get HTML table.
    """
    return_string = file == "string"
    if return_string:
        from io import StringIO
        file = StringIO()

    if title is not None: print(title, file=file)
    if sortby is not None and sortby in frame:
        frame = frame.sort_values(sortby, inplace=False)

    import pandas as pd
    with pd.option_context("display.max_rows", len(frame),
                           "display.max_columns", len(list(frame.keys())),
                           "display.precision", precision,
                           ):
        if display is None:
            print(frame, file=file)
            print(" ", file=file)
            if return_string: return file.getvalue()
        else:
            from IPython.core.display import HTML
            output = getattr(frame, "_repr_%s_" % display)()
            return HTML(output)
