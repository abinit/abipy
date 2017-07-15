"""Utilities for pandas dataframe"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import pandas as pd


def print_frame(frame, title=None, precision=6, sortby=None, file=sys.stdout):
    """
    Print entire pandas DataFrame.

    Args:
        frame: pandas DataFrame.
        title: Optional string to print as initial title.
        precision: Floating point output precision (number of significant digits).
            This is only a suggestion [default: 6] [currently: 6]
        sortby: string name or list of names which refer to the axis items to be sorted (dataframe is not changed)
        file: a file-like object (stream); defaults to the current sys.stdout.
    """
    if title is not None: print(title, file=file)
    if sortby is not None and sortby in frame:
        frame = frame.sort_values(sortby, inplace=False)

    with pd.option_context("display.max_rows", len(frame),
                           "display.max_columns", len(list(frame.keys())),
                           "display.precision", precision,
                           ):
    #with pd.option_context('display.expand_frame_repr', False):
    #with pd.option_context('display.large_repr', 'truncate', 'display.max_columns', 0):
        print(frame, file=file)
        print(" ", file=file)
