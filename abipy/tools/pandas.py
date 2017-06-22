"""Utilities for pandas dataframe"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import pandas as pd


def print_frame(frame, title=None, sortby=None, file=sys.stdout):
    """
    Print entire pandas DataFrame.

    Args:
        frame: pandas DataFrame.
        title: Optional string to print as initial title.
        sortby: string name or list of names which refer to the axis items to be sorted (dataframe is not changed)
        file: a file-like object (stream); defaults to the current sys.stdout.
    """
    if title is not None: print(title, file=file)
    if sortby is not None:
        frame = frame.sort_values(sortby, inplace=False)
    with pd.option_context('display.max_rows', len(frame),
                           'display.max_columns', len(list(frame.keys()))):
        print(frame, file=file)
        print(" ", file=file)
