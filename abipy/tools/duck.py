# coding: utf-8
"""Duck-typing tests"""
from __future__ import print_function, division, unicode_literals, absolute_import

import collections
import numpy as np


def is_string(s):
    """True if s behaves like a string (duck typing test)."""
    try:
        s + " "
        return True
    except TypeError:
        return False


def is_intlike(obj):
    """
    True if obj represents an integer (float such as 1.0 are included as well).
    """
    # isinstance(i, numbers.Integral)
    try:
        return int(obj) == obj
    except (ValueError, TypeError):
        return False


def is_number_like(obj):
    """True if obj represents a number."""
    try:
        obj - 1
        return True
    except TypeError:
        return False


def is_listlike(obj):
    #if isinstance(branch, (list, tuple, np.ndarray)):
    if isinstance(obj, np.ndarray): return True
    if not isinstance(obj, collections.Sequence): return False
    if is_string(obj): return False

    try:
        obj[:0]
        return True
    except TypeError:
        return False
