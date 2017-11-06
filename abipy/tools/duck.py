# coding: utf-8
"""Duck-typing tests"""
from __future__ import print_function, division, unicode_literals, absolute_import

import collections
import warnings
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
        # This to get rid of warnings about casting complex to real.
        if np.iscomplexobj(obj) and np.isreal(obj):
            return int(obj.real) == obj
        else:
            return int(obj) == obj
    except (ValueError, TypeError):
        return False
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


def torange(obj):
    """
    Convert obj into a range. Accepts integer, slice object  or any object
    with an __iter__ method. Note that an integer is converted into range(int, int+1)

    >>> list(torange(1))
    [1]
    >>> list(torange(slice(0, 4, 2)))
    [0, 2]
    >>> list(torange([1, 4, 2]))
    [1, 4, 2]
    """
    if is_intlike(obj):
        return range(obj, obj + 1)

    elif isinstance(obj, slice):
        start = obj.start if obj.start is not None else 0
        step = obj.step if obj.step is not None else 1
        return range(start, obj.stop, step)

    else:
        try:
            return obj.__iter__()
        except:
            raise TypeError("Don't know how to convert %s into a range object" % str(obj))
