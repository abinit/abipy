# coding: utf-8
"""Decorators."""
from __future__ import print_function, division, unicode_literals, absolute_import

import time

from functools import wraps


def return_straceback_ifexc(func):
    """
    Decorator for functions that are supposed to return a string for logging purposes (e.g. str)
    Instead of raising an exception, the decorated function returns a string with the
    traceback so that execution can continue.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            import traceback
            return traceback.format_exc()
    return wrapper


def timeit(method):
    """
    timeit decorator adapted from:
    https://medium.com/pythonhive/python-decorator-to-measure-the-execution-time-of-methods-fa04cb6bb36d
    sets the timing of the routine as an attribute of the class
    """
    def timed(self, *args, **kw):
        ts = time.time()
        result = method(self, *args, **kw)
        te = time.time()

        setattr(self,"time_"+method.__name__, (te - ts) * 1000)
        return result
    return timed
