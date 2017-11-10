# coding: utf-8
"""Decorators."""
from __future__ import print_function, division, unicode_literals, absolute_import

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
