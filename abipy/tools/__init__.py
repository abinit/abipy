# coding: utf-8
"""Helper functions."""
from __future__ import print_function, division, absolute_import, unicode_literals

from .devtools import *
from .iotools  import *
from .numtools import *
from .text import *


class NoDefaultProvided(object):
    pass


def hasattrd(obj, name):
    """
    The arguments are an object and a string.
    The result is True if the string is the name of one of the object’s attributes, False if not.
    Unlike the builtin hasattr, hasattrd supports dot notation e.g. hasattr(int, "__class__.__name__")
    (This is implemented by calling getattrd(object, name) and seeing whether it raises an exception or not.)
    """
    try:
        getattrd(obj, name)
        return True
    except AttributeError:
        return False


def getattrd(obj, name, default=NoDefaultProvided):
    """
    Same as getattr(), but allows dot notation lookup e.g. getattrd(obj, "a.b")

    Raises: AttributeError if ``name`` is not found and ``default`` is not given.

    Discussed in: http://stackoverflow.com/questions/11975781
    """
    from functools import reduce
    try:
        return reduce(getattr, name.split("."), obj)
    except AttributeError:
        if default is not NoDefaultProvided:
            return default
        raise
