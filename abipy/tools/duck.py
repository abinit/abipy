# coding: utf-8
"""Duck-typing tests"""
import collections
import warnings
import numpy as np


def is_string(s):
    """True if s behaves like a string (duck typing test)."""
    try:
        s + " "
        return True
    except Exception:
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
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
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
    if not isinstance(obj, collections.abc.Sequence): return False
    if is_string(obj): return False

    try:
        obj[:0]
        return True
    except TypeError:
        return False


def list_ints(arg):
    """
    Always return a list of int, given a int or list of integers as input.

    :Examples:

    >>> list_ints(1)
    [1]
    """
    l = np.array(arg, dtype=int)
    return [int(l)] if l.size == 1 else l.tolist()


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
        except Exception:
            raise TypeError("Don't know how to convert %s into a range object" % str(obj))


def as_slice(obj):
    """
    Convert an integer, a string or a slice object into slice.

    >>> assert as_slice(5) == slice(5, 6, 1)
    >>> assert as_slice("[1:4]") == slice(1, 4, 1)
    >>> assert as_slice("1::2") == slice(1, None, 2)
    """
    if isinstance(obj, slice) or obj is None: return obj

    try:
        # integer.
        if int(obj) == float(obj): return slice(int(obj), int(obj)+1, 1)
    except Exception:
        # assume string defining a python slice [start:stop:step]
        if not obj: return None
        if obj.count("[") + obj.count("]") not in (0, 2):
            raise ValueError("Invalid string %s" % obj)

        obj = obj.replace("[", "").replace("]", "")
        n = obj.count(":")
        if n == 0:
            obj = int(obj)
            return slice(obj, obj+1)

        tokens = [int(f) if f else None for f in obj.split(":")]
        if len(tokens) == 2: tokens.append(1)
        if tokens[2] is None: tokens[2] = 1

        return slice(*tokens)

    raise ValueError("Cannot convert %s into a slice:\n%s" % (type(obj), obj))


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
