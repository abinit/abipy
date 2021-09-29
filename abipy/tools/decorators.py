# coding: utf-8
"""Decorators."""

import time
import functools
import weakref
import inspect

from textwrap import dedent


def return_straceback_ifexc(func):
    """
    Decorator for functions that are supposed to return a string for logging purposes (e.g. str)
    Instead of raising an exception, the decorated function returns a string with the
    traceback so that execution can continue.
    """
    @functools.wraps(func)
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

        setattr(self,"time_" + method.__name__, (te - ts) * 1000)
        return result
    return timed


def memoized_method(*lru_args, **lru_kwargs):
    """
    Implements lru_cache for class methods.
    It takes the exact same parameters as lru_cache, and works exactly the same.
    However it never passes self to lru_cache and instead uses a per-instance lru_cache.

    Taken from: https://stackoverflow.com/questions/33672412/python-functools-lru-cache-with-class-methods-release-object


    ... example::

            @memoized_method(maxsize=12, typed=False)
            def method(self, a, b):
                ....
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapped_func(self, *args, **kwargs):
            # We're storing the wrapped method inside the instance. If we had
            # a strong reference to self the instance would never die.
            self_weak = weakref.ref(self)

            @functools.wraps(func)
            @functools.lru_cache(*lru_args, **lru_kwargs)
            def cached_method(*args, **kwargs):
                return func(self_weak(), *args, **kwargs)
            setattr(self, func.__name__, cached_method)
            return cached_method(*args, **kwargs)
        return wrapped_func
    return decorator


def dump_args(func):
    """
    Decorator to print function call details.

    This includes parameters names and effective values.
    """

    def wrapper(*args, **kwargs):
        func_args = inspect.signature(func).bind(*args, **kwargs).arguments
        func_args_str = ", ".join(
            "{} = {!r}".format(*item) for item in func_args.items()
        )
        print(f"{func.__module__}.{func.__qualname__} ( {func_args_str} )")
        return func(*args, **kwargs)

    return wrapper


class Appender(object):
    r"""
    A function decorator that appends an addendum to the docstring of the target function.
    This decorator should be robust even if func.__doc__ is None
    (for example, if -OO was passed to the interpreter).
    Usage: construct a docstring.Appender with a string to be joined to
    the original docstring. An optional 'join' parameter may be supplied
    which will be used to join the docstring and addendum. e.g.

    add_copyright = Appender("Copyright (c) 2009", join='\n')

    @add_copyright
    def my_dog(has='fleas'):
        "This docstring will have a copyright below"
        pass

    MG took it from:
    https://github.com/pandas-dev/pandas/blob/3a7f956c30528736beaae5784f509a76d892e229/pandas/util/_decorators.py#L156

    MG: Added dedent and debug args.
    """

    def __init__(self, addendum, join='', indents=0, dedent=True, debug=False):
        if indents > 0:
            self.addendum = indent(addendum, indents=indents)
        else:
            self.addendum = addendum
        self.join = join

        # MG additional stuff
        self.dedent = dedent
        self.debug = debug

    def __call__(self, func):
        func.__doc__ = func.__doc__ if func.__doc__ else ''
        self.addendum = self.addendum if self.addendum else ''

        if self.dedent:
            docitems = [dedent(func.__doc__), dedent(self.addendum)]
        else:
            docitems = [func.__doc__, self.addendum]

        if self.debug:
            print("func.__doc__:\n", func.__doc__)
            print("addendum:\n", self.addendum)

        func.__doc__ = dedent(self.join.join(docitems))

        if self.debug:
            print("new func.__doc__:\n", func.__doc__)

        return func


def indent(text, indents=1):
    if not text or not isinstance(text, str):
        return ''
    jointext = ''.join(['\n'] + ['    '] * indents)
    return jointext.join(text.split('\n'))

