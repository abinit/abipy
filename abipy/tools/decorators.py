# coding: utf-8
"""Decorators."""

import time
import functools
import weakref

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

        setattr(self,"time_"+method.__name__, (te - ts) * 1000)
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
