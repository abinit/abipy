"""This module contains useful decorators for a variety of functions."""
from __future__ import print_function, division, unicode_literals

import sys
import os
import linecache
import functools

# Import all pymatgen decorators in namespace.
from pymatgen.util.decorators import *


def benchmark(func):
    """
    A decorator that computes the time a function takes to execute 
    and stores it in the __etime attribute.
    """
    import time

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t = time.clock()
        res = func(*args, **kwargs)
        etime = time.clock()-t
        print(func.__name__, etime)
        func.__etime = etime
        return res
    return wrapper


def logging(func):
    """
    A decorator that logs the activity of the script.
    (it actually just prints it, but it could be logging!)
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        print(func.__name__, args, kwargs)
        return res
    return wrapper


def verbose(func):
    """
    Decorator for functions. Returns a new function that prints 
    a message when func starts and finishes
    """
    @functools.wraps(func)
    def new_function(*args, **kwargs):
        print("Entering", func.__name__)
        res = func(*args, **kwargs)
        print("Exiting ", func.__name__)
        return res
    return new_function


def linetracing(f):
    """Trace a function (line by line)."""
    def globaltrace(frame, why, arg):
        if why == "call":
            return localtrace
        return None

    def localtrace(frame, why, arg):
        if why == "line":
            # record the file name and line number of every trace
            filename = frame.f_code.co_filename
            lineno = frame.f_lineno

            bname = os.path.basename(filename)
            print("{}({}): {}".format(bname,
                                      lineno,
                                      linecache.getline(filename, lineno)))
        return localtrace

    def _f(*args, **kwds):
        sys.settrace(globaltrace)
        result = f(*args, **kwds)
        sys.settrace(None)
        return result

    return _f
