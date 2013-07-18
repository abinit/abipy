"""This module contains useful decorators for a variety of functions."""
from __future__ import print_function, division

import functools

# Import all pymatgen decorators in namespace.
from pymatgen.util.decorators import *

def benchmark(func):
    """
    A decorator that prints the time a function takes to execute.
    """
    import time
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t = time.clock()
        res = func(*args, **kwargs)
        print(func.__name__, time.clock()-t)
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
        func(*args, **kwargs)
        print("Exiting ", func.__name__)
    return new_function
