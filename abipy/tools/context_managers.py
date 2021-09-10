"""Context managers"""
from __future__ import annotations

import time

from contextlib import contextmanager


class Timer:
    """
    Context manager to time code section.

    .. example::

        with Timer("Timing code section"):
            do_stuff()
    """

    def __init__(self, description: str):
        self.description = description

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, exc_type, value, traceback):
        self.end = time.time()
        print(f"{self.description}.\n\tCompleted in {self.end - self.start:.4f} seconds\n")


@contextmanager
def temporary_change_attributes(something, **kwargs):
    """
    https://stackoverflow.com/questions/38531851/how-to-assign-member-variables-temporarily/38532086

    Usage:

        class Something(object):
            def __init__(self, x, y):
                self.x = x
                self.y = y

            def say_hello(self):
                print("hello", self.x, self.y)

        s = Something(1, 2)
        s.say_hello()  # hello 1 2
        with temporary_change_attributes(s, x=4, y=5):
            s.say_hello()  # hello 4 5
        s.say_hello()  # hello 1 2
    """
    previous_values = {k: getattr(something, k) for k in kwargs}
    for k, v in kwargs.items():
        setattr(something, k, v)
    try:
        yield
    finally:
        for k, v in previous_values.items():
            setattr(something, k, v)
