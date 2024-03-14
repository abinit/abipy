"""Context managers"""
from __future__ import annotations

import signal
import sys

from contextlib import contextmanager
from time import perf_counter


class Timer:
    """
    Context manager to time code section.

    .. example::

        with Timer("Timing code section"):
            do_stuff()

    or

        with Timer(header=f"Begin ABINIT", footer="ABINIT GS") as timer:
            do_stuff()
    """
    def __init__(self, footer=None, header=None, file=sys.stdout):
        self.header = header
        self.footer = footer
        self.file = file

    def __enter__(self):
        self.time = perf_counter()
        if self.header is not None:
            print(self.header, file=self.file)
        return self

    def __str__(self):
        return self.readout

    def __exit__(self, type, value, traceback):
        self.time = perf_counter() - self.time
        self.readout = f'Time: {self.time:.3f} seconds'
        if self.footer is not None:
            msg = f'{self.footer} completed in {self.time:.3f} seconds'.lstrip()
            print(msg, file=self.file)


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


class Timeout:
    """
    Taken from https://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish/22348885#22348885
    """

    def __init__(self, seconds: int, message: str = 'Timeout'):
        self.seconds = int(seconds)
        if self.seconds <= 0:
            raise ValueError(f"seconds should be > 0 while it is: {self.seconds}")
        self.message = message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        # If seconds is zero, any pending alarm is canceled.
        signal.alarm(0)
