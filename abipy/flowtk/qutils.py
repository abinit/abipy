# coding: utf-8
"""
Collection of low-level tools to facilitate the interface with resource managers.

The preferred way of importing this module is:

    import qutils as qu
"""
from __future__ import annotations

import os
import json

from monty.string import is_string
from pymatgen.core.units import Time, Memory
from abipy.tools import duck


def slurm_parse_timestr(s: str) -> Time:
    """
    A slurm time parser. Accepts a string in one the following forms:

        # "days-hours",
        # "days-hours:minutes",
        # "days-hours:minutes:seconds".
        # "minutes",
        # "minutes:seconds",
        # "hours:minutes:seconds",

    Returns: Time in seconds.

    Raises:
        `ValueError` if string is not valid.
    """
    days, hours, minutes, seconds = 0, 0, 0, 0

    if duck.is_number_like(s):
        return Time(s, "s")

    if '-' in s:
        # "days-hours",
        # "days-hours:minutes",
        # "days-hours:minutes:seconds".
        days, s = s.split("-")
        days = int(days)

        if ':' not in s:
            hours = int(float(s))
        elif s.count(':') == 1:
            hours, minutes = map(int, s.split(':'))
        elif s.count(':') == 2:
            hours, minutes, seconds = map(int, s.split(':'))
        else:
            raise ValueError("More that 2 ':' in string!")

    else:
        # "minutes",
        # "minutes:seconds",
        # "hours:minutes:seconds",
        if ':' not in s:
            minutes = int(float(s))
        elif s.count(':') == 1:
            minutes, seconds = map(int, s.split(':'))
        elif s.count(':') == 2:
            hours, minutes, seconds = map(int, s.split(':'))
        else:
            raise ValueError("More than 2 ':' in string!")

    return Time((days*24 + hours)*3600 + minutes*60 + seconds, "s")


def time2slurm(timeval: float, unit="s") -> str:
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the slurm convention: "days-hours:minutes:seconds".

    >>> assert time2slurm(61) == '0-0:1:1' and time2slurm(60*60+1) == '0-1:0:1'
    >>> assert time2slurm(0.5, unit="h") == '0-0:30:0'
    """
    d, h, m, s = 24*3600, 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    days, hours = divmod(timeval, d)
    hours, minutes = divmod(hours, h)
    minutes, secs = divmod(minutes, m)

    return "%d-%d:%d:%d" % (days, hours, minutes, secs)


def time2pbspro(timeval: float, unit="s") -> str:
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the PbsPro convention: "hours:minutes:seconds".

    >>> assert time2pbspro(2, unit="d") == '48:0:0'
    """
    h, m, s = 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    hours, minutes = divmod(timeval, h)
    minutes, secs = divmod(minutes, m)

    return "%d:%d:%d" % (hours, minutes, secs)


def time2loadlever(timeval: float, unit="s") -> str:
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the LoadLever convention. format hh:mm:ss (hours:minutes:seconds)

    >>> assert time2loadlever(2, unit="d") == '48:00:00'
    """
    h, m, s = 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    hours, minutes = divmod(timeval, h)
    minutes, secs = divmod(minutes, m)

    return "%d:%02d:%02d" % (hours, minutes, secs)


def timelimit_parser(s):
    """Convert a float or a string into time in seconds."""
    try:
        return Time(float(s), "s")
    except ValueError:
        return slurm_parse_timestr(s)


def any2mb(s):
    """Convert string or number to memory in megabytes."""
    if is_string(s):
        return int(Memory.from_str(s).to("Mb"))
    else:
        return int(s)


def slurm_get_jobs(username=None) -> dict[int, dict]:
    """
    Invoke squeue, parse output and return list of dictionaries with job info indexed by job id.
    """
    # Based on https://gist.github.com/stevekm/7831fac98473ea17d781330baa0dd7aa
    username = os.getlogin() if username is None else username
    import subprocess as sp
    process = sp.Popen(['squeue', '-u',  username, "-o", '%all'],
                       stdout=sp.PIPE, stderr=sp.PIPE, shell=False, universal_newlines=True)
    proc_stdout, proc_stderr = process.communicate()

    lines = proc_stdout.split('\n')
    header_line = lines.pop(0)
    header_cols = header_line.split('|')
    entries = []
    error_lines = [] # do something with this later
    for line in lines:
        parts = line.split('|')
        if len(parts) != len(header_cols):
            error_lines.append((len(parts), line, parts))
        else:
            d = {}
            for i, key in enumerate(header_cols):
                d[key] = parts[i]
                if key == "JOBID":
                    d[key] = int(d[key])
            entries.append(d)

    return {e["JOBID"]: e for e in entries}
