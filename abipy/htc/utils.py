"""Tools and helper functions for abinit calculations"""
from __future__ import print_function, division #, unicode_literals

import os.path
import collections

from copy import deepcopy
from itertools import chain
from monty.string import list_strings

##########################################################################################

class _ewc_tuple(collections.namedtuple("ewc_tuple", "errors, warnings, comments")):

    def tostream(self, stream):
        "Return a string that can be visualized on stream (with colors if stream support them)."
        str_colorizer = StringColorizer(stream)

        red  = lambda num : str_colorizer(str(num), "red")
        blue = lambda num : str_colorizer(str(num), "blue")

        nums = map(len, [self.errors, self.warnings, self.comments])

        colors = (red, blue, str)

        for (i, n) in enumerate(nums):
            color = colors[i]
            nums[i] = color(n) if n else str(n)

        return "%s errors, %s warnings, %s comments in main output file" % tuple(nums)

##########################################################################################

def parse_ewc(filename, nafter=5):
    """
    Extract errors, warnings and comments from file filename.

    :arg nafter: Save nafter lines of trailing context after matching lines.
    :return: namedtuple instance. The lists of strings with the corresponding messages are
             available in tuple.errors, tuple.warnings, tuple.comments.
    """
    # TODO
    # we have to standardize the abinit WARNING, COMMENT and ERROR  so that we can parse them easily
    # without having to use nafter.

    errors, warnings, comments = [], [], []
    # Note the space after the name.
    exc_cases = ["ERROR ", "BUG ", "WARNING ", "COMMENT "]

    handlers = {
        "ERROR "   : errors.append,
        "BUG "     : errors.append,
        "WARNING " : warnings.append,
        "COMMENT " : comments.append,
    }

    def exc_case(line):
        for e in exc_cases:
            if e in line: return e
        else:
            return None

    with open(filename, "r") as fh:
        lines = fh.readlines()
        nlines = len(lines)
        for (lineno, line) in enumerate(lines):
            handle = handlers.get(exc_case(line))
            if handle is None: continue
            context = lines[lineno: min(lineno+nafter, nlines)]
            handle( "".join([c for c in context]) )

    return _ewc_tuple(errors, warnings, comments)

##########################################################################################

def find_file(files, ext, prefix=None, dataset=None, image=None):
  """
  Given a list of file names, return the file with extension "_" + ext, None if not found.

  The prefix, the dataset index and the image index can be specified

  .. warning::

     There are some border cases that will confuse the algorithm
     since the order of dataset and image is not tested.
     Solving this problem requires the knowledge of ndtset and nimages
     This code, however should work in 99.9% of the cases.
  """
  separator = "_"

  for filename in list_strings(files):
      # Remove Netcdf extension (if any)
      f = filename[:-3] if filename.endswith(".nc") else filename
      if separator not in f: continue
      tokens = f.split(separator)
      if tokens[-1] == ext:
        found = True
        if prefix is not None:  found = found and filename.startswith(prefix)
        if dataset is not None: found = found and "DS" +  str(dataset) in tokens
        if image is not None:   found = found and "IMG" + str(image)   in tokens
        if found: return filename
  else:
      return None

##########################################################################################

def abinit_output_iscomplete(output_file):
    "Returns True if the abinit output file is complete."
    if not os.path.exists(output_file):
        return False

    chunk = 5 * 1024  # Read only the last 5Kb of data.
    nlines = 10       # Check only in the last 10 lines.

    MAGIC = "Calculation completed."

    with open(output_file, 'r') as f:
        size = f.tell()
        f.seek(max(size - chunk, 0))
        try:
            for line in f.read().splitlines()[-nlines:]:
                if MAGIC in line:
                    return True
        except:
            pass

    return False


# =========================================================================== #

def is_number(s):
    """Returns True if the argument can be made a float."""
    try:
        float(s)
        return True
    except:
        return False

def is_iter(obj):
    """Return True if the argument is list-like."""
    return hasattr(obj, '__iter__')

def is_scalar(obj):
    """Return True if the argument is not list-like."""
    return not is_iter

def flatten(iterable):
    """Make an iterable flat, i.e. a 1d iterable object."""
    iterator = iter(iterable)
    array, stack = collections.deque(), collections.deque()
    while True:
        try:
            value = next(iterator)
        except StopIteration:
            if not stack:
                return tuple(array)
            iterator = stack.pop()
        else:
            if not isinstance(value, str) \
               and isinstance(value, collections.Iterable):
                stack.append(iterator)
                iterator = iter(value)
            else:
                array.append(value)

def listify(obj):
    """Return a flat list out of the argument."""
    if not obj:
        obj = list()
    elif is_iter(obj):
        obj = list(flatten(obj))
    else:
        obj = [obj]
    return deepcopy(obj)


class StringColorizer(object):
    COLORS = {"default": "",
               "blue": "\x1b[01;34m",
               "cyan": "\x1b[01;36m",
               "green": "\x1b[01;32m",
               "red": "\x1b[01;31m",
               # lighting colors.
               #"lred":    "\x1b[01;05;37;41m"
    }

    def __init__(self, stream):
        self.has_colours = stream_has_colours(stream)

    def __call__(self, string, colour):
        if self.has_colours:
            code = self.COLORS.get(colour, "")
            if code:
                return code + string + "\x1b[00m"
            else:
                return string
        else:
            return string
