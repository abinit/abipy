"""Utilities for working with strings and text."""
from __future__ import print_function, division

import sys
import os

from pymatgen.util.string_utils import is_string, list_strings, pprint_table, WildCard


#from string import maketrans, punctuation
#_TABLE = maketrans("", "")
#
#def strip_punct(s):
#    """Remove punctuation characters from string s."""
#    return s.translate(_TABLE, punctuation)


def tonumber(s):
    """Convert string to number, raise ValueError if s cannot be converted."""
    # Duck test.
    try:
        stnum = s.upper().replace("D", "E")  # D-01 is not recognized by python: Replace it with E.
        stnum = strip_punct(stnum)           # Remove punctuation chars.
        return float(stnum)                  # Try to convert.

    except ValueError:
        raise

    except:
        raise RuntimeError("Don't know how to handle string: %s" % s)


def nums_and_text(line):
    """split line into (numbers, text)."""
    tokens = line.split()
    text = ""
    numbers = []

    for tok in tokens:
        try:
            numbers.append(tonumber(tok))
        except ValueError:
            text += " " + tok

    return numbers, text
