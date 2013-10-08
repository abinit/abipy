"""Utilities for working with strings and text."""
from __future__ import print_function, division

import sys
import os

from string import maketrans, punctuation
from pymatgen.util.string_utils import is_string, list_strings, pprint_table


_TABLE = maketrans("", "")


def strip_punct(s):
    """Remove punctuation characters from string s."""
    return s.translate(_TABLE, punctuation)


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


#def short_path(path, max_len=50):
#    """Returns the short version of path with maximum max_len characters."""
#    if len(path) <= max_len:
#        return path
#
#    tokens = path.split(os.path.sep)
#
#    l = [tokens[-1]]
#    for tok in tokens[-2::-1]:
#        if len(tok) + sum([len(t) for t in tokens]) <= max_len:
#            l.append(tok)
#        else:
#            break
#
#    s = os.path.join(*l)
#    return s


def marquee(text="", width=78, mark='*'):
    """
    Return the input string centered in a 'marquee'.

    :Examples:

    >>> marquee('A test', width=40)
    '**************** A test ****************'

    >>> marquee('A test', width=40, mark='-')
    '---------------- A test ----------------'

    marquee('A test',40, ' ')
    '                 A test                 '
    """
    if not text:
        return (mark*width)[:width]

    nmark = (width-len(text)-2)//len(mark)//2
    if nmark < 0: 
        nmark = 0

    marks = mark * nmark
    return '%s %s %s' % (marks, text, marks)
