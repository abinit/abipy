"""Utilities for working with strings and text."""
from __future__ import print_function, division

import sys
import os

from string import maketrans, punctuation


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
            numbers.append( tonumber(tok) )
        except ValueError:
            text += " " + tok

    return numbers, text


#def short_path(path, max_len=50):
#    """Returns the shor version of path with maximum 50 characters."""
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


def pprint_table(table, out=sys.stdout, rstrip=False):
    """
    Prints out a table of data, padded for alignment
    Each row must have the same number of columns.

    Args:
        out:    
            Output stream (file-like object)
        table: 
            The table to print. A list of lists.
        rstrip: 
            if True, trailing withespaces are removed from the entries.
    """
    def max_width_col(table, col_idx):
        """Get the maximum width of the given column index"""
        return max([len(row[col_idx]) for row in table])

    if rstrip:
        for row_idx, row in enumerate(table):
            table[row_idx] = [c.rstrip() for c in row]

    col_paddings = []
    ncols = len(table[0])
    for i in range(ncols):
        col_paddings.append(max_width_col(table, i))

    for row in table:
        # left col
        out.write( row[0].ljust(col_paddings[0] + 1) )
        # rest of the cols
        for i in range(1, len(row)):
            col = row[i].rjust(col_paddings[i] + 2)
            out.write(col)
        out.write("\n")


def is_string(s):
    """True if s behaves as a string (duck typing test)."""
    try:
        dummy = s + " "
        return True

    except TypeError:
        return False


def list_strings(arg):
    """
    Always return a list of strings, given a string or list of strings as input.

    :Examples:

    >>> list_strings('A single string')
    ['A single string']

    >>> list_strings(['A single string in a list'])
    ['A single string in a list']

    >>> list_strings(['A','list','of','strings'])
    ['A', 'list', 'of', 'strings']
    """
    if isinstance(arg, str):
        return [arg]
    else:
        return arg


def marquee(txt='', width=78, mark='*'):
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
    if not txt:
        return (mark*width)[:width]

    nmark = (width-len(txt)-2)//len(mark)//2
    if nmark < 0: 
        nmark = 0

    marks = mark * nmark
    return '%s %s %s' % (marks, txt, marks)

