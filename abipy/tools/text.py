# coding: utf-8
"""Utilities for working with strings and text."""
from __future__ import print_function, division, unicode_literals, absolute_import


def tonumber(s):
    """
    Convert string to number, raise ValueError if s cannot be converted.
    """
    # Duck test.
    try:
        stnum = s.upper().replace("D", "E")  # D-01 is not recognized by python: Replace it with E.
        # stnum = strip_punct(stnum)           # Remove punctuation chars.
        return float(stnum)                  # Try to convert.

    except ValueError:
        raise

    except:
        raise RuntimeError("Don't know how to handle type %s: %s" % (type(s), str(s)))


def nums_and_text(line):
    """
    Split line into (numbers, text).
    """
    tokens = line.split()
    text = ""
    numbers = []

    for tok in tokens:
        try:
            numbers.append(tonumber(tok))
        except ValueError:
            text += " " + tok

    return numbers, text
