#!/usr/bin/env python
# coding: utf-8
"""
This script tests the parser implemented in AbinitInputFile.
It tries to read and parse all the input files of the Abinit test suite."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from abipy import abilab

def is_abinit_input(path):
    """True if path is one of the input files used in the Test suite."""
    #beg = "#%%<BEGIN TEST_INFO>"
    #end = "#%%<END TEST_INFO>"
    if not path.endswith(".in"): return False

    with open(path, "rt") as fh:
        for line in fh:
            if "executable" in line and "abinit" in line: return True
        return False


def test_abinit_input_parser(top):

    inputs = []
    nfiles, retcode = 0, 0
    for dirpath, dirnames, filenames in os.walk(top):
        for fname in filenames:
            path = os.path.join(dirpath, fname)
            if not is_abinit_input(path): continue
            nfiles += 1
            try:
                inp = abilab.AbinitInputFile.from_file(path)
            except Exception as exc:
                retcode += 1
                print("[%s]: Exception:\n%s" % (path, exc))

    print("failed: %d/%d [%.1f%%]" % (retcode, nfiles, 100 * retcode/nfiles))
    return retcode


if __name__ == "__main__":
    import sys
    try:
        top = sys.argv[1]
    except IndexError:
        print("Usage: test_abinit_parser.py directory_with_input_files")
        sys.exit(1)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, "DEBUG", None)
    #if not isinstance(numeric_level, int):
    #    raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    sys.exit(test_abinit_input_parser(top))
