#!/usr/bin/env python
"""This script extracts the docstrings from the run_*.py scripts located in this directory."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os

from monty.fnmatch import WildCard

def main():
    # Find (runnable) scripts.
    dir = os.path.dirname(os.path.abspath(__file__))
    wildcard = WildCard("run_*.py")
    scripts = [f.replace(".py", "") for f in wildcard.filter(os.listdir(dir))]

    retcode = 0

    for script in scripts:
        mod = __import__(script)
        if mod.__doc__ is None: retcode += 1
        doc = str(mod.__doc__).lstrip().rstrip()

        print(script + ":\n" + doc + "\n")

    #assert retcode == 0
    return retcode


if __name__ == "__main__":
    sys.exit(main())
