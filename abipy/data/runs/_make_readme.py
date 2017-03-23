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

    missing = []
    with open("README.md", "wt") as fh:
        for script in scripts:
            mod = __import__(script)
            if mod.__doc__ is None: missing.append(script)
            doc = str(mod.__doc__).lstrip().rstrip()
            print("`%s`" % script + ":\n" + doc + "\n", file=fh)

    if missing:
        raise RuntimeError("The following script do not provide a doc string:\n" + "\n".join(missing))

    return len(missing)


if __name__ == "__main__":
    sys.exit(main())
