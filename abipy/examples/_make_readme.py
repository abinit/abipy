#!/usr/bin/env python
"""
This script extracts the docstrings from the run_*.py scripts located in this directory.
"""
import sys
import os

from monty.fnmatch import WildCard


def main():
    # Find (runnable) scripts.
    dirpath = os.path.dirname(os.path.abspath(__file__))
    wildcard = WildCard("run_*.py")
    scripts = [f.replace(".py", "") for f in wildcard.filter(os.listdir(dirpath))]

    missing = []
    with open("README.md", "wt") as fh:
        for script in scripts:
            mod = __import__(script)
            if mod.__doc__ is None: missing.append(script)
            doc = str(mod.__doc__).lstrip().rstrip()
            doc = doc.replace("\n", "\n    ")
            print("``%s``:\n\n    " % script + doc + "\n", file=fh)

    if missing:
        raise RuntimeError("The following script do not provide a doc string:\n" + "\n".join(missing))

    return len(missing)


if __name__ == "__main__":
    sys.exit(main())
