#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals
import sys

def main():
    """Check that the env on the local machine is properly setup"""
    try:
        from abipy.abilab import abicheck
        return abicheck()
    except:
        return 1

if __name__ == "__main__":
    sys.exit(main())
