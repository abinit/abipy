#!/usr/bin/env python

import sys

def main():
    """Check that the env on the local machine is properly setup"""
    try:
        return abicheck()
    except:
        return 1

if __name__ == "__main__":
    sys.exit(main())
