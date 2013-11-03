#!/usr/bin/env python

import sys

def main():
    try:
        import netCDF4
        return 0
    except ImportError:
        raise 


if __name__ == '__main__':
    sys.exit(main())
