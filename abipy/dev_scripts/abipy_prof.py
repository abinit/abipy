#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
from os.path import exists, join as pj


def abipy_import():
    "Import the abipy package"
    #print sys.path
    import abipy 
    print("Imported abipy package", abipy.__file__)


def main(top):
    top = os.path.abspath(top)
    if not exists(pj(top, "abipy")) or not exists(pj(top, "setup.py")): 
        raise ValueError("top %s is not the top-level abipy directory" % top)

    os.chdir(top)

    sys.path.insert(0, top)

    import cProfile
    import pstats

    prof_file = pj(top, ".abipy_prof")

    # Profile the import of the package
    cProfile.run('abipy_import()', prof_file)

    p = pstats.Stats(prof_file)
    #If you were looking to see what functions were looping a lot, and taking a lot of time, you would do:
    p.sort_stats('time').print_stats(10)

    #This will sort all the statistics by file name, and then print out statistics for only the class init methods 
    #(since they are spelled with
    #__init__ in them). 
    #p.sort_stats('file').print_stats('__init__')

    return 0

if __name__ == '__main__':
    top = sys.argv[1]
    sys.exit(main(top))
