#!/usr/bin/env python
"""
This script remove all pyc files and all __pycache__ directory.
"""
from __future__ import print_function, division, unicode_literals


import sys
import os
import shutil


def rm_pycaches(top):
    count = 0
    for dirpath, dirnames, filenames in os.walk(top):
        for d in dirnames:
            if not d ==  "__pycache__": continue
            path = os.path.join(dirpath, d)
            #print("Will remove %s" % path)
            shutil.rmtree(path)
            count += 1

    print("Removed %d __pycache__ directories" % count)


def rm_pycfiles(top):
    count = 0
    for dirpath, dirnames, filenames in os.walk(top):
        for f in filenames:
            if not f.endswith(".pyc"): continue
            path = os.path.join(dirpath, f)
            #print("Will remove %s" % path)
            os.remove(path)
            count += 1

    print("Removed %d .pyc files" % count)
    return count

if __name__ == "__main__":
    if len(sys.argv) == 1:
        top = "."
    elif len(sys.argv) == 2:
        top = sys.argv[1]
    else:
        raise ValueError("Usage: pyclean [top]")

    rm_pycaches(top)
    rm_pycfiles(top)