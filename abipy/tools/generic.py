# coding: utf-8
"""Generic helper functions"""
from __future__ import print_function, division, unicode_literals

import os

from subprocess import Popen, PIPE


def pipe_commands(args1, args2):
    p1 = Popen(args1, stdout=PIPE)
    p2 = Popen(args2, stdin=p1.stdout, stdout=PIPE)
    return p2.communicate() # (stdout, stderr)


def unzip(gz_fname, dest=None):
    """decompress a gz file."""
    import gzip
    if not gz_fname.endswith(".gz"):
        raise ValueError("%s should end with .gz" % gz_fname)

    try:
        gz_fh = gzip.open(gz_fname, 'rb')
        file_content = gz_fh.read()
    finally:
        gz_fh.close() # Cannot use try, except, finally in python2-4

    try:
        if dest is None: dest = gz_fname[:-3]
        out_fh = open(dest, "wb")
        out_fh.write(file_content)
    finally:
        out_fh.close()
