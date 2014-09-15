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


def touch(filename):
    try:
        open(fname,"w").close()
        return 0
    except:
        raise IOError("trying to create file = %s" % filename)


def tail_file(fname, n, aslist=False):
    """Assumes a unix-like system."""

    args = ["tail", "-n " + str(n), fname]

    p = Popen(args, shell=False, stdout=PIPE, stderr=PIPE)
    ret_code = p.wait()

    if ret_code != 0:
        raise RuntimeError("return_code = %s, cmd = %s" %(ret_code, " ".join(args) ))

    if aslist:
        return p.stdout.readlines()
    else:
        return p.stdout.read()


def which(program):
    """
    python version of the unix tool which locates a program file in the user's path

    Returns None if program cannot be found.
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

