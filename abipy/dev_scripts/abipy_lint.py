#!/usr/bin/env python
import sys
import os
from os.path import exists, join as pj

from subprocess import call

def abipy_reindent(top)
    """
    Change Python (.py) files to use 4-space indents and no hard tab characters.
    Also trim excess spaces and tabs from ends of lines, and remove empty lines
    at the end of files.  Also ensure the last line ends with a newline.
    No backup is done
    """
    reindent_script = pj(top, "tools", "reindent.py")

    return call([reindent_script,  top, "-rnv"])


def main(top):
    top = os.path.abspath(top)
    if not exists(pj(top, "abipy")) or not exists(pj(top, "setup.py")): 
        raise ValueError("top %s is not the top-level abipy directory" % top)

    abipy_reindent(top)

    return 0

if __name__ == '__main__':
    top = sys.argv[1]
    sys.exit(main(top))
