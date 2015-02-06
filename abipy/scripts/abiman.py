#!/usr/bin/env python
"""
Interface to the database of ABINIT input variables
"""
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse


def main():
    def str_examples():
        examples = """\
Usage example:\n
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                    help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('varname', help="ABINIT variable")

    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    from abipy.abilab import abinit_help
    abinit_help(options.varname)


if __name__ == "__main__":
    try:
        do_prof = sys.argv[1] == "prof"
        if do_prof: sys.argv.pop(1)
    except: 
        pass

    if do_prof:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()
    else:
        sys.exit(main())
