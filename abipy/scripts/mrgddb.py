#!/usr/bin/env python
"""This script wraps the Fortran executable mrgddb."""
from __future__ import division, print_function

import os
import sys
import argparse

from pymatgen.io.abinitio.wrappers import Mrgddb

##########################################################################################
# Helper functions.

def str_examples():
    examples = \
"""Example usage:\n
    mrgddb -o output_DDB file1_DDB file2_DDB 
    mrgddb -o output_DDB -d "String with description" file1_DDB file2_DDB 
    mrgddb -o output_DDB  *_DDB
"""
    return examples


def show_examples_and_exit(err_msg=None, error_code=0):
    "Display the usage of the script."
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


##########################################################################################


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='verbose, can be supplied multiple times to increase verbosity.')  

    parser.add_argument('-e', '--executable', metavar='STRING', type=str, default="mrgddb", help="Path to the mrgddb executable.")

    parser.add_argument('-d', '--description', metavar='STRING', type=str, default="No description available", 
                        help="Description added to the merged DDB file.")

    parser.add_argument('-o', '--out_ddb', metavar='STRING', type=str, default="", help="Name of the output DDB file")

    parser.add_argument('ddb_files', nargs="+", help="List of DDB files to merge.")

    # Parse the command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    if not options.out_ddb:
        raise ValueError("out_ddb must be specified")

    mrgddb = Mrgddb(executable=options.executable, verbose=options.verbose)

    mrgddb.merge(options.ddb_files, options.out_ddb, options.description, cwd=None)
                                  
    return 0

##########################################################################################

if __name__ == "__main__":
    sys.exit(main())
