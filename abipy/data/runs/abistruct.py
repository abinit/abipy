#!/usr/bin/env python
from __future__ import division, print_function

import sys
import os
import warnings
import argparse

from abipy import abilab


def str_examples():
    examples = """
Usage example:\n
    abirun.py DIRPATH  singleshot              => Fetch the first available task and run it.
    abirun.py DIRPATH  rapidfire               => Keep repeating, stop when no task can be executed
                                                  due to inter-dependency.
    abirun.py DIRPATH gui                      => Open the GUI 
    nohup abirun.py DIRPATH sheduler -s 30 &   => Use a scheduler to schedule task submission
"""
    return examples


def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")

    sys.exit(error_code)


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('filepath', nargs="?", help="Netcdf file with the crystalline structure")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for convert command.
    p_convert = subparsers.add_parser('convert', help="Convert structure to the specified format.")

    p_convert.add_argument('format', nargs="?", default="cif", type=str, help="Format of the output file (ciff, POSCAR, json).")

    # Subparser for single command.
    #p_single = subparsers.add_parser('singleshot', help="Run single task.") #, aliases=["single"])


    # Parse command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    structure = abilab.Structure.from_file(options.filepath)

    if options.command == "convert":
        print(options.format)
        #for format in ["cif", "POSCAR", "cssr", "json"]:
        s = structure.convert(format=options.format)

        #print((" Abinit --> %s " % format).center(80, "*"))
        print(s)

    else:
        raise ValueError("Unsupported command %s" % options.command)

    return 0

if __name__ == "__main__":
    sys.exit(main())
