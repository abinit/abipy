#!/usr/bin/env python
"""Script to export/visualize the crystal structure saved in the netcdf files produced by ABINIT."""
from __future__ import division, print_function

import sys
import os
import warnings
import argparse

from abipy import abilab


def str_examples():
    examples = """
Usage example:\n
    abistruct.py NCFILE convert cif    => Read the structure from the netcdf FILE.nc and create the CIF file.
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

    # Subparser for visualize command.
    #p_visualize = subparsers.add_parser('visualizeshot', help="Visualize the structure with the specified visualizer")

    #p_visualize.add_argument('visualizer', nargs="?", default="xcrysden", type=str, help="Visualizer.")

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
