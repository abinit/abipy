#!/usr/bin/env python
"""Script to export/visualize the crystal structure saved in the netcdf files produced by ABINIT."""
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse

from abipy import abilab
from abipy.iotools.visualizer import Visualizer

def main():

    def str_examples():
        examples = """\
Usage example:\n

    abistruct.py filepath convert cif         => Read the structure from file and print CIF file.
    abistruct.py filepath visualize xcrysden  => Visualize the structure with XcrysDen.
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)


    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('filepath', nargs="?", help="Netcdf file with the crystalline structure")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands, use command --help for help")

    # Subparser for convert command.
    p_convert = subparsers.add_parser('convert', help="Convert structure to the specified format.")
    p_convert.add_argument('format', nargs="?", default="cif", type=str, help="Format of the output file (cif, cssr, POSCAR, json, mson).")

    # Subparser for visualize command.
    p_visualize = subparsers.add_parser('visualize', help="Visualize the structure with the specified visualizer")

    p_visualize.add_argument('visualizer', nargs="?", default="xcrysden", type=str, help=("Visualizer name. "
        "List of visualizer supported: %s" % ", ".join(Visualizer.all_visunames())))

    # Parse command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    structure = abilab.Structure.from_file(options.filepath)

    if options.command == "convert":
        s = structure.convert(format=options.format)
        #print((" Abinit --> %s " % format).center(80, "*"))
        print(s)

    elif options.command == "visualize":
        structure.visualize(options.visualizer)

    else:
        raise ValueError("Unsupported command %s" % options.command)

    return 0

if __name__ == "__main__":
    sys.exit(main())
