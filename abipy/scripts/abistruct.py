#!/usr/bin/env python
"""Script to export/visualize the crystal structure saved in the netcdf files produced by ABINIT."""
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse

from pprint import pprint
from abipy import abilab
from abipy.iotools.visualizer import Visualizer

def main():

    def str_examples():
        examples = """\
Usage example:\n

    abistruct.py convert filepath cif         => Read the structure from file and print CIF file.
    abistruct.py visualize filepath xcrysden  => Visualize the structure with XcrysDen.
    abistruct.py abivars filepath             => Print the ABINIT variables defining the structure.
    abistruct.py pmgdata mp-149               => Get structure from pymatgen database and print its JSON representation.
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for commands that need to know the filepath
    path_selector = argparse.ArgumentParser(add_help=False)
    path_selector.add_argument('filepath', nargs="?", help="File with the crystalline structure")

    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands, use command --help for help")

    # Subparser for convert command.
    p_convert = subparsers.add_parser('convert', parents=[path_selector], help="Convert structure to the specified format.")
    p_convert.add_argument('format', nargs="?", default="cif", type=str, help="Format of the output file (cif, cssr, POSCAR, json, mson).")

    # Subparser for abivars command.
    p_convert = subparsers.add_parser('abivars', parents=[path_selector], help="Print the ABINIT variables defining the structure.")

    # Subparser for visualize command.
    p_visualize = subparsers.add_parser('visualize', parents=[path_selector], help="Visualize the structure with the specified visualizer")
    p_visualize.add_argument('visualizer', nargs="?", default="xcrysden", type=str, help=("Visualizer name. "
        "List of visualizer supported: %s" % ", ".join(Visualizer.all_visunames())))

    # Subparser for pmgid command.
    p_pmgdata = subparsers.add_parser('pmgdata', help="Get structure from the pymatgen database. Requires internet connection and MAPI_KEY")
    p_pmgdata.add_argument("pmgid", type=str, default=None, help="Pymatgen identifier")
    p_pmgdata.add_argument("--mapi-key", default=None, help="Pymatgen MAPI_KEY. Use env variable if not specified.")
    p_pmgdata.add_argument("--host", default="www.materialsproject.org", help="Pymatgen database.")

    # Parse command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    # Read structure form file
    if hasattr(options, "filepath"):
        structure = abilab.Structure.from_file(options.filepath)

    if options.command == "convert":
        s = structure.convert(format=options.format)
        #print((" Abinit --> %s " % format).center(80, "*"))
        print(s)

    elif options.command == "abivars":
        print(structure.abi_string)

    elif options.command == "visualize":
        structure.visualize(options.visualizer)

    elif options.command == "pmgdata":
        # Get the Structure corresponding the a material_id.
        structure = abilab.Structure.from_material_id(options.pmgid, final=True, api_key=options.mapi_key, 
            host=options.host)

        # Convert to json and print it.
        s = structure.convert(format="json")
        #s = structure.convert(format="mson")
        print(s)

    else:
        raise ValueError("Unsupported command: %s" % options.command)

    return 0

if __name__ == "__main__":
    sys.exit(main())
