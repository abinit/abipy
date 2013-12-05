#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import argparse 
import collections

import abipy.gui.wxapps as wxapps 


def str_examples():
    examples = (
      "\n"
      "Usage example:\n\n" 
      "abiopen.py file foo_WFK.nc ==> Visualize all WFK files foo.\n"
      "abiopen.py list dirpath    ==> Visualize all netcdf files (*.nc) in the directory dirpath (flat list mode) .\n"
      "abiopen.py tree dirpath    ==> Visualize all netcdf files (*.nc) in the directory dirpath (tree mode).\n"
      "abiopen.py scan dirpath    ==> Scan the given directory.\n"
    )
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")
    sys.exit(error_code)


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    #parser.add_argument('-f', '--filepaths', nargs="+", default=None, help="List of files.")
    #parser.add_argument("dirpaths", nargs="*", help="List of directories.")

    # Create the parsers for the sub-commands.
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_files = subparsers.add_parser('files', help="Run the specified file(s).") 
    p_files.add_argument("filepaths", nargs="+", help="List of filepaths.")

    # Subparser for list command.
    p_list = subparsers.add_parser('list', help="List files in directory.") 
    p_list.add_argument("dirpaths", nargs="+", help="List of filepaths.")
    p_list.add_argument('-w', '--wildcard', type=str, default="*.nc", help="wildcards. Default *.nc")

    # Subparser for tree command.
    p_tree = subparsers.add_parser('tree', help="Show files in directory tree.") 
    p_tree.add_argument("dirpaths", nargs="+", help="List of filepaths.")
    p_tree.add_argument('-w', '--wildcard', type=str, default="*.nc", help="wildcards. Default *.nc")

    # Subparser for scan command.
    p_scan = subparsers.add_parser('scan', help="Show files in directory tree.") 
    p_scan.add_argument("top", nargs="+", help="Top.")
    p_scan.add_argument('-w', '--wildcard', type=str, default="*.nc", help="wildcards. Default *.nc")
    #p_scan.add_argument('-e', '--exclude', type=str, default="*.nc", help="wildcards. Default *.nc")
    p_scan.add_argument('-r', '--recurse', type=bool, default=False, help="Recursive mode.")
    #p_scan.add_argument('-t', '--type', type=string, default="", help="Recursive mode.")

    # Parse the command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if options.command == "list":
        if not options.dirpaths: options.dirpaths = [os.getcwd()]

        wxapps.wxapp_listbrowser(dirpaths=options.dirpaths, 
                                 #filepaths=options.filepaths, 
                                 wildcard=options.wildcard).MainLoop()

    elif options.command == "tree":
        if not options.dirpaths: options.dirpaths = [os.getcwd()]
        wxapps.wxapp_dirbrowser(dirpath=options.dirpaths[0]).MainLoop()

    elif options.command == "files":
        # Dictionary mapping WX application classes to list of files to open.
        acls2files = collections.defaultdict(list)
        for path in options.filepaths:
            app_class = wxapps.file2appcls(path)
            acls2files[app_class].append(path)
        #print(acls2files)

        for cls, files in acls2files.items():
            if cls is None:
                print("Got None for files %s" % files)
            else:
                cls(files).MainLoop()

    elif options.command == "scan":
        filepaths = select_files(options)


    return 0


def select_files(options):
    wildcard = WildCard(options.wildcard)
    filepaths = []

    if options.recurse:
        for root, dirnames, filenames in os.walk(options.root):
            fnames = [os.path.join(root, f) for f in filenames]
            filepaths += wildcard.filter(fnames)
    else:
        # Select only the files in root.
        fnames = map(os.path.isfile, os.listdir(options.root))
        filepaths += wildcard.filter(fnames)

    return filepaths


if __name__ == "__main__":
    sys.exit(main())
