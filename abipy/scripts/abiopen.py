#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse 
import collections

from abipy.tools.text import WildCard 
import abipy.gui.wxapps as wxapps 

import logging
logger.getLogger(__name__)


def main():

    def str_examples():
        examples = (
          "\n"
          "Usage example:\n\n" 
          "abiopen.py files foo_WFK.nc          ==> Visualize the WFK file foo_WFK.nc\n"
          "                                         (many other Abinit files are supported, just try!).\n"
          "abiopen.py list dirpath              ==> Visualize all files in the directory dirpath (flat list mode) .\n"
          "abiopen.py tree dirpath              ==> Visualize all files in the directory dirpath (tree mode).\n"
          "abiopen.py scan dirpath              ==> Scan all the supported files in the given directory (recursive mode).\n"
          "abiopen.py scan dirpath -w *GSR.nc   ==> Walk the directory tree starting from dirpath \n"
          "                                         and open all the GSR.nc files encountered.\n"
        )
        return examples


    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: 
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)


    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands.
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_files = subparsers.add_parser('files', aliases=["file", "f"], help="Run the specified file(s).") 
    p_files.add_argument("filepaths", nargs="+", help="List of filepaths.")

    # Subparser for list command.
    p_list = subparsers.add_parser('list', aliases=["l"], help="List files in directory.") 
    p_list.add_argument("dirpaths", nargs="+", help="List of filepaths.")
    p_list.add_argument('-w', '--wildcard', type=str, default="*.nc", help="wildcards. Default *.nc")

    # Subparser for tree command.
    p_tree = subparsers.add_parser('tree', aliases=["t"], help="Show files in directory tree.") 
    p_tree.add_argument("dirpaths", nargs="+", help="List of filepaths.")
    p_tree.add_argument('-w', '--wildcard', type=str, default="*.nc", help="wildcards. Default *.nc")

    # Subparser for scan command.
    p_scan = subparsers.add_parser('scan', aliases=["s"], help="Show files in directory tree.") 
    p_scan.add_argument("top", help="Top.")
    p_scan.add_argument('-w', '--wildcard', type=str, default="*.nc", help="wildcards. Default *.nc")
    p_scan.add_argument('--no-walk', default=False, action="store_true", help="Disable walk mode.")
    #p_scan.add_argument('-e', '--extension', type=str, default="GSR.nc", help="File extension to search for")
    #p_scan.add_argument('-t', '--type', type=string, default="", help="Recursive mode.")

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    if options.verbose:
        print("*** Command line options *** ")
        print(options)

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
        acls2files = appclasses_from_files(options.filepaths)
        if options.verbose: print(acls2files)

        for cls, files in acls2files.items():
            cls(files).MainLoop()

    elif options.command == "scan":
        # Select the files to open.
        filepaths = select_files(options)

        # Dictionary mapping WX application classes to list of files to open.
        acls2files = appclasses_from_files(filepaths)
        if options.verbose: print(acls2files)

        for cls, files in acls2files.items():
            cls(files).MainLoop()

    return 0


def appclasses_from_files(filepaths):
    """Returns a dictionary mapping WX application classes to the list of files to open."""
    acls2files, bad_files = collections.defaultdict(list), []

    for path in filepaths:
        app_class = wxapps.file2appcls(path)
        if app_class is None:
            bad_files.append(path)
        else:
            acls2files[app_class].append(path)
    #print(acls2files)

    if bad_files:
        logger.warning("Cannot find wx application for files:\n%s" % bad_files)

    return acls2files


def select_files(options):
    options.top = os.path.abspath(options.top)
    #print("top",options.top)
    wildcard = WildCard(options.wildcard)
    filepaths = []

    if options.no_walk:
        # Select only the files in the top directory.
        fnames = [os.path.join(options.top, f) for f in os.listdir(options.top)]
        fnames = filter(os.path.isfile, fnames)
        #print(fnames)
        filepaths += wildcard.filter(fnames)

    else:
        # Navigate the directory tree starting from top.
        for root, dirnames, filenames in os.walk(options.top):
            fnames = [os.path.join(root, f) for f in filenames]
            filepaths += wildcard.filter(fnames)

    return filepaths


if __name__ == "__main__":
    sys.exit(main())
