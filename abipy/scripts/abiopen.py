#!/usr/bin/env python
from __future__ import print_function, division

import sys
import os
import argparse 

import abipy.gui.wxapps as wxapps 


def str_examples():
    examples = (
      "\n"
      "Usage example:\n\n" 
      "abiopen.py  dirpath    ==> Visualize all netcdf files (*.nc) in the directory dirpath.\n"
      "abiopen.py -f *WFK.nc  ==> Use a glob pattern to match all the WFK files in the current directory.\n"
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

    parser.add_argument('-m', '--view-mode', type=str, default="list",
                        help="View mode (list or tree)")

    parser.add_argument('-w', '--wildcard', type=str, default="*.nc",
                        help="wildcards. Default *.nc")

    parser.add_argument('-f', '--filepaths', nargs="+", default=None,
                        help="List of files.")

    #parser.add_argument('-r', '--recurse', type=bool, default=False,
    #                    help="Recursive mode.")

    parser.add_argument("dirpaths", nargs="*", help="List of directories.")

    # Parse the command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    if not options.dirpaths and not options.filepaths:
        options.dirpaths = [os.getcwd()]

    if options.view_mode in ["list", "l"]:
        app = wxapps.wxapp_listbrowser(dirpaths=options.dirpaths, 
                                       filepaths=options.filepaths, 
                                       wildcard=options.wildcard)

    elif options.view_mode in ["compare", "c"]:
        app = wxapps.wxapp_comparison(dirpaths=options.dirpaths, 
                                      filepaths=options.filepaths, 
                                      wildcard=options.wildcard)

    elif options.view_mode in ["tree", "t"]:
        app = wxapps.wxapp_dirbrowser(dirpath=options.dirpaths[0])

    else:
        raise ValueError("Wrong view_mode %s" % options.view_mode)

    app.MainLoop()

    return 0

if __name__ == "__main__":
    sys.exit(main())
