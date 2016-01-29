#!/usr/bin/env python
"""Script to analyze/plot the timing data of single or multiple runs."""
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse 

from monty.termcolor import cprint
from pymatgen.io.abinit.abitimer import AbinitTimerParser
from abipy.core.release import __version__


def main():
    def str_examples():
        return """\
Usage example:\n

  abitime.py.py scan filename(s)  =>
  abitime.py.py walk . --ext=abo  => 
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: 
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parse for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)

    copts_parser.add_argument('-i', '--ipython', default=False, action='store_true',
                       help='Open interactive ipython terminal to interact with the parser.')

    copts_parser.add_argument('-p', '--plot', default=False, action='store_true',
                       help='Plot data with matplotlib.')

    #copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                   help='verbose, can be supplied multiple times to increase verbosity')

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for scan command.
    p_scan = subparsers.add_parser('scan', parents=[copts_parser], help=".")
    p_scan.add_argument('paths', nargs="+", help="File(s) to be analyzed")

    # Subparser for scan command.
    p_walk = subparsers.add_parser('walk', parents=[copts_parser], help=".")

    p_walk.add_argument('top', nargs="?", default=".", help="Top level directory. Defaults to working dir.")
    p_walk.add_argument("-e", '--ext', type=str, default="abo", help="File extension. Default: abo")

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    #import seaborn as sns
    #sns.set(style='ticks', palette='Set2')
    #sns.set(style="dark", palette="Set2")
    #And to remove "chartjunk", do:
    #sns.despine()
    #plt.tight_layout()
    #sns.despine(offset=30, trim=True)

    if options.command == "scan":
        parser = AbinitTimerParser()
        okfiles = parser.parse(options.paths)

        if okfiles != options.paths:
            badfiles = [f for f in options.paths if f not in okfiles]
            cprint("Cannot parse timing data from the following files:", color="magenta")
            for bad in badfiles: print(bad)
            
    elif options.command == "walk":
        print("Walking directory tree from top:", options.top, "Looking for file extension:", options.ext)
        paths = []
        for root, dirs, files in os.walk(options.top):
            for f in files:
                if f.endswith(options.ext):
                    paths.append(os.path.join(root, f))

        if not paths: 
            cprint("Empty file list!", color="magenta")
            return 1

        print("Found %d files\n" % len(paths))
        parser = AbinitTimerParser()
        okfiles = parser.parse(paths)

        if okfiles != paths:
            badfiles = [f for f in paths if f not in okfiles]
            cprint("Cannot parse timing data from the following files:", color="magenta")
            for bad in badfiles: print(bad)

    else:
        raise RuntimeError("Don't know what to do with command: %s!" % options.command)

    if parser is None:
        cprint("Cannot analyze timing data. parser is None", color="magenta")
        return 1

    if options.ipython:
        cprint("Invoking ipython shell. Use parser to access the object inside ipython", color="blue")
        import IPython
        IPython.start_ipython(argv=[], user_ns={"parser": parser})

    if options.plot:
        parser.plot_all()
        #parser.plot_efficiency()
        #parser.plot_pie()
        #parser.plot_stacked_hist()

    print(parser.summarize())
    #table = parser.make_table()
    #print(table)
    return

    for timer in parser.timers():
        print(timer.get_dataframe())
        #print(timer) 
        #print(timer.totable()) 

    return 0


if __name__ == "__main__":
    sys.exit(main())
