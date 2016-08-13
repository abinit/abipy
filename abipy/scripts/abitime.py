#!/usr/bin/env python
"""Script to analyze/plot the timing data of single or multiple runs."""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse
import abipy.abilab as abilab

from monty.termcolor import cprint


def main():
    def str_examples():
        return """\
Usage example:\n

  abitime.py.py scan [FILES] -p     => Parse timing data in files and plot results
  abitime.py.py walk . --ext=abo    => Scan directory tree from `.`, look for files with extension `abo`
                                       parse timind data and plot results.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parse for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-ipy', '--ipython', default=False, action='store_true',
                              help='Open interactive ipython terminal to interact with the parser.')
    copts_parser.add_argument('-nb', '--notebook', default=False, action="store_true", help='Generate jupyter notebook.')
    copts_parser.add_argument('-p', '--plot', default=False, action='store_true', help='Plot data with matplotlib.')
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                              help='verbose, can be supplied multiple times to increase verbosity')

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + abilab.__version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    parser.add_argument('--seaborn', action="store_true", help="Use seaborn settings")

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

    if options.seaborn:
        import seaborn as sns

    if options.command == "scan":
        parser = abilab.AbinitTimerParser()
        okfiles = parser.parse(options.paths)

        if okfiles != options.paths:
            badfiles = [f for f in options.paths if f not in okfiles]
            cprint("Cannot parse timing data from the following files:", color="magenta")
            for bad in badfiles: print(bad)

    elif options.command == "walk":
        print("Walking directory tree from top:", options.top, "Looking for file extension:", options.ext)
        parser, paths, okfiles = abilab.AbinitTimerParser.walk(top=options.top, ext=options.ext)

        if not paths:
            cprint("Empty file list!", color="magenta")
            return 1

        print("Found %d files\n" % len(paths))
        if okfiles != paths:
            badfiles = [f for f in paths if f not in okfiles]
            cprint("Cannot parse timing data from the following files:", color="magenta")
            for bad in badfiles: print(bad)

    else:
        raise RuntimeError("Don't know what to do with command: %s!" % options.command)

    if parser is None:
        cprint("Cannot analyze timing data. parser is None", color="magenta")
        return 1

    print(parser.summarize())

    if options.plot:
        parser.plot_all()

    if options.verbose:
        for timer in parser.timers():
            print(timer.get_dataframe())

    if options.ipython:
        cprint("Invoking ipython shell. Use parser to access the object inside ipython", color="blue")
        import IPython
        IPython.start_ipython(argv=[], user_ns={"parser": parser})
    elif options.notebook:
        parser.make_and_open_notebook(daemonize=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
