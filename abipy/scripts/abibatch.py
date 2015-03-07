#!/usr/bin/env python
"""
This script allows the user to submit the calculations contained in the `Flow`.
It provides both a command line interface as well as a graphical interfaced based on wxpython.
"""
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse
#import time

#from pprint import pprint
#from monty import termcolor
#from monty.termcolor import cprint, get_terminal_size
#from monty.string import make_banner
from pymatgen.io.abinitio.launcher import  BatchLauncher
#import abipy.abilab as abilab


def main():
    def str_examples():
        examples = """\
Usage example:\n

    abibatch.py load    => Keep repeating, stop when no task can be executed.
    abibatch.py load    => Keep repeating, stop when no task can be executed.

    If FLOWDIR is not given, abirun.py automatically selects the database located within 
    the working directory. An Exception is raised if multiple databases are found.

    Options for developers:

        abirun.py prof ABIRUN_ARGS               => to profile abirun.py
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    #parser.add_argument('flowdir', nargs="?", help=("File or directory containing the ABINIT flow"
    #                                                "If not given, the first flow in the current workdir is selected"))

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for submit.
    p_submit = subparsers.add_parser('submit', help="Find all flows in dir and submit them")
    p_submit.add_argument('paths', nargs="*", default=".", help=("Directories containing the object." 
                          "Use current working directory if not specified"))

    # Subparser for load.
    p_load = subparsers.add_parser('load', help="Load object from pickle file.")
    p_load.add_argument('top', help="File or directory containing the object")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    retcode = 0

    if options.command == "submit":
        print("paths", options.paths)

        launcher = BatchLauncher.from_dir(options.paths[0])
        print(launcher.to_string())
        retcode = launcher.submit(verbose=options.verbose)
        print("retcode", retcode)

    elif options.command == "load":
        launcher = BatchLauncher.pickle_load(options.top)
        print(launcher.to_string())

        launcher.show_summary(verbose=options.verbose)
        #launcher.show_status(verbose=options.verbose)

    else:
        raise RuntimeError("Don't know what to do with command %s!" % options.command)

    return retcode


if __name__ == "__main__":
    retcode = 0
    do_prof, do_tracemalloc = 2* [False]
    try:
        do_prof = sys.argv[1] == "prof"
        do_tracemalloc = sys.argv[1] == "tracemalloc"
        if do_prof or do_tracemalloc: sys.argv.pop(1)
    except: 
        pass

    if do_prof:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()

    elif do_tracemalloc:
        # Requires py3.4
        import tracemalloc
        tracemalloc.start()

        retcode = main()

        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')

        print("[Top 10]")
        for stat in top_stats[:10]:
            print(stat)
    else:
        sys.exit(main())

    #open_hook.print_open_files()
    sys.exit(retcode)
