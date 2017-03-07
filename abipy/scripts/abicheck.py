#!/usr/bin/env python
"""
This script checks that the environment on the local machine is properly configured.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import argparse

from monty import termcolor
from monty.termcolor import cprint
from monty.functools import prof_main
from abipy import abilab


@prof_main
def main():

    def str_examples():
        return """\
Usage example:
    abicheck.py
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + abilab.__version__)
    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity.')
    parser.add_argument('--no-colors', default=False, action="store_true", help='Disable ASCII colors.')

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

    if options.no_colors:
        # Disable colors
        termcolor.enable(False)

    errmsg = abilab.abicheck(verbose=options.verbose)
    if errmsg:
        cprint(errmsg, "red")
    else:
        print()
        cprint("Abipy requirements are properly configured", "green")

    return len(errmsg)


if __name__ == "__main__":
    sys.exit(main())
