#!/usr/bin/env python
"""Check that the env on the local machine is properly setup"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import argparse 

from abipy.core.release import __version__

def main():

    def str_examples():
        examples = (
          "\n"
          "Usage example:\n\n" 
          "  abicheck.py \n"
        )
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: 
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)

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

    #try:
    from abipy.abilab import abicheck
    errmsg = abicheck()
    #except:
    #    retcode = 1

    if errmsg:
        print(errmsg)
    else:
        print("Abipy requirements are properly configured")

    return len(errmsg)


if __name__ == "__main__":
    sys.exit(main())
