#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse 

from abipy import abilab

import logging
logger = logging.getLogger(__name__)


def main():

    def str_examples():
        examples = (
          "\n"
          "Usage example:\n\n" 
          "  abiopen.py out_DDB\n"
          "  abiopen.py out_GSR\n"
          "\nMany other Abinit files are supported. Just try!\n"
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

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                     help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument("filepath", help="File to open.")

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

    # Start ipython shell with namespace 
    abifile = abilab.abiopen(options.filepath)
    import IPython
    # USe embed because I don't know how to show a header with start_ipython.
    IPython.embed(header="The Abinit file is bound to the `abifile` variable.\nTry `print(abifile)`")
    #IPython.start_ipython(argv=options.argv, 
    #                      user_ns={"abifile": abifile},
    #                      banner="hello",
    #                      banner1="hello1",
    #                      header="hello_header",
    #                      )
    # 
    return 0


if __name__ == "__main__":
    sys.exit(main())
