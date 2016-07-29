#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse

from abipy import abilab

def main():
    def str_examples():
        return """\
Usage example:\n

  abidiff.py gs_scf run1.abo run2.abo         => Compare the SCF cycles in two output files.
  abidiff.py dfpt2_scf                        => Compare the DFPT SCF cycles in two output files.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('paths', nargs="+", help="List of files to compare")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + abilab.__version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for gs_scf command.
    p_gs_scf = subparsers.add_parser('gs_scf', parents=[copts_parser], help="Compare ground-state SCF cycles.")

    # Subparser for dfpt2_scf command.
    p_dftp2_scf = subparsers.add_parser('dfpt2_scf', parents=[copts_parser], help="Compare DFPT SCF cycles.")

    # Subparser for gsr command.
    #p_gsr = subparsers.add_parser('gsr', parents=[copts_parser], help="Compare electron bands.")

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
    paths = options.paths

    if options.command == "gs_scf":
        f0 = abilab.AbinitOutputFile(paths[0])
        f0.compare_gs_scf_cycles(paths[1:])

    elif options.command == "dfpt2_scf":
        f0 = abilab.AbinitOutputFile(paths[0])
        f0.compare_d2de_scf_cycles(paths[1:])

    #elif options.command == "gsr":

    else:
        raise RuntimeError("Don't know what to do with command: %s!" % options.command)

    return 0


if __name__ == "__main__":
    sys.exit(main())
