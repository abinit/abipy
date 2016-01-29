#!/usr/bin/env python
"""
This script retrieve information on a job in the queue
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse

from pymatgen.io.abinit.qjobs import QueueJob
from abipy.core.release import __version__

def main():
    def str_examples():
        examples = """\
Usage example:\n

    abiq.py job_id                => Get info on the job from its job identifier
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                    help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('--qtype', default="pbspro", help="Queue Type.")
    parser.add_argument('job_id', help="Job Indentifier.")

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

    job = QueueJob.from_qtype_and_id(options.qtype, int(options.job_id))
    #print(job)

    print("get_info", job.get_info())
    #print("get_stats", job.get_stats())
    #print("estimated_start_time", job.estimated_start_time())
    print(job)

    return retcode


if __name__ == "__main__":
    sys.exit(main())
