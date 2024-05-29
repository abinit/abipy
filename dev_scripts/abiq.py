#!/usr/bin/env python
"""
This script retrieve information on a job in the queue
"""
import sys
import os
import argparse
import abipy.tools.cli_parsers as cli
import abipy.flowtk.qutils as qu

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

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    cli.set_loglevel(options.loglevel)

    jobs_dict = qu.slurm_get_jobs()
    for job_id, dct in jobs_dict.items():
        print(f"{job_id=}", dct)

    return 0


if __name__ == "__main__":
    sys.exit(main())
