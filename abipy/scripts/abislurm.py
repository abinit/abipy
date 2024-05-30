#!/usr/bin/env python
"""
This script retrieve information on Slurm jobs.
"""
import sys
import os
import argparse
import abipy.tools.cli_parsers as cli
import abipy.flowtk.qutils as qu

from abipy.core.release import __version__


def get_epilog() -> str:
   return """\
Usage example:\n

    abislurm.py running                => Get info on all the running jobs
    abislurm.py completed 111 112      => Get info on completed jobs
"""


def get_parser(with_epilog=False):
    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='verbose, can be supplied multiple times to increase verbosity')

    #parser.add_argument('-V', '--version', action='version', version="%(prog)s version " + __version__)
    #parser.add_argument('--loglevel', default="ERROR", type=str,
    #                    help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity')
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    job_ids_parser = argparse.ArgumentParser(add_help=False)
    job_ids_parser.add_argument('job_ids', nargs="+", help="List of job ids.")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help',
        description="Valid subcommands, use command --help for help")

    # Subparser for running command.
    p_running = subparsers.add_parser('running', parents=[copts_parser],
        help="Check info on all the running jobs.")

    # Subparser for completed command.
    p_completed = subparsers.add_parser('completed', parents=[copts_parser, job_ids_parser],
        help="Returning info on completed jobs.")

    return parser


def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = get_parser(with_epilog=True)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    if not options.command:
        show_examples_and_exit(error_code=1)

    cli.set_loglevel(options.loglevel)

    if options.verbose > 2:
        print(options)

    if options.command == "running":
        jobs_dict = qu.slurm_get_jobs()
        for job_id, dct in jobs_dict.items():
            print(f"{job_id=}", dct)

    #elif options.command == "running_from_abilogs":

    elif options.command == "completed":
        for job_id in options.job_ids:
            print(qu.get_completed_job_info(job_id))

    #elif options.command == "completed_from_abilogs":
    #    job_ids = []

    else:
        raise ValueError(f"Unsupported {options.command=}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
