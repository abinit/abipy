#!/usr/bin/env python
"""
This script ...
"""
from __future__ import division, print_function

import sys 
import os
import argparse

import abipy.abilab as abilab
#from abipy.tools.devtools import number_of_cpus

def str_examples():
    examples = """
      Usage example:\n\n
      abirun.py single      => Run single task.
      abirun.py rapidfire   => Run tasks in rapidfire mode.
      abirun.py gui         => Open GUI 
    """
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")

    sys.exit(error_code)


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for single command.
    p_single = subparsers.add_parser('single', help="Run single task.")

    p_single.add_argument('path', help='Directory with __workflow__.pickle database or file.')
    #p_single.add_argument('-xc', '--xc-type', default="GGA", help="XC functional type. (default GGA).")  

    # Subparser for gui command.
    p_gui = subparsers.add_parser('gui', help="Open GUI.")

    p_gui.add_argument('path', help='Directory with __workflow__.pickle database or file.')

    ###################################
    # Parse command line and dispatch #
    ###################################
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    if options.verbose:
        print(options)

    # Read the worflow for the pickle database.
    work = abilab.Workflow.pickle_load(options.path)

    # Recompute the status of each task since tasks that
    # have been submitted previously might be completed.
    #work.recheck_task_status()

    if options.verbose:
        print(work)
        print([task.status for task in work])

    if options.command == "single":
        try:
            task = work.fetch_task_to_run()

            if task is None:
                print("No task to run!")

            else:
                print("got task",task)
                task.start()

                #task.start_and_wait()
                #retcode = task.returncode
                #print("Task returncode", task.returncode)

        except StopIteration as exc:
            print(str(exc))

        finally:
            work.pickle_dump()

    #if options.command == "rapidfire":
    #    while True:
    #        try:
    #            task = work.fetch_task_to_run()

    #            if task is None:
    #                break

    #            print("got task",task)
    #            task.start()

    #        except StopIteration as exc:
    #            print(str(exc))
    #            break

    #   work.pickle_dump()

    if options.command == "gui":
        # TODO
        #wxapps.workflow_viewer(work).MainLoop()
        #work.wxbrowse()
        work.wxshow_outputs()

    return 0

if __name__ == "__main__":
    sys.exit(main())
