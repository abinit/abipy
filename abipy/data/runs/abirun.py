#!/usr/bin/env python
"""
This script ...
"""
from __future__ import division, print_function

import sys 
import os
import argparse
import abipy.abilab as abilab

from pymatgen.io.abinitio.launcher import SimpleResourceManager
from abipy.tools import pprint_table

def str_examples():
    examples = """
Usage example:\n
    abirun.py single      => Run single task.
    abirun.py gui         => Open GUI
"""
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")

    sys.exit(error_code)

def treat_workflow(path, options):

    # Read the worflow from the pickle database.
    retcode = 0
    work = abilab.Workflow.pickle_load(path)

    # Recompute the status of each task since tasks that
    # have been submitted previously might be completed.
    work.recheck_status()

    if options.verbose:
        print(work)
        print([task.status for task in work])

    if options.command == "single":
        try:
            task = work.fetch_task_to_run()

            if task is None:
                print("No task to run!. Possible deadlock")

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

    if options.command == "all":
        retcodes = SimpleResourceManager(work, max_ncpus=1, sleep_time=5).run()
        recode = max(retcodes)
        print("all tasks completed with return code %s" % retcode)

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

    if options.command == "status":
        print(work)

        table = [["Task", "Status", "Errors", "Warnings", "Comments"]]
        for task in work:
            task_name = os.path.basename(task.name)

            # Parse the events in the main output.
            report = task.parse_events()
            events = map(str, 3*[0])
            if report is not None: 
                events = map(str, [report.num_errors, report.num_warnings, report.num_comments])

            table.append([task_name, task.str_status] + events)

        pprint_table(table)

    if options.command == "gui":
        # TODO 
        # Modify the viewer so that we can have multiple tabs for each workflow.
        from abipy.gui.workflow_viewer import wxapp_workflow_viewer
        wxapp_workflow_viewer(work).MainLoop()

    return retcode

def main():
    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Subparser for single command.
    p_single = subparsers.add_parser('single', help="Run single task.")
    p_single.add_argument('path', nargs="+", help='Directory with __workflow__.pickle database or file.')

    # Subparser for all command.
    p_all = subparsers.add_parser('all', help="Run all tasks.")
    p_all.add_argument('path', nargs="+", help='Directory with __workflow__.pickle database or file.')

    # Subparser for status command.
    p_status = subparsers.add_parser('status', help="Show task status.")
    p_status.add_argument('path', nargs="+", help='Directory with __workflow__.pickle database or file.')

    # Subparser for gui command.
    p_gui = subparsers.add_parser('gui', help="Open GUI.")
    p_gui.add_argument('path', nargs="+", help='Directory with __workflow__.pickle database or file.')

    ###################################
    # Parse command line and dispatch #
    ###################################
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    if options.verbose:
        print(options)

    retcode = 0
    for path in options.path:
        retcode = treat_workflow(path, options)
        if retcode !=0:
            return retcode

    return retcode
    

if __name__ == "__main__":
    sys.exit(main())
