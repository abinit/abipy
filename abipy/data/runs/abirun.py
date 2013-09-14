#!/usr/bin/env python
"""
This script ...
"""
from __future__ import division, print_function

import sys 
import os
import argparse
import abipy.abilab as abilab

from pymatgen.io.abinitio.launcher import PyResourceManager
from abipy.tools import pprint_table, StringColorizer

def str_examples():
    examples = """
Usage example:\n
    abirun.py singleshot  => Fetch the first available task and run it.
    abirun.py rapidfire   => Keep repeating, stop when no task can be executed
                             due to workflow dependency.
    abirun.py gui         => Open the GUI 
"""
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")

    sys.exit(error_code)

def treat_workflow(work, options):

    # Read the worflow from the pickle database.
    retcode = 0

    # Recompute the status of each task since tasks that
    # have been submitted previously might be completed.
    work.recheck_status()

    # Dispatch.
    if options.command in ["single", "singleshot"]:
        try:
            task = work.fetch_task_to_run()

            if task is None:
                print("No task to run!. Possible deadlock")

            else:
                print("got task", task)
                task.start()

                #task.start_and_wait()
                #retcode = task.returncode
                #print("Task returncode", task.returncode)

        except StopIteration as exc:
            print(str(exc))

        finally:
            work.pickle_dump()

    if options.command == "pymanager":
        retcodes = PyResourceManager(work, max_ncpus=1, sleep_time=5).run()
        recode = max(retcodes)
        print("all tasks completed with return code %s" % retcode)

        work.pickle_dump()

    if options.command in ["rapid", "rapidfire"]:
        while True:
            try:
                task = work.fetch_task_to_run()
                if task is None:
                    break
                print("got task",task)
                task.start()

            except StopIteration as exc:
                print(str(exc))
                break

        work.pickle_dump()

    if options.command == "status":
        print(work)
        colorizer = StringColorizer(stream=sys.stdout)

        table = [["Task", "Status", "queue_id", "Errors", "Warnings", "Comments"]]
        for task in work:
            task_name = os.path.basename(task.name)

            # Parse the events in the main output.
            report = task.get_event_report()

            events = map(str, 3*["N/A"])
            if report is not None: 
                events = map(str, [report.num_errors, report.num_warnings, report.num_comments])

            str_status = task.str_status

            colour = {
                #task.S_READY: "Ready",
                #task.S_SUB: "Submitted",
                #task.S_RUN: "Running",
                #task.S_DONE: "Done",
                task.S_ERROR: "red",
                task.S_OK: "blue",
            }.get(task.status, None)
            #task_name = colorizer(task_name, colour)

            table.append([task_name, str_status, str(task.queue_id)] + events)

        pprint_table(table)

    return retcode


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_single = subparsers.add_parser('singleshot', help="Run single task.") #, aliases=["single"])
    p_single.add_argument('paths', nargs="+", help="Directories containing ABINIT workflows")

    # Subparser for rapidfire command.
    p_rapid = subparsers.add_parser('rapidfire', help="Run all tasks in rapidfire mode") # aliases=["rapid"])
    p_rapid.add_argument('paths', nargs="+", help="Directories containing ABINIT workflows")

    # Subparser for pymanager command.
    p_pymanager = subparsers.add_parser('pymanager', help="Run all tasks.")
    p_pymanager.add_argument('paths', nargs="+", help="Directories containing ABINIT workflows")

    # Subparser for status command.
    p_status = subparsers.add_parser('status', help="Show task status.")
    p_status.add_argument('paths', nargs="+", help="Directories containing ABINIT workflows")

    # Subparser for gui command.
    p_gui = subparsers.add_parser('gui', help="Open GUI.")
    p_gui.add_argument('paths', nargs="+", help="Directories containing ABINIT workflows")

    # Parse command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    if options.verbose:
        print(options)

    paths = options.paths

    # Walk through each directory in options.paths and find the pickle databases.
    paths = []
    for root in options.paths:

        for dirpath, dirnames, filenames in os.walk(root):
            for fname in filenames:
                if fname == abilab.Workflow.PICKLE_FNAME:
                    paths.append(os.path.join(dirpath, fname))

    options.paths = paths
    retcode = 0

    workflows = [abilab.Workflow.pickle_load(path) for path in options.paths]

    if options.command == "gui":
        from abipy.gui.workflow_viewer import wxapp_workflow_viewer

        #for work in workflows:
        #    work.recheck_status()

        wxapp_workflow_viewer(workflows).MainLoop()

    else:
        for work in workflows:
            retcode = treat_workflow(work, options)
            if retcode != 0:
                return retcode

    return retcode
    

if __name__ == "__main__":
    sys.exit(main())
