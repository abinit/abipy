#!/usr/bin/env python
"""
This script allows the user to submit the calculations contained in the `AbinitFlow`.
It provides both a command line interface as well as a graphical interfaced based on wxpython.
"""
from __future__ import division, print_function

import sys 
import os
import warnings
import argparse
import abipy.abilab as abilab

from pymatgen.io.abinitio.launcher import PyResourceManager, PyLauncher
from abipy.tools import pprint_table, StringColorizer

def str_examples():
    examples = """
Usage example:\n
    abirun.py DIRECTORY singleshot  => Fetch the first available task and run it.
    abirun.py DIRECTORY rapidfire   => Keep repeating, stop when no task can be executed
                                       due to inter-dependency.
    abirun.py DIRECTORY gui         => Open the GUI 
"""
    return examples

def show_examples_and_exit(err_msg=None, error_code=1):
    """Display the usage of the script."""
    sys.stderr.write(str_examples())
    if err_msg: 
        sys.stderr.write("Fatal Error\n" + err_msg + "\n")

    sys.exit(error_code)


def treat_flow(flow, options):
    retcode = 0

    # Dispatch.
    if options.command in ["single", "singleshot"]:
        nlaunch = PyLauncher(flow).single_shot()
        print("Number of tasks launched %d" % nlaunch)

    if options.command in ["rapid", "rapidfire"]:
        nlaunch = PyLauncher(flow).rapidfire()
        print("Number of tasks launched %d" % nlaunch)

    #if options.command == "pymanager":
    #    retcodes = PyResourceManager(work, max_ncpus=1, sleep_time=5).run()
    #    recode = max(retcodes)
    #    print("all tasks completed with return code %s" % retcode)
    #    work.pickle_dump()

    if options.command == "status":
        colorizer = StringColorizer(stream=sys.stdout)

        for i, work in enumerate(flow):
            print(80*"=")
            print("Workflow #%d: %s, Finalized=%s\n" % (i, work, work.finalized) )

            table = [["Task", "Status", "queue_id", "Errors", "Warnings", "Comments", "MPI", "OMP"]]
            for task in work:
                task_name = os.path.basename(task.name)

                # Parse the events in the main output.
                report = task.get_event_report()

                events = map(str, 3*["N/A"])
                if report is not None: 
                    events = map(str, [report.num_errors, report.num_warnings, report.num_comments])

                colour = {
                    #task.S_READY: "Ready",
                    #task.S_SUB: "Submitted",
                    #task.S_RUN: "Running",
                    #task.S_DONE: "Done",
                    task.S_ERROR: "red",
                    task.S_OK: "blue",
                }.get(task.status, None)
                #task_name = colorizer(task_name, colour)

                cpu_info = map(str, [task.mpi_ncpus, task.omp_ncpus])

                table.append([task_name, str(task.status), str(task.queue_id)] + events + cpu_info)

            pprint_table(table)

    return retcode


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('--loglevel', default="ERROR", type=str,
                         help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('paths', nargs="+", help="Directories containing ABINIT workflows")

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_single = subparsers.add_parser('singleshot', help="Run single task.") #, aliases=["single"])

    # Subparser for rapidfire command.
    p_rapid = subparsers.add_parser('rapidfire', help="Run all tasks in rapidfire mode") # aliases=["rapid"])

    # Subparser for pymanager command.
    p_pymanager = subparsers.add_parser('pymanager', help="Run all tasks.")

    # Subparser for status command.
    p_status = subparsers.add_parser('status', help="Show task status.")

    # Subparser for gui command.
    p_gui = subparsers.add_parser('gui', help="Open GUI.")

    # Parse command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    if options.verbose:
        print(options)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    paths = options.paths
    print(paths)

    # Walk through each directory in options.paths and find the pickle databases.
    paths = []
    for root in options.paths:
        #print(root)
        for dirpath, dirnames, filenames in os.walk(root):
            for fname in filenames:

                if fname == abilab.AbinitFlow.PICKLE_FNAME:
                    paths.append(os.path.join(dirpath, fname))

    import cPickle as pickle
    paths = [paths[0]]
    #print("paths", str(paths))

    options.paths = paths
    retcode = 0

    if len(options.paths) == 0:
        warnings.warn("The directories specifies do not contain any valid AbinitFlow")
        return 1

    if len(options.paths) > 1:
        raise ValueError("Multiple flows are not supported")

    path = options.paths[0]

    # Read the worflow from the pickle database.
    with open(path, "rb") as fh:
        flow = pickle.load(fh)

    flow.connect_signals()
    #for w in flow:
    #    print(w)
    #flow.show_dependencies()

    # Recompute the status of each task since tasks that
    # have been submitted previously might be completed.
    flow.check_status()

    if options.command == "gui":
        from abipy.gui.workflow_viewer import wxapp_flow_viewer
        wxapp_flow_viewer(flow).MainLoop()

    else:
        retcode = treat_flow(flow, options)

    return retcode
    

if __name__ == "__main__":
    sys.exit(main())
