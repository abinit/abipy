#!/usr/bin/env python
"""
This script allows the user to submit the calculations contained in the `AbinitFlow`.
It provides both a command line interface as well as a graphical interfaced based on wxpython.
"""
from __future__ import division, print_function

import sys
import os
import argparse
import abipy.abilab as abilab

from pymatgen.io.abinitio.launcher import PyFlowScheduler, PyLauncher


def str_examples():
    examples = """
Usage example:\n
    abirun.py DIRPATH singleshot              => Fetch the first available task and run it.
    abirun.py DIRPATH rapidfire               => Keep repeating, stop when no task can be executed
                                                  due to inter-dependency.
    abirun.py DIRPATH gui                      => Open the GUI 
    nohup abirun.py DIRPATH sheduler -s 30 &   => Use a scheduler to schedule task submission
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

    if options.command == "scheduler":

        opt_names = [
            "weeks",
            "days",
            "hours",
            "minutes",
            "seconds",
        ]

        sched_options = {oname: getattr(options, oname) for oname in opt_names}

        if all(v == 0 for v in sched_options.values()):
            sched = PyFlowScheduler.from_user_config()
        else:
            sched = PyFlowScheduler(**sched_options)

        # Check that the env on the local machine is properly setup 
        # before starting the scheduler.
        abilab.abicheck()

        sched.add_flow(flow)
        print(sched)
        sched.start()

    if options.command == "status":
        flow.show_status()
        #import pstats, cProfile
        #cProfile.runctx("flow.show_status()", globals(), locals(), "Profile.prof")
        #s = pstats.Stats("Profile.prof")
        #s.strip_dirs().sort_stats("time").print_stats()

    if options.command == "open":
        flow.open_files(what=options.what, wti=None, status=None, op="==")

    if options.command == "cancel":
        num_cancelled = flow.cancel()
        print("Number of jobs cancelled %d" % num_cancelled)

    return retcode


def main():
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('path', nargs="?", help=("File or directory containing the ABINIT flow\n" +
                                                 "If not given, the first flow in the current workdir is selected"))

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_single = subparsers.add_parser('singleshot', help="Run single task.") #, aliases=["single"])

    # Subparser for rapidfire command.
    p_rapid = subparsers.add_parser('rapidfire', help="Run all tasks in rapidfire mode") # aliases=["rapid"])

    # Subparser for scheduler command.
    p_scheduler = subparsers.add_parser('scheduler', help="Run all tasks with a Python scheduler.")

    p_scheduler.add_argument('-w', '--weeks', default=0, type=int, help="number of weeks to wait")

    p_scheduler.add_argument('-d', '--days', default=0, type=int, help="number of days to wait")

    p_scheduler.add_argument('-hs', '--hours', default=0, type=int, help="number of hours to wait")

    p_scheduler.add_argument('-m', '--minutes', default=0, type=int, help="number of minutes to wait")

    p_scheduler.add_argument('-s', '--seconds', default=0, type=int, help="number of seconds to wait")

    # Subparser for status command.
    p_status = subparsers.add_parser('status', help="Show task status.")

    # Subparser for scheduler command.
    p_cancel = subparsers.add_parser('cancel', help="Cancel the tasks in the queue.")

    # Subparser for open command.
    p_open = subparsers.add_parser('open', help="Open files (command line interface)")

    p_open.add_argument('what', default="o", 
        help="""\
Specify the files to open. Possible choices:\n
i ==> input_file\n
o ==> output_file\n
f ==> files_file\n              
j ==> job_file\n                
l ==> log_file\n                
e ==> stderr_file\n             
q ==> qerr_file\n
""")

    # Subparser for gui command.
    p_gui = subparsers.add_parser('gui', help="Open GUI.")
    p_gui.add_argument("--chroot", default="", type=str, help="Directory for chroot.")

    p_new_manager = subparsers.add_parser('new_manager', help="Change the TaskManager.")
    p_new_manager.add_argument("manager_file", default="", type=str, help="YAML file with the new manager")

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc: 
        show_examples_and_exit(error_code=1)

    if options.verbose:
        print("options", options)

    # loglevel is bound to the string value obtained from the command line argument. 
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    # Read the flow from the pickle database.

    if options.path is None:
        # "Will try to figure out the location of the Flow"
        options.path = os.getcwd()

    flow = abilab.AbinitFlow.pickle_load(options.path)
    retcode = 0

    if options.command == "gui":
        if options.chroot:
            # Change the workdir of flow.
            print("Will chroot to %s" % options.chroot)
            flow.chroot(options.chroot)

        from abipy.gui.flowviewer import wxapp_flow_viewer
        wxapp_flow_viewer(flow).MainLoop()

    elif options.command == "new_manager":
        # Read the new manager from file.
        new_manager = abilab.TaskManager.from_file(options.manager_file)

        # Change the manager of the errored tasks.
        for task in flow.iflat_tasks(status="S_ERROR"):
            task.reset()
            task.set_manager(new_manager)
            
        # Update the database.
        return flow.build_and_pickle_dump()

    else:
        retcode = treat_flow(flow, options)

    return retcode
    

if __name__ == "__main__":
    sys.exit(main())
    #import pstats, cProfile
    #cProfile.runctx("main()", globals(), locals(), "Profile.prof")
    #s = pstats.Stats("Profile.prof")
    #s.strip_dirs().sort_stats("time").print_stats()
