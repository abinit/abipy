#!/usr/bin/env python
"""
This script allows the user to submit the calculations contained in the `AbinitFlow`.
It provides both a command line interface as well as a graphical interfaced based on wxpython.
"""
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse
import time

from pprint import pprint
from monty import termcolor
from pymatgen.io.abinitio.launcher import PyFlowScheduler, PyLauncher
import abipy.abilab as abilab

# Replace python open to detect open files.
#from abipy.tools import open_hook
#open_hook.install()


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def str_examples():
    examples = """
Usage example:\n
    abirun.py [DIRPATH] singleshot              => Fetch the first available task and run it.
    abirun.py [DIRPATH] rapidfire               => Keep repeating, stop when no task can be executed
                                                  due to inter-dependency.
    abirun.py [DIRPATH] gui                      => Open the GUI 
    nohup abirun.py [DIRPATH] sheduler -s 30 &   => Use a scheduler to schedule task submission

    If DIRPATH is not given, abirun.py selects automatically the database located within 
    the working directory. An Exception is raised if multiple databases are found.
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
        flow.show_status()

    if options.command in ["rapid", "rapidfire"]:
        nlaunch = PyLauncher(flow).rapidfire()
        print("Number of tasks launched %d" % nlaunch)
        flow.show_status()

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

        #open_hook.print_open_files()

    if options.command == "status":
        if options.delay:
            print("Entering infinite loop. Press CTRL+C to exit")
            try:
                while True:
                    print(2*"\n" + time.asctime() + "\n")
                    flow.check_status(show=True)
                    if flow.all_ok: break
                    time.sleep(options.delay)
            except KeyboardInterrupt:
                pass
        else:
            flow.show_status(verbose=options.verbose)

    if options.command == "open":
        flow.open_files(what=options.what, wti=None, status=None, op="==")

    if options.command == "cancel":
        print("Number of jobs cancelled %d" % flow.cancel())
        # Remove directory
        if options.rmtree:
            flow.rmtree()

    if options.command == "restart":
        nlaunch, excs = 0, []
        for task in flow.unconverged_tasks:
            try:
                fired = task.restart()
                if fired: nlaunch += 1
            except Exception:
                excs.append(straceback())

        print("Number of jobs restarted %d" % nlaunch)
        if nlaunch:
            # update database
            flow.pickle_dump()

        if excs:
            print("Exceptions raised\n")
            pprint(excs)

    if options.command == "reset":
        count = 0
        for task, wi, ti in flow.iflat_tasks_wti(status=options.task_status):
            task.reset()
            count += 1	
        print("%d tasks have been resetted" % count)

        #if count:
        #    flow.pickle_dump()

        nlaunch = PyLauncher(flow).rapidfire()
        print("Number of tasks launched %d" % nlaunch)
        flow.show_status()

    if options.command == "tail":
        paths = [t.output_file.path for t in flow.iflat_tasks(status="Running")]
        if not paths:
            print("No job is running. Exiting!")
        else:
            print("Press CTRL+C to interrupt. Will follow %d output files" % len(paths))
            try:
                os.system("tail -f %s" % " ".join(paths))
            except KeyboardInterrupt:
                pass

    return retcode


def main():

    # Decorate argparse classes to add portable support for aliases in add_subparsers
    class MyArgumentParser(argparse.ArgumentParser):
        def add_subparsers(self, **kwargs):
            new = super(MyArgumentParser, self).add_subparsers(**kwargs)
            # Use my class
            new.__class__ = MySubParserAction
            return new

    class MySubParserAction(argparse._SubParsersAction):
        def add_parser(self, name, **kwargs):
            """Allows one to pass the aliases option even if this version of ArgumentParser does not support it."""
            try:
                return super(MySubParserAction, self).add_parser(name, **kwargs)
            except Exception as exc:
                if "aliases" in kwargs: 
                    # Remove aliases and try again.
                    kwargs.pop("aliases")
                    return super(MySubParserAction, self).add_parser(name, **kwargs)
                else:
                    # Wrong call.
                    raise exc

    parser = MyArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                        help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('--no-colors', default=False, help='Disable ASCII colors')

    parser.add_argument('--loglevel', default="ERROR", type=str,
                        help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    parser.add_argument('path', nargs="?", help=("File or directory containing the ABINIT flow\n" +
                                                 "If not given, the first flow in the current workdir is selected"))

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_single = subparsers.add_parser('single', aliases=["singleshot"], help="Run single task.")

    # Subparser for rapidfire command.
    p_rapid = subparsers.add_parser('rapid', aliases=["rapidfire"], help="Run all tasks in rapidfire mode")

    # Subparser for scheduler command.
    p_scheduler = subparsers.add_parser('scheduler', aliases=["sched"], help="Run all tasks with a Python scheduler.")

    p_scheduler.add_argument('-w', '--weeks', default=0, type=int, help="number of weeks to wait")

    p_scheduler.add_argument('-d', '--days', default=0, type=int, help="number of days to wait")

    p_scheduler.add_argument('-hs', '--hours', default=0, type=int, help="number of hours to wait")

    p_scheduler.add_argument('-m', '--minutes', default=0, type=int, help="number of minutes to wait")

    p_scheduler.add_argument('-s', '--seconds', default=0, type=int, help="number of seconds to wait")

    # Subparser for status command.
    p_status = subparsers.add_parser('status', help="Show task status.")
    p_status.add_argument('-d', '--delay', default=0, type=int, help=("If 0, exit after the first analysis.\n" + 
                          "If > 0, enter an infinite loop and delay execution for the given number of seconds."))

    # Subparser for cancel command.
    p_cancel = subparsers.add_parser('cancel', help="Cancel the tasks in the queue.")
    p_cancel.add_argument("-r", "--rmtree", action="store_true", default=False, help="Remove flow directory.")

    # Subparser for restart command.
    p_restart = subparsers.add_parser('restart', help="Restart the tasks of the flow that are not converged.")

    # Subparser for restart command.
    p_reset = subparsers.add_parser('reset', help="Reset the tasks of the flow with the specified status.")
    p_reset.add_argument('task_status', default="QCritical") 

    # Subparser for open command.
    p_open = subparsers.add_parser('open', help="Open files in $EDITOR, type `abirun.py ... open --help` for help)")
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
    p_gui.add_argument("--chroot", default="", type=str, help=("Use chroot as new directory of the flow.\n" +
                       "Mainly used for opening a flow located on a remote filesystem mounted with sshfs.\n" +
                       "In this case chroot is the absolute path to the flow on the **localhost**\n",
                       "Note that it is not possible to change the flow from remote when chroot is used."))

    p_new_manager = subparsers.add_parser('new_manager', help="Change the TaskManager.")
    p_new_manager.add_argument("manager_file", default="", type=str, help="YAML file with the new manager")

    p_tail = subparsers.add_parser('tail', help="Use tail to follow the main output file of the flow.")

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

    if options.no_colors:
        termcolor.enable(False)

    # Read the flow from the pickle database.
    if options.path is None:
        # Will try to figure out the location of the Flow.
        options.path = os.getcwd()

    flow = abilab.AbinitFlow.pickle_load(options.path)
    retcode = 0

    if options.command == "gui":
        if options.chroot:
            # Change the workdir of flow.
            print("Will chroot to %s..." % options.chroot)
            flow.chroot(options.chroot)

        from abipy.gui.flowviewer import wxapp_flow_viewer
        wxapp_flow_viewer(flow).MainLoop()

    elif options.command == "new_manager":
        # Read the new manager from file.
        new_manager = abilab.TaskManager.from_file(options.manager_file)

        # Change the manager of the errored tasks.
        status = "S_QCRITICAL"
        #status = "S_ERROR"
        for task, wi, ti in flow.iflat_tasks_wti(status=status):
            task.reset()
            task.set_manager(new_manager)
            
        # Update the database.
        return flow.build_and_pickle_dump()

    else:
        retcode = treat_flow(flow, options)

    return retcode
    

if __name__ == "__main__":
    # perform profiling if `abirun.py prof ...` else run script.
    do_prof = False
    try:
        do_prof = sys.argv[1] == "prof"
        if do_prof: sys.argv.pop(1)
    except: 
        pass

    if not do_prof:
        sys.exit(main())
    else:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()
