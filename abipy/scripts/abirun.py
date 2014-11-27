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
from monty.termcolor import cprint
from pymatgen.io.abinitio.launcher import PyFlowScheduler, PyLauncher
import abipy.abilab as abilab

# Replace python open to detect open files.
#from abipy.tools import open_hook
#open_hook.install()


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def as_slice(obj):
    """
    Convert an integer, a string or a slice object into slice.

    >>> assert as_slice(5) == slice(5, 6, 1)
    >>> assert as_slice("[1:4]") == slice(1, 4, 1)
    >>> assert as_slice("1::2") == slice(1, None, 2)
    """
    if isinstance(obj, slice): return obj

    try:
        # integer.
        if int(obj) == float(obj): return slice(int(obj), int(obj)+1)
    except:
        # assume string defining a python slice [start:stop:step]
        if not obj: return None
        if obj.count("[") + obj.count("]") not in (0, 2): 
            raise ValueError("Invalid string %s" % obj)

        obj = obj.replace("[", "").replace("]", "")
        n = obj.count(":") 
        if n == 0:
            obj = int(obj)
            return slice(obj, obj+1)

        tokens = [int(f) if f else None for f in obj.split(":")]
        if len(tokens) == 2: tokens.append(1)
        if tokens[2] is None: tokens[2] = 1

        return slice(*tokens)

    raise ValueError("Cannot convert %s into a slice:\n%s" % (type(obj), obj))


def main():

    def str_examples():
        examples = """\
Usage example:\n
    abirun.py [DIRPATH] single                   => Fetch the first available task and run it.
    abirun.py [DIRPATH] rapid                    => Keep repeating, stop when no task can be executed
                                                    due to inter-dependency.
    abirun.py [DIRPATH] gui                      => Open the GUI 
    nohup abirun.py [DIRPATH] sheduler -s 30 &   => Use a scheduler to schedule task submission

    If DIRPATH is not given, abirun.py automatically selects the database located within 
    the working directory. An Exception is raised if multiple databases are found.

    Options for developers:
        abirun.py prof ABIRUN_OPTIONS      to profile abirun.py
        abirun.py tracemalloc ABIRUN_ARGS  to trace memory blocks allocated by Python
"""
        return examples


    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    # Build the parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

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
    p_single = subparsers.add_parser('single', help="Run single task.")

    # Subparser for rapidfire command.
    p_rapid = subparsers.add_parser('rapid', help="Run all tasks in rapidfire mode")

    # Subparser for scheduler command.
    p_scheduler = subparsers.add_parser('scheduler', help="Run all tasks with a Python scheduler.")

    p_scheduler.add_argument('-w', '--weeks', default=0, type=int, help="number of weeks to wait")

    p_scheduler.add_argument('-d', '--days', default=0, type=int, help="number of days to wait")

    p_scheduler.add_argument('-hs', '--hours', default=0, type=int, help="number of hours to wait")

    p_scheduler.add_argument('-m', '--minutes', default=0, type=int, help="number of minutes to wait")

    p_scheduler.add_argument('-s', '--seconds', default=0, type=int, help="number of seconds to wait")

    # Subparser for status command.
    p_status = subparsers.add_parser('status', help="Show task status.")
    p_status.add_argument('-d', '--delay', default=0, type=int, help=("If 0, exit after the first analysis.\n" + 
                          "If > 0, enter an infinite loop and delay execution for the given number of seconds."))
    p_status.add_argument("-w", '--work-slice', default="", type=str, 
                          help=("Select the list of works to analyze (python syntax for slices):\n"
                          "-w1 to select the second workflow, -w:3 for 0,1,2, -w-1 for the last workflow, -w::2 for even indices"))

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
Specify the files to open. Possible choices:
    i ==> input_file
    o ==> output_file
    f ==> files_file
    j ==> job_file
    l ==> log_file
    e ==> stderr_file
    q ==> qerr_file
""")
    p_open.add_argument('--wti', default=None, help="Index of workflow:task")
    p_open.add_argument('--nids', default=None, help="Node identifier(s) used to select the task")

    # Subparser for gui command.
    p_gui = subparsers.add_parser('gui', help="Open GUI.")
    p_gui.add_argument("--chroot", default="", type=str, help=("Use chroot as new directory of the flow.\n" +
                       "Mainly used for opening a flow located on a remote filesystem mounted with sshfs.\n" +
                       "In this case chroot is the absolute path to the flow on the **localhost**\n",
                       "Note that it is not possible to change the flow from remote when chroot is used."))

    p_new_manager = subparsers.add_parser('new_manager', help="Change the TaskManager.")
    p_new_manager.add_argument("manager_file", default="", type=str, help="YAML file with the new manager")

    p_tail = subparsers.add_parser('tail', help="Use tail to follow the main output file of the flow.")
    p_tail.add_argument('what_tail', nargs="?", type=str, default="o", help="What to follow: o for output (default), l for logfile, e for stderr")

    p_qstat = subparsers.add_parser('qstat', help="Show additional info on the jobs in the queue.")
    #p_qstat.add_argument('what_tail', nargs="?", type=str, default="o", help="What to follow: o for output (default), l for logfile, e for stderr")

    p_deps = subparsers.add_parser('deps', help="Show dependencies.")

    p_robot = subparsers.add_parser('robot', help="Use a robot to analyze the results of multiple tasks (requires ipython)")
    p_robot.add_argument('robot_ext', nargs="?", type=str, default="GSR", help="The file extension of the netcdf file")

    p_inspect = subparsers.add_parser('inspect', help="Inspect the tasks")
    #p_inspect.add_argument('--wti', default=None, help="Index of workflow:task")
    p_inspect.add_argument('--nids', default=None, help="Node identifier(s) used to select the task")

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

    elif options.command in ("single", "singleshot"):
        nlaunch = PyLauncher(flow).single_shot()
        flow.show_status()
        print("Number of tasks launched: %d" % nlaunch)

    elif options.command in ("rapid", "rapidfire"):
        nlaunch = PyLauncher(flow).rapidfire()
        flow.show_status()
        print("Number of tasks launched: %d" % nlaunch)

    elif options.command == "scheduler":

        sched_options = {oname: getattr(options, oname) for oname in 
            ("weeks", "days", "hours", "minutes", "seconds")}

        if all(v == 0 for v in sched_options.values()):
            sched = PyFlowScheduler.from_user_config()
        else:
            sched = PyFlowScheduler(**sched_options)

        # Check that the env on the local machine is properly setup 
        # before starting the scheduler.
        abilab.abicheck()
        #errors = flow.loop_before_you_leap()
        #if errors:
        #    raise RuntimeError("look_before_you_leap returned:\n %s" % str(errors))

        sched.add_flow(flow)
        print(sched)
        sched.start()

        #open_hook.print_open_files()

    elif options.command == "status":
        work_slice = as_slice(options.work_slice)
        #print(work_slice, options.work_slice)

        if options.delay:
            cprint("Entering infinite loop. Press CTRL+C to exit", color="magenta", end="", flush=True)
            try:
                while True:
                    print(2*"\n" + time.asctime() + "\n")
                    flow.check_status()
                    flow.show_status(verbose=options.verbose, work_slice=work_slice)
                    if flow.all_ok: break
                    time.sleep(options.delay)
            except KeyboardInterrupt:
                pass
        else:
            flow.show_status(verbose=options.verbose, work_slice=work_slice)
            if flow.manager.has_queue:
                print("Total number of jobs in queue: %s" % flow.manager.get_njobs_in_queue())

    elif options.command == "open":

        if options.nids is not None:
            nids = map(int, options.nids.split(","))
            wti = flow.wti_from_nids(nids)

        elif options.wti is not None:
            wti = [int(t) for t in options.wti.split(":")]
        else:
            raise ValueError("Either nids or wti option must be specified")

        flow.open_files(what=options.what, wti=wti, status=None, op="==")

    elif options.command == "cancel":
        print("Number of jobs cancelled %d" % flow.cancel())
        # Remove directory
        if options.rmtree: flow.rmtree()

    elif options.command == "restart":
        nlaunch, excs = 0, []
        for task in flow.unconverged_tasks:
            try:
                fired = task.restart()
                if fired: nlaunch += 1
            except Exception:
                excs.append(straceback())

        cprint("Number of jobs restarted %d" % nlaunch, "blue")
        if nlaunch:
            # update database
            flow.pickle_dump()

        if excs:
            print("Exceptions raised\n")
            pprint(excs)

    elif options.command == "reset":
        print("Will reset tasks with status: %s" % options.task_status)

        count = 0
        for task, wi, ti in flow.iflat_tasks_wti(status=options.task_status):
            task.reset()
            count += 1	
        cprint("%d tasks have been resetted" % count, "blue")

        nlaunch = PyLauncher(flow).rapidfire()
        flow.show_status()
        print("Number of tasks launched: %d" % nlaunch)

        if nlauch == 0:
            deadlocked, runnables, running = flow.deadlocked_runnables_running()
            print("deadlocked:", deadlocked)
            print("runnables:", runnables)
            print("running:", running)
            if deadlocked and not (runnables or running):
                print("*** Flow is deadlocked ***")

        flow.pickle_dump()

    elif options.command == "tail":
        def get_path(task):
            """Helper function used to select the files of a task."""
            choices = {
                "o": task.output_file,
                "l": task.log_file,
                "e": task.stderr_file,
            }
            return getattr(choices[options.what_tail], "path")

        paths = [get_path(task )for task in flow.iflat_tasks(status="Running")]
        if not paths:
            cprint("No job is running. Exiting!", "red")
        else:
            cprint("Press CTRL+C to interrupt. Number of output files %d" % len(paths), color="magenta", end="", flush=True)
            try:
                os.system("tail -f %s" % " ".join(paths))
            except KeyboardInterrupt:
                pass

    elif options.command == "qstat":
        for task in flow.iflat_tasks(): #status=
            if not task.qjob: continue
            print("qjob", task.qjob)
            print("info", task.qjob.get_info())
            print("e start-time", task.qjob.estimated_start_time())
            print("qstats", task.qjob.get_stats())

    elif options.command == "deps":
        flow.check_status()
        flow.show_dependencies()

    elif options.command == "robot":
        import IPython
        with abilab.abirobot(flow, options.robot_ext) as robot:
            IPython.embed(
                header=str(robot) + "\nType `robot` in the terminal and use <TAB> to list its methods", 
                robot=robot
            )

    elif options.command == "inspect":
        if options.nids is not None:
            nids = map(int, options.nids.split(","))
            tasks = flow.tasks_from_nids(nids)
        else:
            tasks = list(flow.iflat_tasks())

        # Use different thread to inspect the task so that master can catch KeyboardInterrupt and exit.
        # One could use matplotlib non-blocking interface with show(block=False) but this one seems to work well.
        from multiprocessing import Process

        def plot_graphs():
            for task in tasks:
                if hasattr(task, "inspect"):
                    task.inspect()
                else:
                    cprint("Task %s does not provide an inspect method" % task, color="blue")

        p = Process(target=plot_graphs)
        p.start()
        num_tasks = len(tasks)

        if num_tasks == 1:
            p.join()
        else:
            cprint("Will produce %d matplotlib plots. Press CTRL+C to interrupt..." % num_tasks, color="magenta", end="", flush=True)
            try:
                p.join()
            except KeyboardInterrupt:
                print("\nTerminating thread...")
                p.terminate()

    else:
        raise RuntimeError("Don't know what to do with command %s!" % options.command)

    return retcode
    

if __name__ == "__main__":
    # profile if `abirun.py prof ...` else run script.
    do_prof = False
    try:
        do_prof = sys.argv[1] == "prof"
        do_tracemalloc = sys.argv[1] == "tracemalloc"
        if do_prof or do_tracemalloc: sys.argv.pop(1)
    except: 
        pass

    if do_prof:
        import pstats, cProfile
        cProfile.runctx("main()", globals(), locals(), "Profile.prof")
        s = pstats.Stats("Profile.prof")
        s.strip_dirs().sort_stats("time").print_stats()

    elif do_tracemalloc:
        # Requires py3.4
        import tracemalloc
        tracemalloc.start()

        main()

        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')

        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)
    else:
        sys.exit(main())
