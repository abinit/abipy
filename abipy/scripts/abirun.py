#!/usr/bin/env python
"""
This script allows the user to submit the calculations contained in the `Flow`.
It provides a command line interface as well as a graphical interface based on wxpython.
"""
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse
import shlex
import time

from pprint import pprint
from collections import defaultdict
from monty import termcolor
from monty.termcolor import cprint, get_terminal_size
from monty.string import make_banner
from pymatgen.io.abinit.nodes import Status
from pymatgen.io.abinit.events import autodoc_event_handlers, EventsParser
import abipy.abilab as abilab


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
    if isinstance(obj, slice) or obj is None: return obj

    try:
        # integer.
        if int(obj) == float(obj): return slice(int(obj), int(obj)+1, 1)
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


def find_flowdir_wtpos(nodepath):
    """"
    Given a directory `nodepath` containing a node of the `Flow`,
    this function locates the directory of the flow (e.g. the dir with the pickle file)
    and returns the position of the node inside the flow by parsing the directory tree

    Return: flowdir, w_pos, t_pos

        where w_pos and t_pos are the position of the work/task.
        t_pos is set to None, if we have a work.
    """
    nodepath = os.path.abspath(nodepath)

    # Is nodepath a work director?
    back, tail = os.path.split(nodepath)
    p = os.path.join(back, abilab.Flow.PICKLE_FNAME)

    if os.path.exists(p):
        # /root/flow_dir/w[num]
        head, w_dirname = os.path.split(nodepath)
        w_pos = int(w_dirname.replace("w" , ""))
        t_pos = None
        #print("got work", w_pos, t_pos)
        return p, w_pos, t_pos

    # Is nodepath a task directory?
    back, tail = os.path.split(back)
    p = os.path.join(back, abilab.Flow.PICKLE_FNAME)

    if os.path.exists(p):
        # /root/flow_dir/w[num]/t[num]
        head, t_dirname = os.path.split(nodepath)
        head, w_dirname = os.path.split(head)
        w_pos = int(w_dirname.replace("w" , ""))
        t_pos = int(t_dirname.replace("t" , ""))
        #print("got task", w_pos, t_pos)
        return p, w_pos, t_pos

    raise RuntimeError("Cannot locate flowdir from %s" % nodepath)


def selected_nids(flow, options):
    """Return the list of node ids selected by the user via the command line interface."""
    return [task.node_id for task in flow.select_tasks(nids=options.nids, wslice=options.wslice)]


def write_notebook(flow, options):
    """See http://nbviewer.ipython.org/gist/fperez/9716279"""
    from IPython.nbformat import current as nbf
    nb = nbf.new_notebook()

    cells = [
        #nbf.new_text_cell('heading', "This is an auto-generated notebook for %s" % os.path.basename(pseudopath)),
        nbf.new_code_cell("""\
##%%javascript
##IPython.OutputArea.auto_scroll_threshold = 9999;

from __future__ import print_function
from abipy import abilab
%matplotlib inline
mpld3 = abilab.mpld3_enable_notebook()

import pylab
pylab.rcParams['figure.figsize'] = (25.0, 10.0)
import seaborn as sns
#sns.set(style="dark", palette="Set2")
sns.set(style='ticks', palette='Set2')"""),

        nbf.new_code_cell("flow = abilab.Flow.pickle_load('%s')" % flow.workdir),
        nbf.new_code_cell("flow.show_dependencies()"),
        nbf.new_code_cell("flow.check_status(show=True, verbose=0)"),
        nbf.new_code_cell("flow.show_inputs(nids=None, wslice=None)"),
        nbf.new_code_cell("flow.inspect(nids=None, wslice=None)"),
        nbf.new_code_cell("flow.show_abierrors()"),
        nbf.new_code_cell("flow.show_qouts()"),
    ]

    # Now that we have the cells, we can make a worksheet with them and add it to the notebook:
    nb['worksheets'].append(nbf.new_worksheet(cells=cells))

    # Next, we write it to a file on disk that we can then open as a new notebook.
    # Note: This should be as easy as: nbf.write(nb, fname), but the current api is a little more verbose and needs a real file-like object.
    import tempfile
    _, tmpfname = tempfile.mkstemp(suffix='.ipynb', text=True)

    with open(tmpfname, 'w') as fh:
        nbf.write(nb, fh, 'ipynb')

    os.system("ipython notebook %s" % tmpfname)
    #os.execv("/Users/gmatteo/Library/Enthought/Canopy_64bit/User/bin/ipython", ["notebook %s" % tmpfname])
    #os.execv("/Users/gmatteo/Library/Enthought/Canopy_64bit/User/bin/ipython", ["notebook"])
    #os.execv("/Users/gmatteo/Library/Enthought/Canopy_64bit/User/bin/python", ["ipython", "notebook", tmpfname])


def main():

    def str_examples():
        examples = """\
usage example:
    abirun.py [FLOWDIR] rapid                    => Keep repeating, stop when no task can be executed.
    abirun.py [FLOWDIR] scheduler                => Execute flow with the scheduler
    abirun.py [FLOWDIR] events                   => Print ABINIT events (Warning/Error/Comment)
    abirun.py [FLOWDIR] history                  => Print Task history.
    abirun.py [FLOWDIR] gui                      => Open the GUI.
    abirun.py [FLOWDIR] manager slurm            => Document the TaskManager options availabe for Slurm.
    abirun.py [FLOWDIR] manager script           => Show the job script that will be produced.
    nohup abirun.py [FLOWDIR] sheduler -s 30 &   => Start the scheduler to schedule task submission.

    If FLOWDIR is not given, abirun.py automatically selects the database located within 
    the working directory. An Exception is raised if multiple databases are found.

    Note, moreover, that you can also replace FLOWDIR with the directory of a work/task
    to make the command operate on this node of the flow without having to specify --nids.
    To have the list of events of the task in `FLOWDIR/w0/t1` simply use: 

        abirun.py FLOWDIR/w0/t1 events 

    instead of 

        abirun.py FLOWDIR events -n 123

    where 123 is the node identifier of w0/t1.

    Options for developers:

        abirun.py prof ABIRUN_ARGS               => to profile abirun.py
        abirun.py tracemalloc ABIRUN_ARGS        => to trace memory blocks allocated by Python
"""
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    def parse_nids(s):
        """parse nids argument"""
        if s is None: return s
        try:
            if "," in s:
                return [int(t) for t in s.split(",")]
            else:
                # Convert string to slice and return list.
                s = as_slice(s)
                if s.stop is None: raise argparse.ArgumentTypeError("stop must be specified")
                return list(range(s.start, s.stop, s.step))
        except:
            raise argparse.ArgumentTypeError("Invalid nids string %s\n Expecting None or int or comma-separated integers or slice sintax" % s)

    def parse_wslice(s):
        s = as_slice(s)
        if s is None: return s
        if s.stop is None: raise argparse.ArgumentTypeError("stop must be specified")
        return s

    # Parent parser for commands that need to know on which subset of tasks/workflows we have to operate.
    # wslide and nids are mutually exclusive.
    flow_selector_parser = argparse.ArgumentParser(add_help=False)
    group = flow_selector_parser.add_mutually_exclusive_group()
    group.add_argument("-n", '--nids', default=None, type=parse_nids, help=(
        "Node identifier(s) used to select the task. Integer or comma-separated list of integers. Use `status` command to get the node ids."
        "Examples: --nids=12 --nids=12,13,16 --nids=10:12 to select 10 and 11, --nids=2:5:2 to select 2,4"  
        ))

    group.add_argument("-w", '--wslice', default=None, type=parse_wslice, 
                                      help=("Select the list of works to analyze (python syntax for slices):"
                                      "Examples: --wslice=1 to select the second workflow, --wslice=:3 for 0,1,2,"
                                      "--wslice=-1 for the last workflow, --wslice::2 for even indices"))

    group.add_argument("-S", '--task-status', default=None, type=Status.as_status, 
                        help="Select only the tasks with the given status. Default: None i.e. ignored. Possible values: %s" %
                        Status.all_status_strings())
    #group.add_argument("-p", "--task-pos", default=None, type=parse_wslice, help="List of tuples with the position of the tasl in the flow.")

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)

    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                              help='verbose, can be supplied multiple times to increase verbosity')
    copts_parser.add_argument('--remove-lock', default=False, type=bool, help="Remove the lock file of the pickle file storing the flow.")
    copts_parser.add_argument('--no-colors', default=False, help='Disable ASCII colors')
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                               help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('flowdir', nargs="?", help=("File or directory containing the ABINIT flow"
                                                    "If not given, the first flow in the current workdir is selected"))

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    subparsers.add_parser('version', parents=[copts_parser], help='Show version number and exit')

    # Subparser for single command.
    p_single = subparsers.add_parser('single', parents=[copts_parser], help="Run single task.")

    # Subparser for rapidfire command.
    p_rapid = subparsers.add_parser('rapid', parents=[copts_parser], help="Run all tasks in rapidfire mode")

    # Subparser for scheduler command.
    p_scheduler = subparsers.add_parser('scheduler', parents=[copts_parser], help="Run all tasks with a Python scheduler.")

    p_scheduler.add_argument('-w', '--weeks', default=0, type=int, help="number of weeks to wait")
    p_scheduler.add_argument('-d', '--days', default=0, type=int, help="number of days to wait")
    p_scheduler.add_argument('-hs', '--hours', default=0, type=int, help="number of hours to wait")
    p_scheduler.add_argument('-m', '--minutes', default=0, type=int, help="number of minutes to wait")
    p_scheduler.add_argument('-s', '--seconds', default=0, type=int, help="number of seconds to wait")

    # Subparser for batch command.
    p_batch = subparsers.add_parser('batch', parents=[copts_parser], help="Run scheduler in batch script.")
    p_batch.add_argument("-t", '--timelimit', default=None, help=("Time limit for batch script. "
                         "Accept int with seconds or string with time given in the slurm convention: "
                         "`days-hours:minutes:seconds`. If timelimit is None, the default value specified"
                         " in the `batch_adapter` entry of `manager.yml` is used."))

    # Subparser for status command.
    p_status = subparsers.add_parser('status', parents=[copts_parser, flow_selector_parser], help="Show task status.")
    p_status.add_argument('-d', '--delay', default=0, type=int, help=("If 0, exit after the first analysis.\n" + 
                          "If > 0, enter an infinite loop and delay execution for the given number of seconds."))
    p_status.add_argument('-s', '--summary', default=False, action="store_true", help="Print short version with status counters.")

    # Subparser for set_status command.
    p_set_status = subparsers.add_parser('set_status', parents=[copts_parser, flow_selector_parser], 
        help="Change the status of the task. WARNING: Option for developers!")
    p_set_status.add_argument('new_status', help="New value of status. Possible values: %s" % Status.all_status_strings())

    # Subparser for cancel command.
    p_cancel = subparsers.add_parser('cancel', parents=[copts_parser, flow_selector_parser], help="Cancel the tasks in the queue.")
    p_cancel.add_argument("-r", "--rmtree", action="store_true", default=False, help="Remove flow directory.")

    # Subparser for restart command.
    p_restart = subparsers.add_parser('restart', parents=[copts_parser, flow_selector_parser], 
                help="Restart the tasks of the flow. By default, only the task with status==Unconverged are restarted."
                     "Use -S `status` and/or -n node_ids to select particular tasks.")

    # Subparser for reset command.
    p_reset = subparsers.add_parser('reset', parents=[copts_parser, flow_selector_parser], 
                                    help="Reset the tasks of the flow with the specified status.")

    # Subparser for move command.
    p_move = subparsers.add_parser('move', parents=[copts_parser], help="Move the flow to a new directory and change the absolute paths")
    p_move.add_argument('dest', nargs=1) 

    # Subparser for open command.
    p_open = subparsers.add_parser('open', parents=[copts_parser, flow_selector_parser], 
                                   help="Open files in $EDITOR, type `abirun.py FLOWDIR open --help` for help)")
    p_open.add_argument('what', nargs="?", default="o", 
        help="""\
Specify the files to open. Possible choices:
    i ==> input_file
    o ==> output_file
    f ==> files_file
    j ==> job_file
    l ==> log_file
    e ==> stderr_file
    q ==> qout_file
    all ==> all files.
""")

    # Subparser for ncopen.
    p_ncopen = subparsers.add_parser('ncopen', parents=[copts_parser, flow_selector_parser], 
                                      help="Open netcdf files in ipython. Use --help` for more info")
    p_ncopen.add_argument('ncext', nargs="?", default="GSR", help="Select the type of file to open")

    # Subparser for gui command.
    p_gui = subparsers.add_parser('gui', help="Open the GUI (requires wxPython).")
    p_gui.add_argument("--chroot", default="", type=str, help=("Use chroot as new directory of the flow." +
                       "Mainly used for opening a flow located on a remote filesystem mounted with sshfs." +
                       "In this case chroot is the absolute path to the flow on the **localhost**",
                       "Note that it is not possible to change the flow from remote when chroot is used."))

    # Subparser for new_manager.
    p_new_manager = subparsers.add_parser('new_manager', parents=[copts_parser, flow_selector_parser], help="Change the TaskManager.")
    p_new_manager.add_argument("manager_file", default="", type=str, help="YAML file with the new manager")

    # Subparser for tail.
    p_tail = subparsers.add_parser('tail', parents=[copts_parser, flow_selector_parser], help="Use tail to follow the main output files of the flow.")
    p_tail.add_argument('what_tail', nargs="?", type=str, default="o", help="What to follow: o for output (default), l for logfile, e for stderr")

    # Subparser for qstat.
    p_qstat = subparsers.add_parser('qstat', parents=[copts_parser], help="Show additional info on the jobs in the queue.")

    # Subparser for deps.
    p_deps = subparsers.add_parser('deps', parents=[copts_parser], help="Show dependencies.")

    # Subparser for robot.
    p_robot = subparsers.add_parser('robot', parents=[copts_parser, flow_selector_parser], 
                                    help="Use a robot to analyze the results of multiple tasks (requires ipython)")
    p_robot.add_argument('robot_ext', nargs="?", type=str, default="GSR", help="The file extension of the netcdf file")

    # Subparser for plot.
    p_plot = subparsers.add_parser('plot', parents=[copts_parser, flow_selector_parser], help="Plot data. Use --help for more info.")
    p_plot.add_argument("what", nargs="?", type=str, default="ebands", help="Object to plot")

    # Subparser for inspect.
    p_inspect = subparsers.add_parser('inspect', parents=[copts_parser, flow_selector_parser], help="Inspect the tasks")

    # Subparser for inputs.
    p_inputs= subparsers.add_parser('inputs', parents=[copts_parser, flow_selector_parser], help="Show the input files of the tasks")

    # Subparser for manager.
    p_manager = subparsers.add_parser('manager', parents=[copts_parser], help="Document the TaskManager options")
    p_manager.add_argument("qtype", nargs="?", default=None, help=("Write job script to terminal if qtype='script' else" 
        " document the qparams for the given QueueAdapter qtype e.g. slurm"))

    # Subparser for events.
    p_events = subparsers.add_parser('events', parents=[copts_parser, flow_selector_parser], 
                                    help="Show ABINIT events (error messages, warnings, comments)")
    #p_events.add_argument("-t", "event-type", default=)

    # Subparser for corrections.
    p_corrections = subparsers.add_parser('corrections', parents=[copts_parser, flow_selector_parser], help="Show abipy corrections")

    # Subparser for history.
    p_history = subparsers.add_parser('history', parents=[copts_parser, flow_selector_parser], help="Show Node history.")
    p_history.add_argument("-m", "--metadata", action="store_true", default=False, help="Print history metadata")
    #p_history.add_argument("-t", "--task-history", action="store_true", default=True, help=)

    # Subparser for handlers.
    p_handlers = subparsers.add_parser('handlers', parents=[copts_parser], help="Show event handlers installed in the flow")
    p_handlers.add_argument("-d", "--doc", action="store_true", default=False, 
                            help="Show documentation about all the handlers that can be installed.")

    # Subparser for notebook.
    p_notebook = subparsers.add_parser('notebook', parents=[copts_parser], help="Create and open an ipython notebook to interact with the flow.")

    # Subparser for ipython.
    p_ipython = subparsers.add_parser('ipython', parents=[copts_parser], help="Embed IPython. Useful for advanced operations or debugging purposes.")
    p_ipython.add_argument('--argv', nargs="?", default="", type=shlex.split, 
                           help="Command-line options passed to ipython. Must be enclosed by quotes. "
                                "Example: --argv='--matplotlib=wx'")

    # Subparser for tar.
    p_tar = subparsers.add_parser('tar', parents=[copts_parser], help="Create tarball file.")
    p_tar.add_argument("-s", "--max-filesize", default=None, 
                       help="Exclude file whose size > max-filesize bytes. Accept integer or string e.g `1Mb`.")

    def parse_strings(s): return s.split(",") if s is not None else s
    p_tar.add_argument("-e", "--exclude-exts", default=None, type=parse_strings,
                       help="Exclude file extensions. Accept string or comma-separated strings. Ex: -eWFK or --exclude-exts=WFK,GSR")

    p_tar.add_argument("-d", "--exclude-dirs", default=None, type=parse_strings,
                       help="Exclude directories. Accept string or comma-separated strings. Ex: --exlude-dirs=indir,outdir")

    p_tar.add_argument("-l", "--light", default=False, action="store_true",
                       help="Create light-weight version of the tarball for debugging purposes. Other options are ignored.")

    # Subparser for debug.
    p_debug = subparsers.add_parser('debug', parents=[copts_parser, flow_selector_parser], 
                                     help="Scan error files and log files for possible error messages.")

    # Subparser for group.
    p_group = subparsers.add_parser('group', parents=[copts_parser, flow_selector_parser], 
                                     help="Group tasks according to property.")

    # Subparser for diff.
    p_diff = subparsers.add_parser('diff', parents=[copts_parser, flow_selector_parser], 
                                   help="Compare files produced by two or three nodes.")
    p_diff.add_argument('what_diff', nargs="?", type=str, default="i", 
                        help="What to diff: i for input (default), o for output, l for logfile, e for stderr")

    # Subparser for networkx.
    p_networkx = subparsers.add_parser('networkx', parents=[copts_parser], #, flow_selector_parser], 
                                     help="Draw flow and node dependecies with networkx package.")
    p_networkx.add_argument('--nxmode', default="status",
                            help="Type of network plot. Possible values: `status`, `network`")
    p_networkx.add_argument('--edge-labels', action="store_true", default=False, help="Show edge labels")

    # Subparser for listext.
    p_listext = subparsers.add_parser('listext', parents=[copts_parser], 
                                     help="List all the output files with the given extension that have been produced by the nodes of the flow.")
    p_listext.add_argument('listexts', nargs="+", help="List of Abinit file extensions. e.g DDB, GSR, WFK etc")

    # Subparser for timer.
    p_timer = subparsers.add_parser('timer', parents=[copts_parser, flow_selector_parser], 
                                    help=("Read the section with timing info from the main ABINIT output file (requires timopt != 0)"
                                          "Open Ipython terminal to inspect data."))

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

    if options.command == "version":
        from abipy.core.release import version
        print(version)
        return 0

    if options.no_colors:
        # Disable colors
        termcolor.enable(False)

    if options.command == "manager":
        # Document TaskManager options and qparams.
        qtype = options.qtype

        if qtype == "script":
            manager = abilab.TaskManager.from_user_config()
            script = manager.qadapter.get_script_str(
                job_name="job_name", 
                launch_dir="workdir",
                executable="executable",
                qout_path="qout_file.path",
                qerr_path="qerr_file.path",
                stdin="stdin", 
                stdout="stdout",
                stderr="stderr",
            )
            print(script)

        else:
            print(abilab.TaskManager.autodoc())
            from pymatgen.io.abinit.qadapters import show_qparams, all_qtypes
                                                                                                 
            print("qtype supported: %s" % all_qtypes())
            print("Use `abirun.py . manager slurm` to have the list of qparams for slurm.\n")

            if qtype is not None:
                print("QPARAMS for %s" % qtype)
                show_qparams(qtype)

        sys.exit(0)


    patch_nids = False
    if options.flowdir is None:
        # Will try to figure out the location of the Flow.
        options.flowdir = os.getcwd()
    else:
        # Sometimes one wants to inspect a work or a task by just using `abirun.py flow/w0/t0 inspect`
        # without knowing its node id. find_flowdir_wtpos will solve the problem!
        if not os.path.exists(os.path.join(options.flowdir, abilab.Flow.PICKLE_FNAME)):
            print("The directory does not contain a flow. Will get node ids from dirpath.")
            #assert options.nids is None
            patch_nids = True
            options.flowdir, w_pos, t_pos = find_flowdir_wtpos(options.flowdir)

    # Read the flow from the pickle database.
    flow = abilab.Flow.pickle_load(options.flowdir, remove_lock=options.remove_lock)
    #flow.set_spectator_mode(False)

    if patch_nids:
        # Create options.nids here  
        node = flow[w_pos]
        if t_pos is not None: node = node[t_pos]
        options.nids = [node.node_id]

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

        # Default status for new_manager is QCritical
        if options.task_status is None: options.task_status = Status.as_status("QCritical")

        # Change the manager of the errored tasks.
        print("Resetting tasks with status: %s" % options.task_status)
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            task.reset()
            task.set_manager(new_manager)
            
        # Update the database.
        return flow.build_and_pickle_dump()

    elif options.command == "events":
        nrows, ncols = get_terminal_size()

        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            report = task.get_event_report()
            print(make_banner(str(task), width=ncols, mark="="))
            #report = report.filter_types()
            print(report)

    elif options.command == "corrections":
        nrows, ncols = get_terminal_size()
        count = 0
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            if task.num_corrections == 0: continue
            count += 1
            print(make_banner(str(task), width=ncols, mark="="))
            for corr in task.corrections:
                pprint(corr)

        if not count: 
            print("No correction found.")

    elif options.command == "history":
        nrows, ncols = get_terminal_size()

        works_done = []
        # Loop on the tasks and show the history of the work is not in works_done
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            work = task.work
            if work not in works_done:
                works_done.append(work)
                print(make_banner(str(work), width=ncols, mark="="))
                print(work.history.to_string(metadata=options.metadata))

            print(make_banner(str(task), width=ncols, mark="="))
            print(task.history.to_string(metadata=options.metadata))

        # Print the history of the flow.
        print(make_banner(str(flow), width=ncols, mark="="))
        print(flow.history.to_string(metadata=options.metadata))

    elif options.command == "handlers":
        if options.doc:
            autodoc_event_handlers()
        else:
            flow.show_event_handlers()

    elif options.command  == "single":
        nlaunch = flow.single_shot()
        print("Number of tasks launched: %d" % nlaunch)
        if nlaunch: flow.show_status()

    elif options.command == "rapid":
        nlaunch = flow.rapidfire()
        print("Number of tasks launched: %d" % nlaunch)
        if nlaunch: flow.show_status()

    elif options.command == "scheduler":
        # Check that the env on the local machine is properly configured before starting the scheduler.
        abilab.abicheck()

        sched_options = {oname: getattr(options, oname) for oname in 
            ("weeks", "days", "hours", "minutes", "seconds")}

        if all(v == 0 for v in sched_options.values()):
            sched = flow.make_scheduler()
        else:
            sched = flow.make_scheduler(**sched_options)

        print(sched)
        return sched.start()

    elif options.command == "batch":
        return flow.batch(timelimit=options.timelimit)

    elif options.command == "status":
        # Select the method to call.
        show_func = flow.show_status if not options.summary else flow.show_summary

        if options.delay:
            cprint("Entering infinite loop. Press CTRL+C to exit", color="magenta", end="", flush=True)
            try:
                while True:
                    print(2*"\n" + time.asctime() + "\n")
                    flow.check_status()
                    show_func(verbose=options.verbose, nids=selected_nids(flow, options))
                    if flow.all_ok: break
                    time.sleep(options.delay)
            except KeyboardInterrupt:
                pass
        else:
            show_func(verbose=options.verbose, nids=selected_nids(flow, options))
            if options.verbose and flow.manager.has_queue:
                print("Total number of jobs in queue: %s" % flow.manager.get_njobs_in_queue())

    elif options.command == "set_status":
        # Default status for reset is QCritical
        if options.task_status is None: options.task_status = Status.as_status("QCritical")
        new_status = Status.as_status(options.new_status)
        print("Will set all tasks with status: ", options.task_status, " to new_status", new_status)

        count = 0
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            task.set_status(new_status, msg="Changed by abirun from %s to %s" % (task.status, new_status))
            count += 1

        print("Number of tasks modified: %s" % count)
        if count:
            # update database
            flow.pickle_dump()

    elif options.command == "open":
        flow.open_files(what=options.what, status=None, op="==", nids=selected_nids(flow, options))

    elif options.command == "ncopen":
        # The name of the method associated to this netcdf file.
        methname = "open_" + options.ncext.lower()
        # List of netcdf file objects.
        ncfiles = [getattr(task, methname)() for task in flow.select_tasks(nids=options.nids, wslice=options.wslice) 
                    if hasattr(task, methname)]
        
        if ncfiles:
            # Start ipython shell with namespace 
            import IPython
            if len(ncfiles) == 1:
                IPython.start_ipython(argv=[], user_ns={"ncfile": ncfiles[0]})
            else:
                IPython.start_ipython(argv=[], user_ns={"ncfiles": ncfiles})
        else:
            cprint("Cannot find any netcdf file with extension %s" % options.ncext, color="magenta")

    elif options.command == "cancel":
        print("Number of jobs cancelled %d" % flow.cancel(nids=selected_nids(flow, options)))
        # Remove directory
        if options.rmtree: flow.rmtree()

    elif options.command == "restart":
        # Default status for reset is Unconverged if no option is provided by the user.
        if options.task_status is None and options.nids is None:
            options.task_status = Status.as_status("Unconverged")

        nlaunch, excs = 0, []
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            if options.verbose: 
                print("Will try to restart %s, with status %s" % (task, task.status))
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
        # Default status for reset is QCritical
        if options.task_status is None: options.task_status = Status.as_status("QCritical")
        print("Will reset tasks with status: %s" % options.task_status)

        count = 0
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            print("Resetting task %s" % task)
            failed = task.reset()
            if failed:
                print("Task %s couldn't be reset" % task)
            else:
                count += 1

        cprint("%d tasks have been reset" % count, "blue")
        nlaunch = flow.rapidfire()
        flow.show_status()
        print("Number of tasks launched: %d" % nlaunch)

        if nlaunch == 0:
            g = flow.find_deadlocks()
            #print("deadlocked:", gdeadlocked)
            #print("runnables:", grunnables)
            #print("running:", g.running)
            if g.deadlocked and not (g.runnables or g.running):
                print("*** Flow is deadlocked ***")

        flow.pickle_dump()

    elif options.command == "move":
        print("Will move flow to %s..." % options.dest)
        flow.chroot(options.dest)
        flow.move(options.dest)

    elif options.command == "tail":
        def get_path(task):
            """Helper function used to select the files of a task."""
            choices = {
                "o": task.output_file,
                "l": task.log_file,
                "e": task.stderr_file,
            }
            return getattr(choices[options.what_tail], "path")

        # Default status for tail is Running
        if options.task_status is None: options.task_status = Status.as_status("Running")

        paths = [get_path(task) for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options))]

        if not paths:
            cprint("No job is running. Exiting!", "red")
        else:
            cprint("Press CTRL+C to interrupt. Number of output files %d" % len(paths), color="magenta", end="", flush=True)
            try:
                os.system("tail -f %s" % " ".join(paths))
            except KeyboardInterrupt:
                pass

    elif options.command == "qstat":
        #for task in flow.select_tasks(nids=options.nids, wslice=options.wslice):
        for task in flow.iflat_tasks():
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
        with abilab.abirobot(flow, options.robot_ext, nids=selected_nids(flow, options)) as robot:
            #IPython.embed(header=str(robot) + "\nType `robot` in the terminal and use <TAB> to list its methods",  robot=robot)
            IPython.start_ipython(argv=[], user_ns={"robot": robot})

    elif options.command == "plot":
        fext = dict(
            ebands="gsr",
        )[options.what]

        open_method = "open_" + fext
        plot_method = "plot_" + options.what

        for task in flow.select_tasks(nids=options.nids, wslice=options.wslice):
            try:
                with getattr(task, open_method)() as ncfile: 
                    getattr(ncfile, plot_method)()
            except Exception as exc:
                print(exc)

    elif options.command == "inspect":
        tasks = flow.select_tasks(nids=options.nids, wslice=options.wslice)

        # Use different thread to inspect the task so that master can catch KeyboardInterrupt and exit.
        # One could use matplotlib non-blocking interface with show(block=False) but this one seems to work well.
        from multiprocessing import Process

        def plot_graphs():
            for task in tasks:
                if hasattr(task, "inspect"):
                    try:
                        task.inspect()
                    except Exception as exc:
                        cprint("%s: inspect method raised %s " % (task, exc), color="blue")
                        
                else:
                    cprint("Task %s does not provide an inspect method" % task, color="blue")

        plot_graphs()

        # This works with py3k but not with py2
        #p = Process(target=plot_graphs)
        #p.start()
        #num_tasks = len(tasks)

        #if num_tasks == 1:
        #    p.join()
        #else:
        #    cprint("Will produce %d matplotlib plots. Press CTRL+C to interrupt..." % num_tasks, color="magenta", end="", flush=True)
        #    try:
        #        p.join()
        #    except KeyboardInterrupt:
        #        print("\nTerminating thread...")
        #        p.terminate()

    elif options.command == "inputs":
        flow.show_inputs(nids=selected_nids(flow, options))

    elif options.command == "notebook":
        write_notebook(flow, options)

    elif options.command == "ipython":
        import IPython
        #IPython.embed(header="")
        #print("options:", options.argv)
        IPython.start_ipython(argv=options.argv, user_ns={"flow": flow})# , header="flow.show_status()")

    elif options.command == "tar":
        if not options.light:
            tarfile = flow.make_tarfile(name=None, 
                                        max_filesize=options.max_filesize, 
                                        exclude_exts=options.exclude_exts, 
                                        exclude_dirs=options.exclude_dirs,
                                        verbose=options.verbose)
            print("Created tarball file %s" % tarfile)
        else:
            tarfile = flow.make_light_tarfile()
            print("Created light tarball file %s" % tarfile)

    elif options.command == "debug":
        nrows, ncols = get_terminal_size()

        # Test for scheduler exceptions first.
        sched_excfile = os.path.join(flow.workdir, "_exceptions")
        if os.path.exists(sched_excfile):
            with open(sched_excfile, "r") as fh:
                cprint(fh.read(), color="red")
                return 0

        if options.task_status is not None: 
            tasks = list(flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)))
        else:
            errors = list(flow.iflat_tasks(status=flow.S_ERROR, nids=selected_nids(flow, options)))
            qcriticals = list(flow.iflat_tasks(status=flow.S_QCRITICAL, nids=selected_nids(flow, options)))
            abicriticals = list(flow.iflat_tasks(status=flow.S_ABICRITICAL, nids=selected_nids(flow, options)))
            tasks = errors + qcriticals + abicriticals

        # For each task selected:
        #
        #     1) Check the error files of the task. If not empty, print the content to stdout and we are done.
        #     2) If error files are empty, look at the master log file for possible errors 
        #     3) If also this check failes, scan all the process log files.
        #        TODO: This check is not needed if we introduce a new __abinit_error__ file 
        #        that is created by the first MPI process that invokes MPI abort!
        #     
        ntasks = 0
        for task in tasks:
            print(make_banner(str(task), width=ncols, mark="="))
            ntasks += 1

            #  Start with error files.
            for efname in ["qerr_file", "stderr_file",]:
                err_file = getattr(task, efname)
                if err_file.exists:
                    s = err_file.read()
                    if not s: continue
                    print(make_banner(str(err_file), width=ncols, mark="="))
                    cprint(s, color="red")
                    #count += 1 

            # Check main log file.
            try:
                report = task.get_event_report()
                if report and report.num_errors: 
                    print(make_banner(os.path.basename(report.filename), width=ncols, mark="="))
                    s = "\n".join(str(e) for e in report.errors)
                else:
                    s = None
            except Exception as exc:
                s = str(exc)

            count = 0 # count > 0 means we found some useful info that could explain the failures.
            if s is not None:
                cprint(s, color="red")
                count += 1

            if not count:
                # Inspect all log files produced by the other nodes.
                log_files = task.tmpdir.list_filepaths(wildcard="*LOG_*")
                if not log_files:
                    cprint("No *LOG_* file in tmpdir. This usually happens if you are running with many CPUs", color="magenta")

                for log_file in log_files:
                    try:
                        report = EventsParser().parse(log_file)
                        if report.errors:
                            print(report)
                            count += 1
                            break
                    except Exception as exc:
                        cprint(str(exc), color="red")
                        count += 1
                        break

            if not count:
                cprint("Houston, we could not find any error message that can explain the problem", color="magenta")

        print("Number of tasks analyzed: %d" % ntasks)

    elif options.command == "group":
        d = defaultdict(list)
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            key = task.status
            d[key].append(task.node_id)

        print("Mapping status --> List of node identifiers")
        for k, v in d.items():
            print("   ",k, " --> ", v)

    elif options.command == "diff":
        if options.nids is None:
            raise ValueError("nids must be specified when using diff command")

        tasks = list(flow.iflat_tasks(nids=selected_nids(flow, options)))

        if len(tasks) not in (2, 3):
            if len(tasks) == 1: 
                cprint("task == task, returning\n" , color="magenta", end="", flush=True)
                return 0
            else:
                raise ValueError("Don't know how to compare files produced by %d tasks" % len(tasks))

        # Build list of lists. Each sub-list contains the files associated to the i-th task.
        files_for_task = [None] * len(tasks)
        for i, task in enumerate(tasks):
            files_for_task[i] = task.select_files(options.what_diff)

        for diff_files in zip(*files_for_task):
            print("Comparing", ", ".join(os.path.relpath(p) for p in diff_files))
            args = " ".join(os.path.relpath(p) for p in diff_files)
            # TODO: I should have written a Differ object somewhere!
            os.system("vimdiff %s" % args)

    elif options.command == "networkx":
        flow.plot_networkx(mode=options.nxmode, 
                           with_edge_labels=options.edge_labels)

    elif options.command == "listext":
        for ext in options.listexts:
            flow.listext(ext)
            print("")

    elif options.command == "timer":
        print("Warning this option is still under development")
        timer = flow.get_abitimer()
        if timer is None:
            cprint("Cannot parse time data!", color="magenta", end="", flush=True)
            return 1

        import IPython
        IPython.start_ipython(argv=[], user_ns={"timer": timer})

    else:
        raise RuntimeError("Don't know what to do with command %s!" % options.command)

    return retcode

if __name__ == "__main__":
    # Replace python open to detect open files.
    #from abipy.tools import open_hook
    #open_hook.install()
    retcode = 0
    do_prof, do_tracemalloc = 2* [False]
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

        retcode = main()

        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')

        print("[Top 10]")
        for stat in top_stats[:10]:
            print(stat)
    else:
        sys.exit(main())

    #open_hook.print_open_files()
    sys.exit(retcode)
