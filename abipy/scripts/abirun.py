#!/usr/bin/env python
"""
This script allows the user to submit the calculations contained in the `Flow`.
It provides a command line interface as well as a graphical interface based on wxpython.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import argparse
import shlex
import time
import platform
import tempfile
import numpy as np
import abipy.flowtk as flowtk
import abipy.abilab as abilab

from pprint import pprint
from collections import defaultdict, OrderedDict
from socket import gethostname
from monty import termcolor
from monty.os.path import which
from monty.functools import prof_main
from monty.termcolor import cprint, colored, get_terminal_size
from monty.string import boxed, list_strings, make_banner
from pymatgen.util.io_utils import ask_yesno
from abipy.tools import duck
from abipy.flowtk import Status
from abipy.core.structure import dataframes_from_structures


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def parse_strings(s):
    """Parse comma separated values. Return None if s is None."""
    return s.split(",") if s is not None else s


def flowdir_wname_tname(dirname):
    """"
    Given a initial directory `dirname` containing a node of the `Flow`,
    this function locates the directory of the flow (e.g. the directory with the pickle file)
    and returns the name of the work and/or of the node.

    Return: flowdir, wname, tname

    where flowdir is the directory containing the pickle file,
    wname and tname are the basenames of the work/task.

    If dirname contains the pickle file we have (wname, tname) == (None, None)
    If dirname is a work --> wname is it's basename and tname is None
    If dirname is a task --> os.path.join(flowdir, wname, tname) == task.workdir.
    """
    if dirname is None: dirname = os.getcwd()
    dirname = os.path.abspath(dirname)
    if os.path.exists(os.path.join(dirname, flowtk.Flow.PICKLE_FNAME)):
        return dirname, None, None

    # Handle works or tasks.
    head = dirname
    wname, tname = None, None
    for i in range(2):
        head, tail = os.path.split(head)
        if i == 0: tail_1 = tail
        if os.path.exists(os.path.join(head, flowtk.Flow.PICKLE_FNAME)):
            if i == 0:
                # We have a work: /root/flow_dir/w[num]
                wname = tail
            if i == 1:
                # We have a task: /root/flow_dir/w[num]/t[num]
                wname = tail
                tname = tail_1

            #print("wname", wname, "tname", tname)
            return head, wname, tname

    raise RuntimeError("Cannot locate flowdir from %s" % dirname)


def selected_nids(flow, options):
    """
    Return the list of node ids selected by the user via the command line interface.
    """
    task_ids = [task.node_id for task in
                flow.select_tasks(nids=options.nids, wslice=options.wslice, task_class=options.task_class)]

    # Have to add the ids of the works containing the tasks.
    if options.nids is not None:
        work_ids = [work.node_id for work in flow if work.node_id in options.nids]
    else:
        work_ids = [work.node_id for work in flow]

    return set(work_ids + task_ids)


def cli_abiopen(options, filepath):
    """Code taken from abiopen.py."""
    # TODO: One should have a single API but make_and_open_notebook are different so
    # for the time being we use two different versions.
    filepath = options.flowdir
    if filepath is None:
        raise RuntimeError("filepath must be specified")

    if not os.path.exists(filepath):
        raise RuntimeError("%s: no such file" % filepath)

    if not options.notebook:
        # Start ipython shell with namespace
        abifile = abilab.abiopen(filepath)
        if options.print:
            if hasattr(abifile, "to_string"):
                print(abifile.to_string(verbose=options.verbose))
            else:
                print(abifile)
            return 0

        import IPython
        # Use embed because I don't know how to show a header with start_ipython.
        IPython.embed(header="The Abinit file is bound to the `abifile` variable.\nTry `print(abifile)`")

    else:
        # Call specialized method if the object is a NotebookWriter
        # else generate simple notebook by calling `make_and_open_notebook`
        cls = abilab.abifile_subclass_from_filename(filepath)
        if hasattr(cls, "make_and_open_notebook"):
            if hasattr(cls, "__exit__"):
                with abilab.abiopen(filepath) as abifile:
                    return abifile.make_and_open_notebook(foreground=options.foreground)
            else:
                abifile = abilab.abiopen(filepath)
                return abifile.make_and_open_notebook(foreground=options.foreground)
        else:
            return make_and_open_notebook(options)

# TODO: These should become flow methods.
def flow_write_open_notebook(flow, options):
    """
    Generate an ipython notebook and open it in the browser.
    Return system exit code.
    """
    import nbformat
    nbf = nbformat.v4
    nb = nbf.new_notebook()

    nb.cells.extend([
        #nbf.new_markdown_cell("This is an auto-generated notebook for %s" % os.path.basename(pseudopath)),
        nbf.new_code_cell("""\
from __future__ import print_function, division, unicode_literals, absolute_import

import sys, os
import numpy as np

%matplotlib notebook
from IPython.display import display

# This to render pandas DataFrames with https://github.com/quantopian/qgrid
#import qgrid
#qgrid.nbinstall(overwrite=True)  # copies javascript dependencies to your /nbextensions folder

from abipy import abilab

# Tell AbiPy we are inside a notebook and use seaborn settings for plots.
# See https://seaborn.pydata.org/generated/seaborn.set.html#seaborn.set
abilab.enable_notebook(with_seaborn=True)

# AbiPy widgets for pandas and seaborn plot APIs
#import abipy.display.seabornw import snw
#import abipy.display.pandasw import pdw
"""),

        nbf.new_code_cell("flow = abilab.Flow.pickle_load('%s')" % flow.workdir),
        nbf.new_code_cell("if flow.num_errored_tasks: flow.debug()"),
        nbf.new_code_cell("flow.check_status(show=True, verbose=0)"),
        nbf.new_code_cell("flow.show_dependencies()"),
        nbf.new_code_cell("flow.plot_networkx();"),
        nbf.new_code_cell("#flow.get_graphviz();"),
        nbf.new_code_cell("flow.show_inputs(nids=None, wslice=None)"),
        nbf.new_code_cell("flow.show_history()"),
        nbf.new_code_cell("flow.show_corrections()"),
        nbf.new_code_cell("flow.show_event_handlers()"),
        nbf.new_code_cell("flow.inspect(nids=None, wslice=None)"),
        nbf.new_code_cell("flow.show_abierrors()"),
        nbf.new_code_cell("flow.show_qouts()"),
    ])

    import tempfile, io
    _, nbpath = tempfile.mkstemp(suffix='.ipynb', text=True)

    with io.open(nbpath, 'wt', encoding="utf8") as fh:
        nbformat.write(nb, fh)

    if which("jupyter") is None:
        raise RuntimeError("Cannot find jupyter in $PATH. Install it with `pip install`")

    if options.foreground:
        return os.system("jupyter notebook %s" % nbpath)
    else:
        fd, tmpname = tempfile.mkstemp(text=True)
        print(tmpname)
        cmd = "jupyter notebook %s" % nbpath
        print("Executing:", cmd)
        print("stdout and stderr redirected to %s" % tmpname)
        import subprocess
        process = subprocess.Popen(cmd.split(), shell=False, stdout=fd, stderr=fd)
        cprint("pid: %s" % str(process.pid), "yellow")


def flow_compare_structures(flow, nids=None, with_spglib=False, what="io", verbose=0,
                            precision=3, printout=False, with_colors=False):
    """
    Analyze structures of the tasks (input and output structures if it's a relaxation
    task. Print pandas DataFrame

    Args:
        nids: List of node identifiers. By defaults all nodes are shown
        with_spglib: If True, spglib is invoked to get the spacegroup symbol and number
        what (str): "i" for input structures, "o" for output structures.
        precision: Floating point output precision (number of significant digits).
            This is only a suggestion
        printout: True to print dataframe.
        with_colors: True if task status should be colored.
    """
    #flow.check_status()
    structures, index, status, max_forces, pressures, task_classes = [], [], [], [], [], []

    def push_data(post, task, structure, cart_forces, pressure):
        """Helper function to fill lists"""
        index.append(task.pos_str + post)
        structures.append(structure)
        status.append(task.status.colored if with_colors else str(task.status))
        if cart_forces is not None:
            fmods = np.sqrt([np.dot(f, f) for f in cart_forces])
            max_forces.append(fmods.max())
        else:
            max_forces.append(None)
        pressures.append(pressure)
        task_classes.append(task.__class__.__name__)

    for task in flow.iflat_tasks(nids=nids):
        if "i" in what:
            push_data("_in", task, task.input.structure, cart_forces=None, pressure=None)

        if "o" not in what:
            continue

        # Add final structure, pressure and max force if relaxation task or GS task
        if task.status in (task.S_RUN, task.S_OK):
            if hasattr(task, "open_hist"):
                # Structural relaxations produce HIST.nc and we can get
                # the final structure or the structure of the last relaxation step.
                try:
                    with task.open_hist() as hist:
                        final_structure = hist.final_structure
                        stress_cart_tensors, pressures_hist = hist.reader.read_cart_stress_tensors()
                        forces = hist.reader.read_cart_forces(unit="eV ang^-1")[-1]
                        push_data("_out", task, final_structure, forces, pressures_hist[-1])
                except Exception as exc:
                    cprint("Exception while opening HIST.nc file of task: %s\n%s" % (task, str(exc)), "red")

            elif hasattr(task, "open_gsr") and task.status == task.S_OK and task.input.get("iscf", 7) >= 0:
                with task.open_gsr() as gsr:
                    forces = gsr.reader.read_cart_forces(unit="eV ang^-1")
                    push_data("_out", task, gsr.structure, forces, gsr.pressure)

    dfs = dataframes_from_structures(structures, index=index, with_spglib=with_spglib, cart_coords=False)

    if any(f is not None for f in max_forces):
        # Add pressure and forces to the dataframe
        dfs.lattice["P [GPa]"] = pressures
        dfs.lattice["Max|F| eV/ang"] = max_forces

    # Add columns to the dataframe.
    status = [str(s) for s in status]
    dfs.lattice["task_class"] = task_classes
    dfs.lattice["status"] = dfs.coords["status"] = status

    if printout:
        abilab.print_dataframe(dfs.lattice, title="Lattice parameters:", precision=precision)
        if verbose:
            abilab.print_dataframe(dfs.coords, title="Atomic positions (columns give the site index):")
        else:
            print("Use `--verbose` to print atoms.")

    return dfs


def flow_compare_ebands(flow, nids=None, with_spglib=False, verbose=0,
                        precision=3, printout=False, with_colors=False, plot_mode=None):
    """
    Analyze electron bands produced by the tasks. Print pandas DataFrame

    Args:
        nids: List of node identifiers. By default, all nodes are shown
        with_spglib: If True, spglib is invoked to get the spacegroup symbol and number
        precision: Floating point output precision (number of significant digits).
            This is only a suggestion
        printout: True to print dataframe.
        with_colors: True if task status should be colored.
        plot_mode: Plot results if not None. Allowed value in ["gridplot", "combiplot"]
    """
    #flow.check_status()
    ebands_list, index, status, ncfiles, task_classes, task_nids = [], [], [], [], [], []

    # Cannot use robots because ElectronBands can be found in different filetypes
    for task in flow.iflat_tasks(nids=nids, status=flow.S_OK):
        # Read ebands from GSR or SIGRES files.
        for ext in ("gsr", "sigres"):
            task_open_ncfile = getattr(task, "open_%s" % ext, None)
            if task_open_ncfile is not None: break
        else:
            continue

        try:
            with task_open_ncfile() as ncfile:
                ebands_list.append(ncfile.ebands)
                index.append(task.pos_str)
                status.append(task.status.colored if with_colors else str(task.status))
                ncfiles.append(os.path.relpath(ncfile.filepath))
                task_classes.append(task.__class__.__name__)
                task_nids.append(task.node_id)
        except Exception as exc:
            cprint("Exception while opening HIST.nc file of task: %s\n%s" % (task, str(exc)), "red")

    if not ebands_list: return
    df = abilab.dataframe_from_ebands(ebands_list, index=index, with_spglib=with_spglib)

    # Add columns to the dataframe.
    status = [str(s) for s in status]
    df["task_class"] = task_classes
    df["ncfile"] = ncfiles
    df["node_id"] = task_nids
    df["status"] = status

    if printout:
        abilab.print_dataframe(df, title="KS electronic bands:", precision=precision)

    if plot_mode is not None:
        plotter = abilab.ElectronBandsPlotter(key_ebands=zip(ncfiles, ebands_list))
        plotfunc = getattr(plotter, plot_mode, None)
        if plotfunc is None:
            raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)
        plotfunc(tight_layout=True)

    return df


def flow_compare_hist(flow, nids=None, with_spglib=False, verbose=0,
                      precision=3, printout=False, with_colors=False, plot_mode=None):
    """
    Analyze HIST nc files produced by the tasks. Print pandas DataFrame with final results.

    Args:
        nids: List of node identifiers. By defaults all nodes are shown
        with_spglib: If True, spglib is invoked to get the spacegroup symbol and number
        precision: Floating point output precision (number of significant digits).
            This is only a suggestion
        printout: True to print dataframe.
        with_colors: True if task status should be colored.
        plot_mode: Plot results if not None. Allowed value in ["gridplot", "combiplot"]
    """
    #flow.check_status()
    hist_paths, index, status, ncfiles, task_classes, task_nids = [], [], [], [], [], []

    for task in flow.iflat_tasks(nids=nids):
        if task.status not in (flow.S_OK, flow.S_RUN): continue
        hist_path = task.outdir.has_abiext("HIST")
        if not hist_path: continue

        hist_paths.append(hist_path)
        index.append(task.pos_str)
        status.append(task.status.colored if with_colors else str(task.status))
        ncfiles.append(os.path.relpath(hist_path))
        task_classes.append(task.__class__.__name__)
        task_nids.append(task.node_id)

    if not hist_paths: return
    robot = abilab.HistRobot.from_files(hist_paths, labels=hist_paths)
    df = robot.get_dataframe(index=index, with_spglib=with_spglib)

    # Add columns to the dataframe.
    status = [str(s) for s in status]
    df["task_class"] = task_classes
    df["ncfile"] = ncfiles
    df["node_id"] = task_nids
    df["status"] = status

    if printout:
        title = "Table with final structures, pressures in GPa and force stats in eV/Ang:\n"
        abilab.print_dataframe(df, title=title, precision=precision)

    #print("plot_mode", plot_mode)
    if plot_mode is not None:
        if len(robot) == 1:
            robot.abifiles[0].plot()
        else:
            if plot_mode == "gridplot":
                for what in robot.what_list:
                    robot.gridplot(what=what)
            elif plot_mode == "combiplot":
                robot.combiplot()
            else:
                raise ValueError("Invalid value of plot_mode: %s" % str(plot_mode))

    return df


def flow_get_dims_dataframe(flow, nids=None, printout=False, with_colors=False):
    """
    Analyze output files produced by Abinit tasks. Print pandas DataFrame with dimensions.

    Args:
        nids: List of node identifiers. By defaults all nodes are shown
        printout: True to print dataframe.
        with_colors: True if task status should be colored.
    """
    abo_paths, index, status, abo_relpaths, task_classes, task_nids = [], [], [], [], [], []

    for task in flow.iflat_tasks(nids=nids):
        if task.status not in (flow.S_OK, flow.S_RUN): continue
        if not task.is_abinit_task: continue

        abo_paths.append(task.output_file.path)
        index.append(task.pos_str)
        status.append(task.status.colored if with_colors else str(task.status))
        abo_relpaths.append(os.path.relpath(task.output_file.relpath))
        task_classes.append(task.__class__.__name__)
        task_nids.append(task.node_id)

    if not abo_paths: return

    # Get dimensions from output files as well as walltime/cputime
    robot = abilab.AboRobot.from_files(abo_paths)
    df = robot.get_dims_dataframe(with_time=True, index=index)

    # Add columns to the dataframe.
    status = [str(s) for s in status]
    df["task_class"] = task_classes
    df["relpath"] = abo_relpaths
    df["node_id"] = task_nids
    df["status"] = status

    if printout:
        abilab.print_dataframe(df, title="Table with Abinit dimensions:\n")

    return df


def flow_compare_abivars(flow, varnames, nids=None, wslice=None, printout=False, with_colors=False):
    """
    Print the input of the tasks to the given stream.

    Args:
        varnames:
            List of Abinit variables. If not None, only the variable in varnames
            are selected and printed.
        nids:
            List of node identifiers. By defaults all nodes are shown
        wslice:
            Slice object used to select works.
        printout: True to print dataframe.
        with_colors: True if task status should be colored.
    """
    varnames = [s.strip() for s in list_strings(varnames)]
    index, rows = [], []
    for task in flow.select_tasks(nids=nids, wslice=wslice):
        index.append(task.pos_str)
        dstruct = task.input.structure.as_dict(fmt="abivars")

        od = OrderedDict()
        for vname in varnames:
            value = task.input.get(vname, None)
            if value is None: # maybe in structure?
                value = dstruct.get(vname, None)
            od[vname] = value

        od["task_class"] = task.__class__.__name__
        od["status"] = task.status.colored if with_colors else str(task.status)
        rows.append(od)

    import pandas as pd
    df = pd.DataFrame(rows, index=index)
    if printout:
        abilab.print_dataframe(df, title="Input variables:")

    return df


def flow_debug_reset_tasks(flow, nids=None, verbose=0):
    """
    Analyze error files produced by reset tasks for possible error messages

    Args:
        nids: List of node identifiers. By defaults all nodes that have been resetted are analyzed.
        verbose: Verbosity level.
    """
    # TODO: Improve implementation.
    ntasks = 0
    nrows, ncols = get_terminal_size()
    for task in flow.select_tasks(nids=nids):
        # See task.reset_from_scratch
        if task.num_restarts == 0: continue
        reset_dir = os.path.join(task.workdir, "_reset")
        reset_file = os.path.join(reset_dir, "_counter")
        if not os.path.isdir(reset_dir) and not os.path.isfile(reset_file):
            continue

        ntasks += 1
        with open(reset_file, "rt") as fh:
            num_reset = int(fh.read())

        for i in range(num_reset):
            #("output_file", "log_file", "stderr_file", "qout_file", "qerr_file", "mpiabort_file")
            for fname in ("stderr_file", "qerr_file", "mpiabort_file"):
                path = os.path.join(reset_dir, fname + "_" + str(i))
                with open(path, "rt") as fh:
                    s = fh.read()
                    if not s: continue
                    print(2 * "\n")
                    print(make_banner(os.path.relpath(path), width=ncols, mark="="))
                    cprint(s, color="red")
                    print(2 * "\n")

    print("Number of tasks analyzed: %d" % ntasks)


def flow_watch_status(flow, delay=5, nids=None, verbose=0, func_name="show_func"):
    """
    Enter an infinite loop and delay execution for the given number of seconds. (default: 5 secs).

    Args:
        delay: delay execution for the given number of seconds. (default: 5 secs).
        nids: List of node identifiers. By defaults all nodes that have been resetted are analyzed.
        verbose: Verbosity level.
        func_name: Name of the function used to show the status of the flow.
    """
    cprint("Entering infinite loop (delay: %d s). Only changes are shown\nPress <CTRL+C> to exit" %
           delay, color="magenta", end="", flush=True)

    show_func = getattr(flow, func_name)
    assert callable(show_func)

    # Total counter and dicts used to detect changes.
    tot_count = 0
    before_task2stat, now_task2stat = {}, {}
    # Progressbar setup
    from tqdm import tqdm
    pbar, pbar_count, pbar_total = None, 0, 100

    exit_code = 0
    def exit_now():
        """
        Function used to test if we have to exit from the infinite loop below.
        Return: != 0 if we must exit. > 0 if some error occurred.
        """
        if flow.all_ok:
            cprint("Flow reached all_ok", "green")
            return -1
        if any(st.is_critical for st in before_task2stat.values()):
            cprint(boxed("Found tasks with critical status"), "red")
            return 1
        return 0

    try:
        while True:
            tot_count += 1
            flow.check_status()

            # Here I test whether there's been some change in the flow
            # before printing the status table.
            # Note that the flow in memory could not correspond to the one that
            # is being executed by the scheduler. This is the reason why we
            # reload it when we reach pbar_count.
            if tot_count == 1:
                for task in flow.iflat_tasks(nids=nids):
                    before_task2stat[task] = task.status
            else:
                for task in flow.iflat_tasks(nids=nids):
                    now_task2stat[task] = task.status

                if (len(before_task2stat) == len(now_task2stat) and
                    all(now_task2stat[t] == before_task2stat[t] for t in now_task2stat)):
                    # In principle this is not needed but ...
                    exit_code = exit_now()
                    if exit_code: break

                    # Progress bar section.
                    if pbar is None:
                        print("No change detected in the flow. Won't print status table till next change...")
                        pbar = tqdm(total=pbar_total)

                    if pbar_count <= pbar_total:
                        pbar_count += 1
                        pbar.update(1)
                    else:
                        pbar_count = 0
                        pbar.close()
                        pbar = tqdm(total=pbar_total)
                        flow.reload()

                    time.sleep(delay)
                    continue

                # copy now --> before
                before_task2stat = now_task2stat.copy()

            # Print status table. Exit if success or critical errors.
            print(2*"\n" + time.asctime() + "\n")
            show_func(verbose=verbose, nids=nids)
            # Add summary table to status table.
            if show_func is flow.show_status: flow.show_summary()

            exit_code = exit_now()
            if exit_code: break
            time.sleep(delay)

        # Print status table if something bad happened.
        if exit_code == 1:
            flow.show_status()

    except KeyboardInterrupt:
        cprint("Received KeyboardInterrupt from user\n", "yellow")

def get_epilog():
    usage = """\

Usage example:

###########
# Execution
###########

  abirun.py FLOWDIR rapid                 => Keep repeating, stop when no task can be executed.
  abirun.py FLOWDIR scheduler             => Execute flow with the scheduler.
  abirun.py FLOWDIR status                => Show status table.
  abirun.py FLOWDIR events                => Print ABINIT events (Warnings/Errors/Comments) found in log files.
  abirun.py FLOWDIR history               => Print Task histories.
  abirun.py FLOWDIR cancel                => Cancel jobs in the queue.
  abirun.py FLOWDIR debug                 => Analyze error files and log files for possible error messages.
  abirun.py FLOWDIR corrections           => Show AbiPy corrections performed at runtime.
  abirun.py FLOWDIR handlers              => Show event handlers installed in the flow.

##########
# Analysis
##########

  abirun.py FLOWDIR inputs                => Print input files.
  abirun.py FLOWDIR abivars -vn ecut,nband  => Print table with these input variables.
  abirun.py FLOWDIR structures            => Compare input/output structures of the tasks.
  abirun.py FLOWDIR ebands -t NscfTask    => Print table with electronic properties computed in NscfTask
  abirun.py FLOWDIR hist -p               => Print table with last iteration in hist files and plot results.
  abirun.py FLOWDIR cycles -p             => Print (and plot) SCF/relaxation cycles extracted from the output of the tasks.
  abirun.py FLOWDIR robot hist            => Build robot for HIST files, start ipython shell to interact with the robot.
  abirun.py FLOWDIR dims                  => Print table with dimensions extracted from the output of the tasks.
  abirun.py FLOWDIR group -g task_class   => Print table with node ids and tasks grouped by class
  abirun.py FLOWDIR inspect               => Call matplotlib to inspect the tasks.
  abirun.py FLOWDIR tail                  => Use Unix tail to follow the main output files of the flow.
  abirun.py FLOWDIR deps                  => Show task dependencies.
  abirun.py FLOWDIR listext GSR SIGRES    => Show output files with the given extension(s).
  abirun.py FILE abiopen                  => Open FILE with abiopen. Same interface as abiopen.py

###############
# Miscelleanous
###############

  abirun.py FLOWDIR ipython               => Open flow in ipython terminal.
  abirun.py FLOWDIR notebook              => Generate jupyter notebook.
  abirun.py FLOWDIR networkx              => Plot flow graph with networkx.
  abirun.py FLOWDIR graphviz              => Plot flow graph with graphviz
                                            (can also select tasks/works and show files with -d).
  abirun.py abibuild                      => Show ABINIT build information and exit.

###############
# Documentation
###############

  abirun.py FLOWDIR doc_manager slurm     => Document the TaskManager options availabe for Slurm.
  abirun.py . doc_manager script          => Show the job script that will be produced with the current settings.
  abirun.py . doc_scheduler               => Document the options available in scheduler.yml.
"""

    notes = """\

Notes:

    If FLOWDIR is not given, abirun.py automatically selects the database located within
    the working directory. An Exception is raised if multiple databases are found.

    Note, moreover, that one can replace FLOWDIR with the directory of a work/task
    to make the command operate on this node of the flow without having to specify the node ids with --nids.
    For example, to have the list of events of the task in `FLOWDIR/w0/t1` use:

        abirun.py FLOWDIR/w0/t1 events

    instead of:

        abirun.py FLOWDIR events -n 123

    where 123 is the node identifier associated to w0/t1.

    To start the scheduler with a time interval of 30 seconds, use:

        nohup abirun.py FLOWDIR scheduler -s 30 > sched.log 2> sched.err

    Alternatively one can specify the scheduler options via the `scheduler.yml` file.
    Remember that AbiPy will first look for `scheduler.yml` and `manager.yml` files
    in the current working directory and then inside $HOME/.abinit/abipy/

    Use `abirun.py --help` for help and `abirun.py COMMAND --help` to get the documentation for `COMMAND`.
    Use `-v` to increase verbosity level (can be supplied multiple times e.g -vv).
"""

    developers = """\

############
# Developers
############

    abirun.py prof ABIRUN_ARGS               => to profile abirun.py
    abirun.py tracemalloc ABIRUN_ARGS        => to trace memory blocks allocated by Python"""

    return notes + usage # + developers


def get_parser(with_epilog=False):

    def parse_nids(s):
        """parse nids argument"""
        if s is None: return s
        try:
            if "," in s:
                return [int(t) for t in s.split(",")]
            else:
                # Convert string to slice and return list.
                s = duck.as_slice(s)
                if s.stop is None: raise argparse.ArgumentTypeError("stop must be specified")
                return list(range(s.start, s.stop, s.step))
        except Exception:
            raise argparse.ArgumentTypeError(
                    "Invalid nids string %s\n Expecting None or int or comma-separated integers or slice sintax" % s)

    def parse_wslice(s):
        s = duck.as_slice(s)
        if s is None: return s
        if s.stop is None: raise argparse.ArgumentTypeError("stop must be specified")
        return s

    # Parent parser for commands that need to know on which subset of tasks/workflows we have to operate.
    # wslide and nids are mutually exclusive.
    flow_selector_parser = argparse.ArgumentParser(add_help=False)
    group = flow_selector_parser.add_mutually_exclusive_group()
    group.add_argument("-n", '--nids', default=None, type=parse_nids, help=(
        "Node identifier(s) used to select the task. Accept single integer, comma-separated list of integers or python slice.\n"
        "Use `status` command to get the node ids.\n"
        "Examples: --nids=12 --nids=12,13,16 --nids=10:12 to select 10 and 11 (slice syntax), --nids=2:5:2 to select 2,4."
        ))

    group.add_argument("-w", '--wslice', default=None, type=parse_wslice,
        help=("Select the list of works to analyze (python syntax for slices): "
              "Examples: --wslice=1 to select the second workflow, --wslice=:3 for 0,1,2, "
              "--wslice=-1 for the last workflow, --wslice::2 for even indices."))
    group.add_argument("-S", '--task-status', default=None, type=Status.as_status,
        help="Select only the tasks with the given status. Default: None i.e. ignored. Possible values: %s." %
             Status.all_status_strings())
    group.add_argument("-t", '--task-class', type=str, default=None,
        help="Select only tasks with the given class e.g. `-t NscfTask`.")

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity.')

    copts_parser.add_argument('--no-colors', default=False, action="store_true", help='Disable ASCII colors.')
    copts_parser.add_argument('--no-logo', default=False, action="store_true", help='Disable AbiPy logo.')
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG.")
    copts_parser.add_argument('--remove-lock', default=False, action="store_true",
        help="Remove the lock on the pickle file used to save the status of the flow.")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('flowdir', nargs="?", help=("File or directory containing the ABINIT flow/work/task. "
                                                    "If not given, the flow in the current workdir is selected."))
    parser.add_argument('-V', '--version', action='version', version=abilab.__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_single = subparsers.add_parser('single', parents=[copts_parser], help="Run single task and exit.")

    # Subparser for rapid command.
    p_rapid = subparsers.add_parser('rapid', parents=[copts_parser], help="Run all tasks in rapidfire mode.")
    p_rapid.add_argument('-m', '--max-nlaunch', default=10, type=int,
        help="Maximum number of launches. default: 10. Use -1 for no limit.")

    # Subparser for scheduler command.
    p_scheduler = subparsers.add_parser('scheduler', parents=[copts_parser],
        help="Run all tasks with a Python scheduler. Requires scheduler.yml either in $PWD or ~/.abinit/abipy.")
    p_scheduler.add_argument('-w', '--weeks', default=0, type=int, help="Number of weeks to wait.")
    p_scheduler.add_argument('-d', '--days', default=0, type=int, help="Number of days to wait.")
    p_scheduler.add_argument('-hs', '--hours', default=0, type=int, help="Number of hours to wait.")
    p_scheduler.add_argument('-m', '--minutes', default=0, type=int, help="Number of minutes to wait.")
    p_scheduler.add_argument('-s', '--seconds', default=0, type=int, help="Number of seconds to wait.")

    # Subparser for batch command.
    p_batch = subparsers.add_parser('batch', parents=[copts_parser], help="Run scheduler in batch script.")
    p_batch.add_argument("-t", '--timelimit', default=None, help=("Time limit for batch script. "
                         "Accept int with seconds or string with time given in the slurm convention: "
                         "`days-hours:minutes:seconds`. If timelimit is None, the default value specified"
                         " in the `batch_adapter` entry of `manager.yml` is used."))

    # Subparser for status command.
    p_status = subparsers.add_parser('status', parents=[copts_parser, flow_selector_parser], help="Show status table.")
    p_status.add_argument('-d', '--delay', nargs="?", const=5, default=0, type=int,
        help="Enter an infinite loop and delay execution for the given number of seconds. (default: 5 secs).")
    p_status.add_argument('-s', '--summary', default=False, action="store_true",
        help="Print short version with status counters.")

    # Subparser for set_status command.
    p_set_status = subparsers.add_parser('set_status', parents=[copts_parser, flow_selector_parser],
        help="Change the status of the task. WARNING: Option for developers!")
    p_set_status.add_argument('new_status', help="New value of status. Possible values: %s." % Status.all_status_strings())

    # Subparser for cancel command.
    p_cancel = subparsers.add_parser('cancel', parents=[copts_parser, flow_selector_parser],
        help="Cancel the tasks in the queue. Not available if qtype is shell.")
    p_cancel.add_argument("-r", "--rmtree", action="store_true", default=False, help="Remove flow directory.")

    # Subparser for restart command.
    p_restart = subparsers.add_parser('restart', parents=[copts_parser, flow_selector_parser],
        help=("Restart the tasks of the flow. By default, only the task whose status==Unconverged are restarted. "
              "Use -S `status` and/or -n node_ids to select particular tasks."))

    # Subparser for reset command.
    p_reset = subparsers.add_parser('reset', parents=[copts_parser, flow_selector_parser],
        help="Reset the tasks of the flow with the specified status.")

    # Subparser for move command.
    p_move = subparsers.add_parser('move', parents=[copts_parser],
        help="Move the flow to a new directory and change the absolute paths.")
    p_move.add_argument('dest', nargs=1)

    # Subparser for open command.
    p_open = subparsers.add_parser('open', parents=[copts_parser, flow_selector_parser],
        help="Open files in $EDITOR, type `abirun.py FLOWDIR open --help` for help).")
    p_open.add_argument('what', nargs="?", default="o",
        choices=["i", "o", "f", "j", "l", "e", "q", "all"],
        help=r"""\
Specify the files to open. Possible choices:
    i ==> input_file
    o ==> output_file
    f ==> files_file
    j ==> job_file
    l ==> log_file
    e ==> stderr_file
    q ==> qout_file
    all ==> all files.
Default: o
""")

    # Subparser for abiopen.
    p_abiopen = subparsers.add_parser('abiopen', parents=[copts_parser], #, flow_selector_parser],
        help="Open FILE with abiopen. Same interface as abiopen.py")
    p_abiopen.add_argument('-nb', '--notebook', action='store_true', default=False, help="Open file in jupyter notebook")
    p_abiopen.add_argument('--foreground', action='store_true', default=False,
                        help="Run jupyter notebook in the foreground.")
    p_abiopen.add_argument('-p', '--print', action='store_true', default=False, help="Print python object and return.")

    # Subparser for abibuild
    p_abibuild = subparsers.add_parser('abibuild', parents=[copts_parser, flow_selector_parser],
        help="Show ABINIT build information and exit.")

    # Subparser for doc_scheduler
    p_docsched = subparsers.add_parser('doc_scheduler', parents=[copts_parser],
        help="Document the options available in scheduler.yml.")

    # Subparser for gui command.
    #p_gui = subparsers.add_parser('gui', parents=[copts_parser], help="Open the GUI (requires wxPython).")
    #p_gui.add_argument("--chroot", default="", type=str, help=("Use chroot as new directory of the flow. " +
    #                   "Mainly used for opening a flow located on a remote filesystem mounted with sshfs. " +
    #                   "In this case chroot is the absolute path to the flow on the **localhost** " +
    #                   "Note that it is not possible to change the flow from remote when chroot is used."))

    # Subparser for new_manager.
    p_new_manager = subparsers.add_parser('new_manager', parents=[copts_parser, flow_selector_parser],
        help="Change the TaskManager.")
    p_new_manager.add_argument("manager_file", default="", type=str, help="YAML file with the new manager.")

    # Subparser for tail.
    p_tail = subparsers.add_parser('tail', parents=[copts_parser, flow_selector_parser],
        help="Use unix tail to follow the main output files of the flow.")
    p_tail.add_argument('what_tail', nargs="?", type=str, default="o", choices=["o", "l", "e"],
        help="What to follow: `o` for output (default), `l` for logfile, `e` for stderr.")

    # Subparser for qstat.
    # TODO: finalize the implementation
    #p_qstat = subparsers.add_parser('qstat', parents=[copts_parser], help="Show additional info on the jobs in the queue.")

    # Subparser for deps.
    p_deps = subparsers.add_parser('deps', parents=[copts_parser], help="Show dependencies.")

    # Subparser for robot.
    p_robot = subparsers.add_parser('robot', parents=[copts_parser, flow_selector_parser],
        help="Use a robot to analyze the results of multiple tasks (requires ipython).")
    p_robot.add_argument('robot_ext', nargs="?", type=str, default="GSR",
        help=("The file extension of the output file (case insensitive). "
              "Must be in: %s" % str(abilab.Robot.get_supported_extensions())))
    p_robot.add_argument('-nb', '--notebook', action='store_true', default=False, help="Generate jupyter notebook.")
    p_robot.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")
    p_robot.add_argument('-p', '--print-dataframe', action='store_true', default=False,
        help="Invoke robot.get_dataframe() and printa table to terminal.")

    # Subparser for cycles.
    p_cycles = subparsers.add_parser('cycles', parents=[copts_parser, flow_selector_parser],
        help=("Print self-consistent/relaxation cycles extracted from the output of the tasks."
              "Use `-t ScfTask` or `FLOWDIR/w0/` or `FLOWDIR/w0/t1/` to select tasks by class, by work or by workdir."))
    p_cycles.add_argument("-wok", "--exclude-ok-tasks", action='store_true', default=False,
        help="Exclude Tasks that have reached S_OK.")
    p_cycles.add_argument("-p", "--plot-mode", nargs="?", default=None, const="combiplot",
        choices=["combiplot", "slideshow"], #"gridplot",
        help=("Plot multiple cycles on the same figure if arg is specified. Use `-p` for gridplot. "
              "Use `-p slideshow` to iterate."))

    # Subparser for dims.
    p_dims = subparsers.add_parser('dims', parents=[copts_parser, flow_selector_parser],
        help="Print table with dimensions extracted from the output of the tasks.")

    # Subparser for inspect.
    p_inspect = subparsers.add_parser('inspect', parents=[copts_parser, flow_selector_parser],
        help="Call matplotlib to inspect the tasks (execute task.inspect method)")

    # Subparser for inputs.
    p_inputs = subparsers.add_parser('inputs', parents=[copts_parser, flow_selector_parser],
        help="Show the input files of the tasks.")
    p_inputs.add_argument("-vn", "--varnames", nargs="?", default=None, type=parse_strings,
        help="Comma-separated variable names. Can be used to print only these variables.")

    # Subparser for abivars.
    p_abivars = subparsers.add_parser('abivars', parents=[copts_parser, flow_selector_parser],
        help="Show pandas dataframe with Abinit input variables.")
    p_abivars.add_argument("-vn", "--varnames", required=True, type=parse_strings,
        help="Comma-separated variable names e.g. `-vn ecut,nband,ngkpt`.")

    # Subparser for structures command.
    p_structures = subparsers.add_parser('structures', parents=[copts_parser, flow_selector_parser],
        help="Compare input/output structures of the tasks. Print max force and pressure if available.")
    p_structures.add_argument("--what", type=str, default="io",
        help="'i' for input structures, 'o' for output, 'io' for both.")

    # Subparser for ebands command.
    p_ebands = subparsers.add_parser('ebands', parents=[copts_parser, flow_selector_parser],
        help="Compare electronic bands produced by the tasks.")
    p_ebands.add_argument("-p", "--plot-mode", nargs="?", default=None, const="gridplot",
        choices=["gridplot", "combiplot", "boxplot", "combiboxplot", "animate"],
        help="Plot multiple bands if arg is specified. Use -p for gridplot. Supports multiple formats.")

    # Subparser for hist command.
    p_hist = subparsers.add_parser('hist', parents=[copts_parser, flow_selector_parser],
        help="Compare HIST.nc files produced by the tasks.")
    p_hist.add_argument("-p", "--plot-mode", nargs="?", default=None, const="combiplot",
        choices=["gridplot", "combiplot"],
        help="Plot multiple HIST files if arg is present. Use -p for combiplot. Supports multiple formats.")

    # Subparser for manager.
    p_manager = subparsers.add_parser('doc_manager', parents=[copts_parser], help="Document the TaskManager options.")
    p_manager.add_argument("qtype", nargs="?", default=None, help=("Write job script to terminal if qtype='script' else "
        "document the qparams for the given QueueAdapter qtype e.g. slurm."))

    # Subparser for events.
    p_events = subparsers.add_parser('events', parents=[copts_parser, flow_selector_parser],
        help="Show ABINIT events (error messages, warnings, comments).")
    #p_events.add_argument("-t", "event-type", default=)

    # Subparser for corrections.
    p_corrections = subparsers.add_parser('corrections', parents=[copts_parser, flow_selector_parser],
        help="Show AbiPy corrections performed at runtime.")

    # Subparser for history.
    p_history = subparsers.add_parser('history', parents=[copts_parser, flow_selector_parser], help="Show Node history.")
    p_history.add_argument("-m", "--metadata", action="store_true", default=False, help="Print history metadata.")
    p_history.add_argument("-f", "--full-history", action="store_true", default=False,
        help="Print full history set, including nodes with an empty history.")

    # Subparser for handlers.
    p_handlers = subparsers.add_parser('handlers', parents=[copts_parser],
        help="Show event handlers installed in the flow.")
    p_handlers.add_argument("-d", "--doc", action="store_true", default=False,
        help="Show documentation about all the handlers that can be installed.")

    # Subparser for notebook.
    p_notebook = subparsers.add_parser('notebook', parents=[copts_parser],
        help="Create and open an ipython notebook to interact with the flow.")
    p_notebook.add_argument('--foreground', action='store_true', default=False,
        help="Run jupyter notebook in the foreground.")

    # Subparser for ipython.
    p_ipython = subparsers.add_parser('ipython', parents=[copts_parser],
        help="Embed IPython. Useful for advanced operations or debugging purposes.")
    p_ipython.add_argument('--argv', nargs="?", default="", type=shlex.split,
        help="Command-line options passed to ipython. Must be enclosed by quotes. "
             "Example: --argv='--matplotlib=wx'")

    # Subparser for tar.
    p_tar = subparsers.add_parser('tar', parents=[copts_parser], help="Create tarball file.")
    p_tar.add_argument("-s", "--max-filesize", default=None,
        help="Exclude file whose size > max-filesize bytes. Accept integer or string e.g `1Mb`.")
    p_tar.add_argument("-e", "--exclude-exts", default=None, type=parse_strings,
        help="Exclude file extensions. Accept string or comma-separated strings. Ex: -eWFK or --exclude-exts=WFK,GSR")
    p_tar.add_argument("-d", "--exclude-dirs", default=None, type=parse_strings,
        help="Exclude directories. Accept string or comma-separated strings. Ex: --exlude-dirs=indir,outdir")
    p_tar.add_argument("-l", "--light", default=False, action="store_true",
        help="Create light-weight version of the tarball for debugging purposes. Other options are ignored.")

    # Subparser for tricky.
    p_tricky = subparsers.add_parser('tricky', parents=[copts_parser],
        help=("Show tricky tasks i.e. tasks that have been restarted, "
              "launched more than once or tasks that have been corrected."))

    # Subparser for debug.
    p_debug = subparsers.add_parser('debug', parents=[copts_parser, flow_selector_parser],
        help="Analyze error files and log files for possible error messages.")

    # Subparser for debug_reset.
    p_debug_reset = subparsers.add_parser('debug_reset', parents=[copts_parser, flow_selector_parser],
        help="Analyze error files and log files produced by reset tasks for possible error messages.")

    # Subparser for clone_task.
    p_clone_task = subparsers.add_parser('clone_task', parents=[copts_parser, flow_selector_parser],
        help="Clone task, change input variables and add new tasks to the flow. Requires clone_task.py.")

    # Subparser for group.
    p_group = subparsers.add_parser('group', parents=[copts_parser, flow_selector_parser],
        help="Group tasks according to property.")
    p_group.add_argument("-g", '--groupby', default="status", type=str, choices=["status", "task_class"],
        help="Select the attribute used to group tasks. Default: status.")

    # Subparser for diff.
    p_diff = subparsers.add_parser('diff', parents=[copts_parser, flow_selector_parser],
        help="Compare files produced by two or three nodes.")
    p_diff.add_argument('what_diff', nargs="?", type=str, default="i",
        help="What to diff: `i` for input (default), `o` for output, `l` for logfile, `e` for stderr.")

    # Subparser for networkx.
    p_networkx = subparsers.add_parser('networkx', parents=[copts_parser],
        help="Draw flow and node dependencies with networkx package.")
    p_networkx.add_argument('--nxmode', default="status",
        help="Type of network plot. Possible values: `status`, `network`. Default: `status`.")
    p_networkx.add_argument('--edge-labels', action="store_true", default=False, help="Show edge labels.")

    # Subparser for graphviz.
    p_graphviz = subparsers.add_parser('graphviz', parents=[copts_parser],
        help=("Draw flow and node dependencies with graphviz package. Accept (FLOWDIR|WORKDIR|TASKDIR). "
             "See https://graphviz.readthedocs.io/."))
    p_graphviz.add_argument("-e", "--engine", type=str, default="automatic",
        help=("graphviz engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']. "
            "Default: automatic i.e. the engine is automatically selected. See http://www.graphviz.org/pdf/dot.1.pdf "
            "Use `conda install python-graphviz` or `pip install graphviz` to install the python package"))
    p_graphviz.add_argument("-d", '--dirtree', default=False, action="store_true",
        help='Visualize files and directories in workdir instead of tasks/works.')

    # Subparser for listext.
    p_listext = subparsers.add_parser('listext', parents=[copts_parser],
        help="List all the output files with the given extension that have been produced by the nodes.")
    p_listext.add_argument('listexts', nargs="*", default=[], help="List of Abinit file extensions. e.g DDB, GSR, WFK etc")

    # Subparser for timer.
    # TODO
    #p_timer = subparsers.add_parser('timer', parents=[copts_parser, flow_selector_parser],
    #    help=("Read the section with timing info from the main ABINIT output file (requires timopt != 0) "
    #          "Open Ipython terminal to inspect data."))

    return parser


@prof_main
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

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if options.verbose > 2:
        print(options)

    # Documentation options that do not need a flow.
    # Print docs and exit immediately.
    if options.command == "doc_manager":
        # Document TaskManager options and qparams.
        qtype = options.qtype

        if qtype == "script":
            manager = flowtk.TaskManager.from_user_config()
            script = manager.qadapter.get_script_str(
                job_name="job_name", launch_dir="workdir", executable="executable",
                qout_path="qout_file.path", qerr_path="qerr_file.path",
                stdin="stdin", stdout="stdout", stderr="stderr")
            print(script)

        else:
            print(flowtk.TaskManager.autodoc())
            print("qtype supported: %s" % flowtk.all_qtypes())
            print("Use `abirun.py . manager slurm` to have the list of qparams for slurm.\n")

            if qtype is not None:
                print("QPARAMS for %s" % qtype)
                flowtk.show_qparams(qtype)

        return 0

    if options.command == "doc_scheduler":
        print("Options that can be specified in scheduler.yml:")
        print(flowtk.PyFlowScheduler.autodoc())
        return 0

    if options.command == "abibuild":
        abinit_build = flowtk.AbinitBuild()
        print()
        print(abinit_build)
        print()
        if not options.verbose:
            print("Use --verbose for additional info")
        else:
            print(abinit_build.info)
        return 0

    # Abiopen does not need flow
    if options.command == "abiopen":
        return cli_abiopen(options, options.flowdir)

    # After this point we start to operate on the flow.
    # 0) Print logo
    # 1) Read flow from pickle file and construct nids set if needed.
    # 2) Operate on the flow depending on the options specified by the user on the CLI.
    if options.no_colors:
        # Disable colors
        termcolor.enable(False)

    if not options.no_logo:
        nrows, ncols = get_terminal_size()
        if ncols > 100: cprint(abilab.abipy_logo1(), "yellow")

        system, node, release, version, machine, processor = platform.uname()
        cprint("Running on %s -- system %s -- Python %s -- %s" % (
               gethostname(), system, platform.python_version(), "abirun" + "-" + abilab.__version__),
               'yellow', attrs=['underline'])

    wname, tname = None, None
    if options.flowdir is None:
        # Will try to figure out the location of the Flow.
        options.flowdir = os.getcwd()
    else:
        # Sometimes one wants to inspect a work or a task by just using `abirun.py flow/w0/t0 inspect`
        # without knowing its node id. flowdir_wname_tname will solve the problem!
        options.flowdir, wname, tname = flowdir_wname_tname(options.flowdir)

    # Read the flow from the pickle database.
    flow = flowtk.Flow.pickle_load(options.flowdir, remove_lock=options.remove_lock)
    #flow.show_info()

    # If we have selected a work/task, we have to convert wname/tname into node ids (nids)
    if wname or tname:
        if wname and tname:
            # Task
            for w_pos, work in enumerate(flow):
                if os.path.basename(work.workdir) == wname: break
            else:
                raise RuntimeError("Cannot find work from name %s" % wname)

            for t_pos, task in enumerate(flow[w_pos]):
                if os.path.basename(task.workdir) == tname: break
            else:
                raise RuntimeError("Cannot find task from name %s" % tname)

            # Create options.nids here
            options.nids = set([flow[w_pos].node_id, flow[w_pos][t_pos].node_id])

        else:
            # Work
            for w_pos, work in enumerate(flow):
                if os.path.basename(work.workdir) == wname: break
            else:
                raise RuntimeError("Cannot find work from name %s" % wname)

            # Create options.nids here
            options.nids = set([flow[w_pos].node_id] + [task.node_id for task in flow[w_pos]])

    retcode = 0

    #if options.command == "gui":
    #    if options.chroot:
    #        # Change the workdir of flow.
    #        print("Will chroot to %s..." % options.chroot)
    #        flow.chroot(options.chroot)
    #    from abipy.gui.flowviewer import wxapp_flow_viewer
    #    wxapp_flow_viewer(flow).MainLoop()

    if options.command == "new_manager":
        # Read the new manager from file.
        new_manager = flowtk.TaskManager.from_file(options.manager_file)

        # Default status for new_manager is QCritical
        if options.task_status is None:
            options.task_status = Status.as_status("QCritical")

        # Change the manager of the errored tasks.
        print("Resetting tasks with status: %s" % options.task_status)
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            task.reset()
            task.set_manager(new_manager)

        # Update the database.
        return flow.build_and_pickle_dump()

    elif options.command == "events":
        flow.show_events(status=options.task_status, nids=selected_nids(flow, options))

    elif options.command == "corrections":
        flow.show_corrections(status=options.task_status, nids=selected_nids(flow, options))

    elif options.command == "history":
        flow.show_history(status=options.task_status, nids=selected_nids(flow, options),
                          full_history=options.full_history, metadata=options.metadata)

    elif options.command == "handlers":
        if options.doc:
            flowtk.autodoc_event_handlers()
        else:
            flow.show_event_handlers(verbose=options.verbose)

    elif options.command  == "single":
        nlaunch = flow.single_shot()
        if nlaunch: flow.show_status()
        cprint("Number of tasks launched: %d" % nlaunch, "yellow")

    elif options.command == "rapid":
        nlaunch = flow.rapidfire(max_nlaunch=options.max_nlaunch, max_loops=1, sleep_time=5)
        if nlaunch: flow.show_status()
        cprint("Number of tasks launched: %d" % nlaunch, "yellow")

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
            flow_watch_status(flow, delay=options.delay, verbose=options.verbose,
                              nids=selected_nids(flow, options), func_name=show_func.__name__)
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

    elif options.command == "cancel":
        print("Number of jobs cancelled %d" % flow.cancel(nids=selected_nids(flow, options)))
        # Remove directory
        if options.rmtree: flow.rmtree()

    elif options.command == "restart":
        # Default status for restart is Unconverged if no option is provided by the user.
        if options.task_status is None and options.nids is None:
            options.task_status = Status.as_status("Unconverged")

        nlaunch, excs = 0, []
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            #if options.verbose:
            print("Will try to restart %s, with status %s" % (task, task.status))
            try:
                fired = task.restart()
                if fired:
                    nlaunch += 1
                    print("\tTask restarted")
            except Exception:
                excs.append(straceback())

        cprint("Total number of jobs restarted %d" % nlaunch, "blue")
        if nlaunch:
            # update database
            flow.pickle_dump()

        if excs:
            print("Exceptions raised\n")
            pprint(excs)

    elif options.command == "reset":
        if options.nids is None and options.task_status is None:
            # Default status for reset command is QCritical
            options.task_status = Status.as_status("QCritical")

        if options.task_status is not None:
            print("Resetting tasks with status: %s" % options.task_status)
        else:
            print("Resetting tasks with node ids: %s" % str(options.nids))

        count = 0
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            print("Resetting task %s... " % task, end="")
            failed = task.reset()
            if failed:
                cprint("[FAILED]", "red")
            else:
                cprint("[OK]", "green")
                count += 1
        cprint("%d tasks have been reset" % count, "blue")

        flow.show_status()

        g = flow.find_deadlocks()
        #print("deadlocked:", gdeadlocked, "\nrunnables:", grunnables, "\nrunning:", g.running)
        if g.deadlocked and not (g.runnables or g.running):
            cprint("*** Flow is deadlocked ***", "red")

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
            cprint("No job is running. Exiting!", "magenta")
        else:
            cprint("Press <CTRL+C> to interrupt. Number of output files %d\n" % len(paths),
                    color="magenta", end="", flush=True)
            try:
                os.system("tail -f %s" % " ".join(paths))
            except KeyboardInterrupt:
                cprint("Received KeyboardInterrupt from user\n", "yellow")

    #elif options.command == "qstat":
    #    print("Warning: this option is still under development.")
    #    #for task in flow.select_tasks(nids=options.nids, wslice=options.wslice):
    #    for task in flow.iflat_tasks():
    #        if not task.qjob: continue
    #        print("qjob", task.qjob)
    #        print("info", task.qjob.get_info())
    #        print("estimated_start-time", task.qjob.estimated_start_time())
    #        print("qstats", task.qjob.get_stats())

    elif options.command == "deps":
        flow.check_status()
        flow.show_dependencies()

    elif options.command == "robot":
        print("Building robot for file extension:", options.robot_ext, "with task_class:", options.task_class)
        robot = abilab.Robot.from_flow(flow, outdirs="all", nids=selected_nids(flow, options),
                                       ext=options.robot_ext)
        if len(robot) == 0:
            cprint("Robot couldn't find files", "yellow")
            return 1

        print("\n%s\n" % robot.get_label_files_str())

        if options.print_dataframe:
            # Print dataframe if robot provides get_dataframe method.
            if hasattr(robot, "get_dataframe"):
                try:
                    df = robot.get_dataframe()
                    abilab.print_dataframe(df, title="Output of robot.get_dataframe():", precision=4)
                except Exception as exc:
                    cprint("Exception:\n%s\n\nwhile invoking get_dataframe. Falling back to to_string" % str(exc), "red")
                    print(robot.to_string(verbose=options.verbose))
            else:
                cprint("%s does not provide `get_dataframe` method. Using `to_string`" % (
                        robot.__class__.__name__), "yellow")
                print(robot.to_string(verbose=options.verbose))

        elif options.notebook:
            # Create notebook.
            return robot.make_and_open_notebook(foreground=options.foreground)

        else:
            # Start ipython to interact with the robot.
            import IPython
            header = robot.get_label_files_str() + "\n\nType `robot` in the terminal and use <TAB> to list its methods."
            IPython.embed(header=header, robot=robot)

    elif options.command == "cycles":
        # Print cycles.
        from pymatgen.io.abinit.abiinspect import CyclesPlotter
        cls2plotter = OrderedDict()
        for task, cycle in flow.get_task_scfcycles(nids=selected_nids(flow, options),
                                                   exclude_ok_tasks=options.exclude_ok_tasks):
            print()
            cprint(repr(task), **task.status.color_opts)
            print("\n", cycle, "\n")

            label = repr(task) if options.verbose else "{} {}".format(task.__class__.__name__, task.relworkdir)
            # Could have different kind of cycles e.g. Scf, Relax, DFPT ...
            # so we group them by building multiple plotters indexed by task class.
            # Plots are produced outside of the loop.
            pkey = task.__class__.__name__
            if pkey not in cls2plotter: cls2plotter[pkey] = CyclesPlotter()
            cls2plotter[pkey].add_label_cycle(label, cycle)

        # Here comes the plot.
        if options.plot_mode is not None and cls2plotter:
            for pkey, plotter in cls2plotter.items():
                try:
                    getattr(plotter, options.plot_mode)(title="Plotter for %s class" % pkey)
                except Exception as exc:
                    cprint("Exception while invoking %s method of %s.\n%s" % (
                           options.plot_mode, plotter.__class__.__name__, str(exc)), "red")

    elif options.command == "dims":
        flow_get_dims_dataframe(flow, nids=selected_nids(flow, options),
                                printout=True, with_colors=not options.no_colors)

    elif options.command == "inspect":
        tasks = flow.select_tasks(nids=selected_nids(flow, options))
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

    elif options.command == "inputs":
        flow.show_inputs(varnames=options.varnames, nids=selected_nids(flow, options))

    elif options.command == "abivars":
        flow_compare_abivars(flow, varnames=options.varnames, nids=selected_nids(flow, options),
                             printout=True, with_colors=not options.no_colors)

    elif options.command == "structures":
        flow_compare_structures(flow, nids=selected_nids(flow, options), what=options.what,
                                verbose=options.verbose, with_spglib=False, printout=True,
                                with_colors=not options.no_colors)

    elif options.command == "ebands":
        flow_compare_ebands(flow, nids=selected_nids(flow, options), verbose=options.verbose,
                            with_spglib=False, printout=True, with_colors=not options.no_colors,
                            plot_mode=options.plot_mode)

    elif options.command == "hist":
        flow_compare_hist(flow, nids=selected_nids(flow, options), verbose=options.verbose,
                          with_spglib=False, printout=True, with_colors=not options.no_colors,
                          plot_mode=options.plot_mode)

    elif options.command == "notebook":
        return flow_write_open_notebook(flow, options)

    elif options.command == "ipython":
        import IPython
        print("Invoking Ipython, `flow` object will be available in the Ipython terminal")
        IPython.start_ipython(argv=options.argv, user_ns={"flow": flow})

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

    elif options.command == "tricky":
        flow.show_tricky_tasks(verbose=options.verbose)

    elif options.command == "debug":
        flow.debug(status=options.task_status, nids=selected_nids(flow, options))

    elif options.command == "debug_reset":
        flow_debug_reset_tasks(flow, nids=selected_nids(flow, options), verbose=options.verbose)

    # TODO
    #elif options.command == "debug_restart":
    #    flow_debug_restart_tasks(flow, nids=selected_nids(flow, options), verbose=options.verbose)

    elif options.command == "clone_task":
        if wname is None and tname is None:
            cprint("Use e.g. `abirun.py FLOWDIR/w0/t0` to select the task to clone.", "yellow")
            return 1

        if flow.has_scheduler:
            cprint("Cannot add new tasks when there's a scheduler running in background.", "yellow")
            return 1

        task_dirpath = os.path.join(flow.workdir, wname, tname)
        for task in flow.iflat_tasks():
            if task.workdir == task_dirpath:
                task_id = task.node_id
                print("Will clone task: ", repr(task))
                break
        else:
            raise ValueError("Cannot find task associated to workdir `%s`" % task_dirpath)

        #print(task.deps, type(task.deps))
        py_file = "clone_task.py"
        if not os.path.exists(py_file):
            cprint("clone_task requires %s in the current working directory" % py_file, "yellow")
            cprint("Will generate template file. Please edit it and rerun", "yellow")
            template = r"""
def list_of_dict_with_vars(task):
    "
    This function is called by `abirun.py clone_task` to build new tasks.
    It receives the task to be cloned and retunn a list of dictionaries.
    Each dictionary contains the Abinit variables that will be added to the initial input.
    To build e.g. two new tasks with a different value of ecut use:

    .. example:

        return [
            {"ecut": 20},
            {"ecut": 30},
        ]
    "
    #return [
    #    {"ecut": 20, nband: 10},
    #    {"ecut": 30, nband: 20},
    #]"""
            with open(py_file, "wt") as fh:
                fh.write(template)
            return 1
        else:
            print("Importing `list_of_dict_with_vars` from ", pyfile)
            import imp
            mod = imp.load_source(pyfile.replace(".py", ""), pyfile)

        dict_list = mod.list_of_dict_with_vars(task)
        if not dict_list:
            cprint("list_of_dict_with_vars returned empty list", "red")
            return 1

        for d in dict_list:
            print("Registering new task with vars:", d)
            task.work.register(task.input.new_with_vars(**d),
                               deps=task.deps, task_class=task.__class__)

        task.work.finalized = False
        flow.allocate()
        if ask_yesno("Do you want to rebuild the flow? [Y/n]"):
            flow.build_and_pickle_dump()

    elif options.command == "group":
        d = defaultdict(list)
        print("\nMapping `%s` ---> List of node identifiers" % options.groupby)
        for task in flow.iflat_tasks(status=options.task_status, nids=selected_nids(flow, options)):
            if options.groupby == "status":
                k = task.status
            elif options.groupby == "task_class":
                k = task.__class__.__name__
            else:
                raise ValueError("Invalid groupby: `%s`" % options.groupby)
            d[k].append(task)

        for k, tasks in d.items():
            s = ", ".join([colored(str(task.node_id), **task.status.color_opts) for task in tasks])
            print("   ",k, " ---> ", s)

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
        flow.plot_networkx(mode=options.nxmode, with_edge_labels=options.edge_labels)

    elif options.command == "graphviz":
        # Select node to visualize.
        node = flow
        if wname or tname:
            if wname and tname: # Task
                node = flow[w_pos][t_pos]
            else:  # Work
                node = flow[w_pos]

        #dpi = 300
        #node.graphviz_imshow(dpi=dpi)
        #node.graphviz_imshow(figsize=(800/dpi, 800/dpi), dpi=dpi)
        #return 0

        directory = tempfile.mkdtemp()
        print("Producing source files in:", directory)

        if options.dirtree:
            graph = node.get_graphviz_dirtree(engine=options.engine)
        else:
            graph = node.get_graphviz(engine=options.engine)

        graph.view(directory=directory, cleanup=False)

    elif options.command == "listext":
        if not options.listexts:
            print("\nPlease specify the file extension(s), e.g. GSR SIGRES.\nList of available extensions:\n")
            print(abilab.abiopen_ext2class_table())
            return 0

        for ext in options.listexts:
            print("")
            flow.listext(ext)
            print("")

    #elif options.command == "timer":
    #    print("Warning this option is still under development")
    #    timer = flow.parse_timing()
    #    if timer is None:
    #        cprint("Cannot parse timer data!", color="magenta", end="", flush=True)
    #        return 1

    #    import IPython
    #    IPython.start_ipython(argv=[], user_ns={"timer": timer})

    else:
        raise RuntimeError("Don't know what to do with command %s!" % options.command)

    return retcode

if __name__ == "__main__":
    sys.exit(main())
