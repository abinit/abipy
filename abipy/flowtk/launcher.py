# coding: utf-8
"""Tools for the submission of Tasks."""
from __future__ import annotations

import abc
import os
import time
import datetime
import pandas as pd
import apscheduler

from collections import deque
from io import StringIO
from queue import Queue, Empty
from typing import List, Optional
from monty.io import get_open_fds
from monty.string import boxed, is_string
from monty.os.path import which
from monty.collections import AttrDict #, dict2namedtuple
from monty.termcolor import cprint
from monty.functools import lazy_property
from abipy.tools.iotools import yaml_safe_load, ask_yesno
from .utils import as_bool


import logging
logger = logging.getLogger(__name__)

try:
    has_sched_v3 = apscheduler.version >= "3.0.0"
except AttributeError:
    has_sched_v3 = False


__all__ = [
    "ScriptEditor",
    "PyLauncher",
    "PyFlowScheduler",
    "MultiFlowScheduler",
]


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()



class ScriptEditor(object):
    """
    Simple editor to simplify the writing of shell scripts
    """

    _shell = '/bin/bash'

    def __init__(self):
        self._lines = []

    @property
    def shell(self):
        return self._shell

    def _add(self, text, pre=""):
        if is_string(text):
            self._lines.append(pre + text)
        else:
            self._lines.extend([pre + t for t in text])

    def reset(self):
        """Reset the editor."""
        try:
            del self._lines
        except AttributeError:
            pass

    def shebang(self):
        """Adds the shebang line."""
        self._lines.append('#!' + self.shell)

    def declare_var(self, key, val):
        """Declare a env variable. If val is None the variable is unset."""
        if val is not None:
            line = "export " + key + '=' + str(val)
        else:
            line = "unset " + key

        self._add(line)

    def declare_vars(self, d):
        """Declare the variables defined in the dictionary d."""
        for k, v in d.items():
            self.declare_var(k, v)

    def export_envar(self, key, val):
        """Export an environment variable."""
        line = "export " + key + "=" + str(val)
        self._add(line)

    def export_envars(self, env):
        """Export the environment variables contained in the dict env."""
        for k, v in env.items():
            self.export_envar(k, v)

    def add_emptyline(self):
        """Add an empty line."""
        self._add("", pre="")

    def add_comment(self, comment):
        """Add a comment"""
        self._add(comment, pre="# ")

    def load_modules(self, modules):
        """Load the list of specified modules."""
        for module in modules:
            self.load_module(module)

    def load_module(self, module):
        self._add('module load ' + module + " 2>> mods.err")

    def add_line(self, line):
        self._add(line)

    def add_lines(self, lines):
        self._add(lines)

    def get_script_str(self, reset=True):
        """Returns a string with the script and reset the editor if reset is True"""
        s = "\n".join(l for l in self._lines)
        if reset:
            self.reset()
        return s


class PyLauncherError(Exception):
    """Error class for PyLauncher."""


class PyLauncher(object):
    """
    This object handle the submission of the tasks contained in a |Flow|.
    """

    Error = PyLauncherError

    def __init__(self, flow, **kwargs):
        """
        Initialize the object

        Args:
            flow: |Flow| object
            max_njobs_inqueue: The launcher will stop submitting jobs when the
                number of jobs in the queue is >= Max number of jobs
        """
        self.flow = flow
        self.max_njobs_inqueue = kwargs.get("max_njobs_inqueue", 200)

        #self.flow.check_pid_file()

    def single_shot(self):
        """
        Run the first :class:`Task` than is ready for execution.

        Returns:
            Number of jobs launched.
        """
        num_launched = 0

        # Get the tasks that can be executed in each workflow.
        tasks = []
        for work in self.flow:
            try:
                task = work.fetch_task_to_run()

                if task is not None:
                    tasks.append(task)
                else:
                    # No task found, this usually happens when we have dependencies.
                    # Beware of possible deadlocks here!
                    logger.debug("No task to run! Possible deadlock")

            except StopIteration:
                logger.info("All tasks completed.")

        # Submit the tasks and update the database.
        if tasks:
            tasks[0].start()
            num_launched += 1

            self.flow.pickle_dump()

        return num_launched

    def rapidfire(self, max_nlaunch=-1, max_loops=1, sleep_time=5):
        """
        Keeps submitting `Tasks` until we are out of jobs or no job is ready to run.

        Args:
            max_nlaunch: Maximum number of launches. default: no limit.
            max_loops: Maximum number of loops
            sleep_time: seconds to sleep between rapidfire loop iterations

        Returns:
            The number of tasks launched.
        """
        num_launched, do_exit, launched = 0, False, []

        for count in range(max_loops):
            if do_exit: break
            if count > 0: time.sleep(sleep_time)

            tasks = self.fetch_tasks_to_run()

            # I don't know why but we receive duplicated tasks.
            if any(task in launched for task in tasks):
                logger.critical("numtasks %d already in launched list:\n%s" % (len(tasks), launched))

            # Preventive test.
            tasks = [t for t in tasks if t not in launched]

            if not tasks: continue

            for task in tasks:
                fired = task.start()
                if fired:
                    launched.append(task)
                    num_launched += 1

                if num_launched >= max_nlaunch > 0:
                    logger.info('num_launched >= max_nlaunch, breaking submission loop')
                    do_exit = True
                    break

        # Update the database.
        self.flow.pickle_dump()

        return num_launched

    def fetch_tasks_to_run(self):
        """
        Return the list of tasks that can be submitted.
        Empty list if no task has been found.
        """
        tasks_to_run = []

        for work in self.flow:
            tasks_to_run.extend(work.fetch_alltasks_to_run())

        return tasks_to_run


class PyFlowSchedulerError(Exception):
    """Exceptions raised by `PyFlowScheduler`."""


class BaseScheduler(metaclass=abc.ABCMeta):
    """
    This object schedules the submission of the tasks in a |Flow|.
    There are two types of errors that might occur during the execution of the jobs:

        #. Python exceptions
        #. Errors in the ab-initio code

    Python exceptions are easy to detect and are usually due to a bug in the python code or random errors such as IOError.
    The set of errors in the ab-initio is much much broader. It includes wrong input data, segmentation
    faults, problems with the resource manager, etc. The flow tries to handle the most common cases
    but there's still a lot of room for improvement.
    Note, in particular, that `PyFlowScheduler` will shutdown automatically in the following cases:

        #. The number of python exceptions is > max_num_pyexcs

        #. The number of task errors (i.e. the number of tasks whose status is S_ERROR) is > max_num_abierrs

        #. The number of jobs launched becomes greater than (`safety_ratio` * total_number_of_tasks).

        #. The scheduler will send an email to the user (specified by `mailto`) every `remindme_s` seconds.
           If the mail cannot be sent, the scheduler will automatically shutdown.
           This check prevents the scheduler from being trapped in an infinite loop.
    """
    # Configuration file.
    YAML_FILE = "scheduler.yml"

    USER_CONFIG_DIR = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")

    Error = PyFlowSchedulerError

    @classmethod
    def autodoc(cls) -> str:
        """Return String with scheduler options."""
        i = cls.__init__.__doc__.index("Args:")
        return cls.__init__.__doc__[i+5:]

    def __init__(self, **kwargs):
        """
        Args:
            weeks: number of weeks to wait (DEFAULT: 0).
            days: number of days to wait (DEFAULT: 0).
            hours: number of hours to wait (DEFAULT: 0).
            minutes: number of minutes to wait (DEFAULT: 0).
            seconds: number of seconds to wait (DEFAULT: 0).
            mailto: The scheduler will send an email to `mailto` every `remindme_s` seconds.
                (DEFAULT: None i.e. not used).
            verbose: (int) verbosity level. (DEFAULT: 0)
            use_dynamic_manager: "yes" if the |TaskManager| must be re-initialized from
                file before launching the jobs. (DEFAULT: "no")
            max_njobs_inqueue: Limit on the number of jobs that can be present in the queue. (DEFAULT: 200)
            max_ncores_used: Maximum number of cores that can be used by the scheduler.
            remindme_s: The scheduler will send an email to the user specified
                by `mailto` every `remindme_s` seconds. (int, DEFAULT: 1 day).
            max_num_pyexcs: The scheduler will exit if the number of python exceptions is > max_num_pyexcs
                (int, DEFAULT: 0)
            max_num_abierrs: The scheduler will exit if the number of errored tasks is > max_num_abierrs
                (int, DEFAULT: 0)
            safety_ratio: The scheduler will exits if the number of jobs launched becomes greater than
               `safety_ratio` * total_number_of_tasks_in_flow. (int, DEFAULT: 5)
            max_nlaunches: Maximum number of tasks launched in a single iteration of the scheduler.
                (DEFAULT: -1 i.e. no limit)
            debug: Debug level. Use 0 for production (int, DEFAULT: 0)
            fix_qcritical: "yes" if the launcher should try to fix QCritical Errors (DEFAULT: "no")
            rmflow: If "yes", the scheduler will remove the flow directory if the calculation
                completed successfully. (DEFAULT: "no")
            killjobs_if_errors: "yes" if the scheduler should try to kill all the running jobs
                before exiting due to an error. (DEFAULT: "yes")
        """
        self.init_kwargs = kwargs.copy()

        # Options passed to the apscheduler scheduler.
        self.sched_options = AttrDict(
            weeks=kwargs.pop("weeks", 0),
            days=kwargs.pop("days", 0),
            hours=kwargs.pop("hours", 0),
            minutes=kwargs.pop("minutes", 0),
            seconds=kwargs.pop("seconds", 0),
        )
        if all(not v for v in self.sched_options.values()):
            raise self.Error("Please specify at least one option among: seconds, minutes, hours ...")

        self.mailto = kwargs.pop("mailto", None)
        self.verbose = int(kwargs.pop("verbose", 0))
        self.use_dynamic_manager = as_bool(kwargs.pop("use_dynamic_manager", False))
        self.max_njobs_inqueue = kwargs.pop("max_njobs_inqueue", 200)
        self.max_ncores_used = kwargs.pop("max_ncores_used", None)

        self.remindme_s = float(kwargs.pop("remindme_s", 1 * 24 * 3600))
        self.max_num_pyexcs = int(kwargs.pop("max_num_pyexcs", 0))
        self.max_num_abierrs = int(kwargs.pop("max_num_abierrs", 0))
        self.safety_ratio = int(kwargs.pop("safety_ratio", 5))
        #self.max_etime_s = kwargs.pop("max_etime_s", 14 * 3600)
        self.max_nlaunches = kwargs.pop("max_nlaunches", -1)
        self.debug = kwargs.pop("debug", 0)
        self.fix_qcritical = as_bool(kwargs.pop("fix_qcritical", False))
        self.rmflow = as_bool(kwargs.pop("rmflow", False))
        self.killjobs_if_errors = as_bool(kwargs.pop("killjobs_if_errors", True))

        # TODO: Add abinit_options, anaddb_options, exec_options

        if kwargs:
            raise self.Error("Unknown arguments `%s`" % str(kwargs))

        # Register the callaback in the scheduler
        if has_sched_v3:
            logger.warning("Using scheduler v >= 3.0.0")
            from apscheduler.schedulers.blocking import BlockingScheduler
            self.sched = BlockingScheduler()
            self.sched.add_job(self.callback, "interval", **self.sched_options)
        else:
            from apscheduler.scheduler import Scheduler
            self.sched = Scheduler(standalone=True)
            self.sched.add_interval_job(self.callback, **self.sched_options)

        self.nlaunch = 0
        self.num_reminders = 1
        self.start_time = None

        # Used to keep track of the exceptions raised while the scheduler is running
        self.exceptions = deque(maxlen=self.max_num_pyexcs + 10)

        # Used to push additional info during the execution.
        self.history = deque(maxlen=200)

    @classmethod
    def from_file(cls, filepath: str) -> BaseScheduler:
        """Read the configuration parameters from a Yaml file."""
        with open(filepath, "rt") as fh:
            return cls(**yaml_safe_load(fh))

    @classmethod
    def from_string(cls, s: str) -> BaseScheduler:
        """Create an istance from string s containing a YAML dictionary."""
        stream = StringIO(s)
        stream.seek(0)
        return cls(**yaml_safe_load(stream))

    @classmethod
    def from_user_config(cls) -> BaseScheduler:
        """
        Initialize the :class:`PyFlowScheduler` from the YAML file 'scheduler.yml'.
        Search first in the working directory and then in the configuration directory of abipy.

        Raises:
            `RuntimeError` if file is not found.
        """
        # Try in the current directory.
        path = os.path.join(os.getcwd(), cls.YAML_FILE)
        if os.path.exists(path):
            return cls.from_file(path)

        # Try in the configuration directory.
        path = os.path.join(cls.USER_CONFIG_DIR, cls.YAML_FILE)
        if os.path.exists(path):
            return cls.from_file(path)

        raise cls.Error("Cannot locate %s neither in current directory nor in %s" % (cls.YAML_FILE, path))

    def __str__(self):
        """String representation."""
        lines = [self.__class__.__name__ + ", Pid: %d" % self.pid]
        app = lines.append
        app("Scheduler options:\n%s" % str(self.sched_options))

        return "\n".join(lines)

    @abc.abstractmethod
    def add_flow(self, flow, **kwargs) -> None:
        """
        Add a flow to the scheduler.
        """

    @abc.abstractmethod
    def start(self) -> int:
        """
        Starts the scheduler. Returns 0 if success.
        """

    @abc.abstractmethod
    def callback(self):
        """The function that will be executed by the scheduler."""

    @lazy_property
    def pid(self) -> int:
        """The pid of the process associated to the scheduler."""
        return os.getpid()

    @property
    def num_excs(self) -> int:
        """Number of exceptions raised so far."""
        return len(self.exceptions)

    def get_delta_etime(self):
        """Returns a `timedelta` object representing with the elapsed time."""
        return datetime.timedelta(seconds=(time.time() - self.start_time))

    def cancel_jobs_if_requested(self, flow):

        if not self.killjobs_if_errors: return
        cprint("killjobs_if_errors set to 'yes'. Killing jobs before aborting the flow.", "yellow")
        try:
            num_cancelled = 0
            for task in flow.iflat_tasks():
                num_cancelled += task.cancel()
            cprint("Killed %d tasks" % num_cancelled, "yellow")
        except Exception as exc:
            cprint("Exception while trying to kill jobs:\n%s" % str(exc), "red")

    def _accept_flow(self, flow) -> None:

        # Check if we are already using a scheduler to run this flow
        flow.check_pid_file()
        flow.set_spectator_mode(False)

        # Build dirs and files (if not yet done)
        flow.build()

        errors = flow.look_before_you_leap()
        if errors:
            raise self.Error(str(errors))

    def restart_unconverged(self, flow, max_nlaunch, excs):
        if max_nlaunch <= 0: return 0

        for task in flow.unconverged_tasks:
            try:
                logger.info("Trying to restart task: `%s`" % repr(task))
                fired = task.restart()
                if fired:
                    self.nlaunch += 1
                    max_nlaunch -= 1
                    if max_nlaunch <= 0:
                        logger.info("Restart: too many jobs in the queue, returning")
                        flow.pickle_dump()
                        return 0

            except task.RestartError:
                excs.append(straceback())

        return max_nlaunch

    def try_to_fix_flow(self, flow):

        # Temporarily disabled by MG because I don't know if fix_critical works after the
        # introduction of the new qadapters
        # reenabled by MsS disable things that do not work at low level
        # fix only prepares for restarting, and sets to ready
        if self.fix_qcritical:
            nfixed = flow.fix_queue_critical()
            if nfixed: print("Fixed %d QCritical error(s)" % nfixed)

        nfixed = flow.fix_abicritical()
        if nfixed: print("Fixed %d AbiCritical error(s)" % nfixed)

    def check_deadlocks(self, flow) -> List[str]:
        err_lines = []

        g = flow.find_deadlocks()
        if g.deadlocked:
            # Check the flow again so that status are updated.
            flow.check_status()

            g = flow.find_deadlocks()
            print("deadlocked:", len(g.deadlocked), ", runnables:", len(g.runnables), ", running:", len(g.running))
            if g.deadlocked and not g.runnables and not g.running:
                err_lines.append("No runnable job with deadlocked tasks:\n%s." % str(g.deadlocked))

        if not g.runnables and not g.running:
            # Check the flow again so that status are updated.
            flow.check_status()
            g = flow.find_deadlocks()
            if not g.runnables and not g.running:
                err_lines.append("No task is running and cannot find other tasks to submit.")

        return err_lines

    def send_email(self, msg, tag=None) -> int:
        """
        Send an e-mail before completing the shutdown.
        Returns 0 if success. Relies on _send_email method provided by subclass.
        """
        try:
            return self._send_email(msg, tag)
        except Exception:
            self.exceptions.append(straceback())
            return -2


class PyFlowScheduler(BaseScheduler):

    @property
    def pid_file(self) -> str:
        """
        Absolute path of the file with the pid.
        The file is located in the workdir of the flow
        """
        return self._pid_file

    @property
    def flow(self):
        """`Flow`."""
        try:
            return self._flow
        except AttributeError:
            return None

    def add_flow(self, flow):
        """
        Add a flow to the scheduler.
        """
        if hasattr(self, "_flow"):
            raise self.Error(f"Only one Flow can be added to a {self.__class__.__name__} scheduler.")

        self._accept_flow(flow)

        with open(flow.pid_file, "wt") as fh:
            fh.write(str(self.pid))

        self._pid_file = flow.pid_file
        self._flow = flow

    def start(self) -> int:
        """
        Starts the scheduler. Returns 0 if success.
        This method blocks until the flow is completed or some exception is raised.
        """
        self.history.append("Started on %s" % time.asctime())
        self.start_time = time.time()

        flow = self.flow
        if flow is not None:
            # Try to run the job immediately.
            # If something goes wrong return without initializing the scheduler.
            self._runem_all()

            if self.exceptions:
                self.cleanup()
                self.send_email(msg="Error while trying to run the flow for the first time!\n %s" % self.exceptions)
                return 1

        # Start the scheduler loop.
        try:
            self.sched.start()
            return 0

        except KeyboardInterrupt:
            self.shutdown(msg="KeyboardInterrupt from user")
            if ask_yesno("Do you want to cancel all the jobs in the queue? [Y/n]"):
                print("Number of jobs cancelled:", flow.cancel())

            flow.pickle_dump()
            return -1

    def _runem_all(self):
        """
        This function checks the status of all tasks,
        tries to fix tasks that went unconverged, abicritical, or queuecritical
        and tries to run all the tasks that can be submitted.
        """
        excs = []
        flow = self.flow
        if flow is None: return

        if self.use_dynamic_manager:
            # Allow to change the manager at run-time
            from .tasks import TaskManager
            flow.change_manager(TaskManager.from_user_config())

        # Here we just count the number of tasks in the flow that are in RUNNING or SUBMITTED status.
        # This logic clearly breaks down if there are multiple schedulers running on the same machine
        # but it's easy to implement without having to contact the resource manager
        # by calling `flow.get_njobs_in_queue()` that is quite expensive (and the sysadmin won't be happy!)
        nqjobs = (len(list(flow.iflat_tasks(status=flow.S_RUN))) +
                  len(list(flow.iflat_tasks(status=flow.S_SUB))))

        if nqjobs >= self.max_njobs_inqueue:
            print(f"Too many jobs in the queue: {nqjobs} >= {self.max_njobs_inqueue}.\n",
                  "No job will be submitted.")
            flow.check_status(show=False)
            return

        max_nlaunch = self.max_njobs_inqueue - nqjobs if self.max_nlaunches == -1 else \
                      min(self.max_njobs_inqueue - nqjobs, self.max_nlaunches)

        # check status.
        flow.check_status(show=False)

        # This check is not perfect, we should make a list of tasks to submit
        # and then select a subset so that we don't exceeed max_ncores_used
        # Many sections of this code should be rewritten though.
        #if self.max_ncores_used is not None and flow.ncores_used > self.max_ncores_used:
        if self.max_ncores_used is not None and flow.ncores_allocated > self.max_ncores_used:
            print("Cannot exceed max_ncores_used %s" % self.max_ncores_used,
                  ", ncores_allocated:", flow.ncores_allocated)
            return

        # Try to restart unconverged tasks.
        max_nlaunch = self.restart_unconverged(flow, max_nlaunch, excs)
        if max_nlaunch <= 0: return

        self.try_to_fix_flow(flow)

        # Update the pickle file.
        flow.pickle_dump()

        # Submit the tasks that are ready.
        try:
            nlaunch = PyLauncher(flow).rapidfire(max_nlaunch=max_nlaunch, sleep_time=10)
            self.nlaunch += nlaunch
            if nlaunch:
                cprint("[%s] Number of launches: %d" % (time.asctime(), nlaunch), "yellow")

        except Exception:
            excs.append(straceback())

        flow.show_status()

        if excs:
            logger.critical("*** Scheduler exceptions:\n *** %s" % "\n".join(excs))
            self.exceptions.extend(excs)

    def callback(self):
        """The function that will be executed by the scheduler."""
        try:
            return self._callback()
        except Exception:
            # All exceptions raised here will trigger the shutdown!
            s = straceback()
            self.exceptions.append(s)
            # This is useful when debugging
            #self.cancel_jobs_if_requested(self.flow)
            self.shutdown(msg="Exception raised in callback!\n" + s)

    def _callback(self):
        """The actual callback."""
        if self.debug: print(">>>>> _callback: Number of open file descriptors: %s" % get_open_fds())

        self._runem_all()

        flow = self.flow
        all_ok = flow.all_ok

        if all_ok:
            # Mission accomplished. Shutdown the scheduler.
            print("Calling flow.finalize() ...")
            flow.finalize()
            if self.rmflow:
                print("Flow directory will be removed...")
                try:
                    flow.rmtree()
                except Exception:
                    logger.warning("Ignoring exception while trying to remove flow dir.")

            return self.shutdown(msg="All tasks have reached S_OK. Will shutdown the scheduler and exit")

        # Handle failures.
        err_lines = []

        # Shall we send a reminder to the user?
        delta_etime = self.get_delta_etime()

        if delta_etime.total_seconds() > self.num_reminders * self.remindme_s:
            self.num_reminders += 1
            msg = ("Just to remind you that the scheduler with pid %s, flow %s\n has been running for %s " %
                  (self.pid, flow, delta_etime))
            retcode = self.send_email(msg, tag="[REMINDER]")

            if retcode:
                msg += ("\nThe scheduler tried to send an e-mail to remind the user\n" +
                        " but send_email returned %d. Error is not critical though!" % retcode)
                print(msg)

        #if delta_etime.total_seconds() > self.max_etime_s:
        #    err_lines.append("\nExceeded max_etime_s %s. Will shutdown the scheduler and exit" % self.max_etime_s)

        # Too many exceptions. Shutdown the scheduler.
        if self.num_excs > self.max_num_pyexcs:
            msg = "Number of exceptions %s > %s. Will shutdown the scheduler and exit" % (
                self.num_excs, self.max_num_pyexcs)
            err_lines.append(boxed(msg))

        # Paranoid check: disable the scheduler if we have submitted
        # too many jobs (it might be due to some bug or other external reasons
        # such as race conditions between difference callbacks!)
        if self.nlaunch > self.safety_ratio * flow.num_tasks:
            msg = "Too many jobs launched %d. Total number of tasks = %s, Will shutdown the scheduler and exit" % (
                self.nlaunch, flow.num_tasks)
            err_lines.append(boxed(msg))

        # Count the number of tasks with status == S_ERROR.
        if flow.num_errored_tasks > self.max_num_abierrs:
            msg = "Number of tasks with ERROR status %s > %s. Will shutdown the scheduler and exit" % (
                flow.num_errored_tasks, self.max_num_abierrs)
            err_lines.append(boxed(msg))

        # Test for deadlocks.
        errors = self.check_deadlocks(flow)
        if errors:
            err_lines.extend(errors)

        # Something wrong. Quit
        if err_lines:
            self.cancel_jobs_if_requested(flow)
            self.shutdown("\n".join(err_lines))

        return len(self.exceptions)

    def cleanup(self) -> None:
        """Cleanup routine: remove the pid file and save the pickle database"""
        try:
            os.remove(self.pid_file)
        except OSError as exc:
            logger.critical("Could not remove pid_file: %s", exc)

        # Save the final status of the flow.
        self.flow.pickle_dump()

    def shutdown(self, msg: str) -> None:
        """Shutdown the scheduler."""
        flow = self.flow
        all_ok = flow.all_ok

        try:
            self.cleanup()
            self.history.append("Completed on: %s" % time.asctime())
            self.history.append("Elapsed time: %s" % self.get_delta_etime())

            if self.debug:
                print(">>>>> shutdown: Number of open file descriptors: %s" % get_open_fds())

            retcode = self.send_email(msg)
            if self.debug: print("send_mail retcode", retcode)

            # Write file with the list of exceptions:
            if self.exceptions:
                dump_file = os.path.join(flow.workdir, "_exceptions")
                with open(dump_file, "wt") as fh:
                    fh.writelines(self.exceptions)
                    fh.write("Shutdown message:\n%s" % msg)

            lines = []
            app = lines.append
            app("Submitted on: %s" % time.ctime(self.start_time))
            app("Completed on: %s" % time.asctime())
            app("Elapsed time: %s" % str(self.get_delta_etime()))

            if all_ok:
                app("Flow completed successfully")
            else:
                app("Flow %s didn't complete successfully" % repr(flow.workdir))
                app("use `abirun.py FLOWDIR debug` to analyze the problem.")
                app("Shutdown message:\n%s" % msg)

            print("")
            print("\n".join(lines))
            print("")

        finally:
            # Shutdown the scheduler thus allowing the process to exit.
            logger.debug('This should be the shutdown of the scheduler')

            # Unschedule all the jobs before calling shutdown
            #self.sched.print_jobs()
            if not has_sched_v3:
                #self.sched.print_jobs()
                for job in self.sched.get_jobs():
                    self.sched.unschedule_job(job)
                self.sched.shutdown()
            else:
                self.sched.shutdown(wait=False)

            # Uncomment the line below if shutdown does not work!
            #os.system("kill -9 %d" % os.getpid())

    def _send_email(self, msg, tag):
        if self.mailto is None: return -1

        header = msg.splitlines()
        app = header.append
        app("Submitted on: %s" % time.ctime(self.start_time))
        app("Completed on: %s" % time.asctime())
        app("Elapsed time: %s" % str(self.get_delta_etime()))

        # Add the status of the flow.
        flow = self.flow
        app("Number of errored tasks: %d" % flow.num_errored_tasks)
        app("Number of unconverged tasks: %d" % flow.num_unconverged_tasks)
        strio = StringIO()
        strio.writelines("\n".join(header) + 4 * "\n")
        flow.show_status(stream=strio)

        if self.exceptions:
            # Report the list of exceptions.
            strio.writelines(self.exceptions)

        if tag is None:
            tag = " [ALL OK]" if flow.all_ok else " [WARNING]"

        return sendmail(subject=flow.name + tag, text=strio.getvalue(), mailto=self.mailto)


class MultiFlowScheduler(BaseScheduler):

    # TODO: history, logging, shutdown better treatment of exceptions....

    def __init__(self, sqldb_path, **kwargs):
        super().__init__(**kwargs)
        self.flows = []

        self.completed_flows = []
        self.errored_flows = []
        self._errored_flow_ids = []

        self.incoming_flow_queue = Queue()

        import threading
        self._lock = threading.Lock()

        self.sqldb_path = sqldb_path
        self.create_sqldb()

    def add_flow(self, flow, user_message, priority=None):
        """
        Add a flow to the scheduler.
        """
        # TODO: Should check the pid file
        self._accept_flow(flow)
        flow.set_user_message(user_message)
        self.incoming_flow_queue.put(flow)

    def register_flow_exception(self, flow_idx, exc) -> None:
        flow = self.flows[flow_idx]
        self.history.append(f"Exception for {repr(flow)}")
        self.history.append(straceback())
        self._errored_flow_ids.append(flow_idx)

    def handle_flow_exception(self):
        if not self._errored_flow_ids: return
        for idx in self._errored_flow_ids:
            flow = self.flows[idx]
            self.errored_flows.append(flow)

        self.flows = [flow for (idx, flow) in self.flows if idx not in set(self._errored_flow_ids)]
        self._errored_flow_ids = []

    def start(self) -> int:
        """
        Starts the scheduler. Returns 0 if success.
        This method blocks until the flow is completed or some exception is raised.
        """
        self.history.append("Started on %s" % time.asctime())
        self.start_time = time.time()
        self.sched.start()

    # TODO
    #def stop(self):
    #def restart(self):

    def sql_connect(self):
        import sqlite3
        con = sqlite3.connect(self.sqldb_path, check_same_thread=True,
                              detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES)
        con.row_factory = sqlite3.Row
        return con

    def create_sqldb(self) -> None:
        if os.path.exists(self.sqldb_path):
            return

        with self._lock, self.sql_connect() as con:
            # Create table
            cur = con.cursor()
            cur.execute("""CREATE TABLE flows (
                        status TEXT NOT NULL,
                        formula TEXT NOT NULL,
                        workdir TEXT NOT NULL,
                        pyfile TEXT NOT NULL,
                        flow_id PRIMARY KEY,
                        upload_date TIMESTAMP,
                        end_date TIMESTAMP DEFAULT NULL,
                        user_message TEXT NOT NULL
                        );
                        """)
        con.close()

    def get_incoming_flows(self) -> List:
        flows = []
        while True:
            try:
                in_flow = self.incoming_flow_queue.get_nowait()
                flows.append(in_flow)
            except Empty:
                break

        if flows:
            def get_record(flow):
                formula = flow[0][0].input.structure.formula
                return (str(flow.status), formula, flow.workdir, os.path.basename(flow.pyfile),
                        flow.node_id, datetime.datetime.now(), None, flow.user_message)

            with self.sql_connect() as con:
                cur = con.cursor()
                values = [get_record(flow) for flow in flows]
                cur.executemany("INSERT INTO flows VALUES (?, ?, ?, ?, ?, ?, ?, ?)", values)
            con.close()

        return flows

    def get_dataframe(self):
        with self.sql_connect() as con:
            df = pd.read_sql_query("SELECT * FROM flows", con)
            #print("dtype", df["upload_date"].dtype)
            return df

    def get_json_status(self):
        # https://stackoverflow.com/questions/25455067/pandas-dataframe-datetime-index-doesnt-survive-json-conversion-and-reconversion
        status = dict(
            dataframe=self.get_dataframe().to_json() #, date_format='iso'#date_unit='ns'),
        )
        return status

    def get_flow_and_status_by_nodeid(self, node_id):
        from abipy.flowtk.flows import Flow
        with self.sql_connect() as con:
            cur = con.cursor()
            cur.execute("SELECT workdir, status FROM flows WHERE flow_id = ?", [node_id])
            row = cur.fetchone()
        con.close()

        if row and os.path.exists(row["workdir"]):
            return Flow.from_file(row["workdir"]), row["status"]
        else:
            return None, None

    def get_sql_rows_with_node_ids(self, node_id_list):
        with self.sql_connect() as con:
            cur = con.cursor()
            query = "SELECT * FROM flows WHERE flow_id IN (%s)" % ','.join('?' * len(node_id_list))
            cur.execute(query, node_id_list)
            rows = cur.fetchall()
        con.close()
        return rows

    def groupby_status(self):
        from collections import defaultdict
        d = defaultdict(list)

        with self.sql_connect() as con:
            cur = con.cursor()
            cur.execute("SELECT * FROM flows")
            rows = cur.fetchall()

        con.close()

        for row in rows:
            d[row["status"]].append(row)

        return d

    def remove_flows_with_status(self, status):
        if status == "Running":
           raise ValueError("You cannot remove a flow that is in `Running` mode!")

        count = 0
        with self._lock, self.sql_connect() as con:
            cur = con.cursor()
            cur.execute("DELETE FROM flows WHERE status = ?", [status])
            rows = cur.fetchall()
            if rows is None:
                con.close()
                return 0

            for row in rows:
                workdir = row["workdir"]
                if not os.path.exists(workdir): continue
                os.rmdir(workdir)
                count += 1

        con.close()
        return count

    def callback(self):
        """The function that will be executed by the scheduler."""
        #locked = self._lock.acquire(blocking=False, timeout=-1)
        #if not locked: return
        #with self._lock:

        new_flows = self.get_incoming_flows()
        if new_flows:
            self.flows.extend(new_flows)

        if not self.flows: return
        excs = []

        # check status.
        for flow in self.flows:
            flow.check_status(show=False)

        # Here we just count the number of tasks in the flow that are in RUNNING or SUBMITTED status.
        # This logic clearly breaks down if there are multiple schedulers running on the same machine
        # but it's easy to implement without having to contact the resource manager
        # by calling `flow.get_njobs_in_queue()` that is quite expensive (and the sysadmin won't be happy!)
        nqjobs = 0
        for flow in self.flows:
            nqjobs += (len(list(flow.iflat_tasks(status=flow.S_RUN))) +
                       len(list(flow.iflat_tasks(status=flow.S_SUB))))

        if nqjobs >= self.max_njobs_inqueue:
            print(f"Too many jobs in the queue: {nqjobs} >= {self.max_njobs_inqueue}.\n",
                  "No job will be submitted.")
            for flow in self.flows:
                flow.check_status(show=False)
            return

        max_nlaunch = self.max_njobs_inqueue - nqjobs if self.max_nlaunches == -1 else \
                      min(self.max_njobs_inqueue - nqjobs, self.max_nlaunches)

        # This check is not perfect, we should make a list of tasks to submit
        # and then select a subset so that we don't exceeed max_ncores_used
        # Many sections of this code should be rewritten though.
        #if self.max_ncores_used is not None and flow.ncores_used > self.max_ncores_used:
        ncores_allocated = sum(flow.ncores_allocated for flow in self.flows)
        if self.max_ncores_used is not None and ncores_allocated > self.max_ncores_used:
            print("Cannot exceed max_ncores_used %s" % self.max_ncores_used,
                  ", ncores_allocated:", ncores_allocated)
            return

        # Try to restart unconverged tasks.
        for flow in self.flows:
            max_nlaunch = self.restart_unconverged(flow, max_nlaunch, excs)
            if max_nlaunch <= 0: return

        for flow in self.flows:
            self.try_to_fix_flow(flow)
            # Update the pickle file.
            flow.pickle_dump()

        #with self.handle_flow_exceptions:

        for i, flow in enumerate(self.flows):
            # Submit the tasks that are ready.
            try:
                if max_nlaunch > 0:
                    nlaunch = PyLauncher(flow).rapidfire(max_nlaunch=max_nlaunch, sleep_time=10)
                    self.nlaunch += nlaunch
                    max_nlaunch -= nlaunch
                    if nlaunch:
                        cprint("[%s] Number of launches: %d" % (time.asctime(), nlaunch), "yellow")

            except Exception as exc:
                self.register_flow_exception(i, exc)

        self.handle_flow_exception()

        for flow in self.flows:
            flow.show_status()

        #if max_nlaunch <= 0: return

        self.update_flows_and_slqdb()

    def update_flows_and_slqdb(self):

        done = []
        for i, flow in enumerate(self.flows):
            if flow.all_ok:
                flow.finalize()
                done.append(i)

        if done:
            for i in done:
                flow = self.flows[i]
                self.completed_flows.append(flow)
            self.flows = [self.flows[i] for i in range(len(self.flows)) if i not in set(done)]

        locked = []
        for i, flow in enumerate(self.flows):
            errs = self.check_deadlocks(flow)
            if errs:
                msg = "\n".join(errs)
                locked.append(i)
                flow.set_status(flow.S_ERROR, msg)

        if locked:
            for i in locked:
                flow = self.flows[i]
                self.errored_flows.append(flow)

            self.flows = [self.flows[i] for i in range(len(self.flows)) if i not in set(locked)]

        with self._lock, self.sql_connect() as con:
            now = datetime.datetime.now()
            cur = con.cursor()
            query = "UPDATE flows SET status = ?, end_date = ? WHERE flow_id = ?"

            values = [(str(flow.S_RUN), now, flow.node_id) for flow in self.flows]

            if self.completed_flows:
               values.extend([(str(flow.status), now, flow.node_id) for flow in self.completed_flows])
               self.completed_flows = []

            if self.errored_flows:
               values.extend([(str(flow.S_ERROR), now, flow.node_id) for flow in self.errored_flows])
               self.errored_flows = []

            cur.executemany(query, values)

        con.close()


def print_flowsdb_file(filepath: str):
    """
    Print flows.db file to terminal.
    """
    import sqlite3
    from abipy.tools.printing import print_dataframe
    with sqlite3.connect(filepath) as con:
        df = pd.read_sql_query("SELECT * FROM flows", con)
        #print(type(df["upload_date"]))
        print_dataframe(df, title=filepath)


def sendmail(subject: str, text: str, mailto: str,
             sender: Optional[str] = None) -> int:
    """
    Sends an e-mail with unix sendmail.

    Args:
        subject: String with the subject of the mail.
        text: String with the body of the mail.
        mailto: String or list of string with the recipients.
        sender: string with the sender address. If sender is None, username@hostname is used.

    Returns: Exit status
    """
    def user_at_host():
        from socket import gethostname
        return os.getlogin() + "@" + gethostname()

    # Body of the message.
    try:
        sender = user_at_host() if sender is None else sender
    except OSError:
        sender = 'abipyscheduler@youknowwhere'

    if is_string(mailto): mailto = [mailto]

    from email.mime.text import MIMEText
    mail = MIMEText(text)
    mail["Subject"] = subject
    mail["From"] = sender
    mail["To"] = ", ".join(mailto)

    msg = mail.as_string()

    # sendmail works much better than the python interface.
    # Note that sendmail is available only on Unix-like OS.
    from subprocess import Popen, PIPE

    _sendmail = which("sendmail")
    if _sendmail is None: return -1
    # msg is string not bytes so must use universal_newlines
    p = Popen([_sendmail, "-t"], stdin=PIPE, stderr=PIPE, universal_newlines=True)

    outdata, errdata = p.communicate(msg)
    return len(errdata)


#def __test_sendmail():
#    retcode = sendmail("sendmail_test", text="hello\nworld", mailto="nobody@nowhere.com")
#    print("Retcode", retcode)
#    assert retcode == 0
