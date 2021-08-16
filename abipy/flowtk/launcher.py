# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""Tools for the submission of Tasks."""

import os
import time
import ruamel.yaml as yaml
import pickle
import apscheduler
has_sched_v3 = apscheduler.version >= "3.0.0"

from collections import deque
from datetime import timedelta
from io import StringIO
from monty.io import get_open_fds
from monty.string import boxed, is_string
from monty.os.path import which
from monty.collections import AttrDict, dict2namedtuple
from monty.termcolor import cprint
from pymatgen.util.io_utils import ask_yesno
from .utils import as_bool, File, Directory
from . import qutils as qu

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "ScriptEditor",
    "PyLauncher",
    "PyFlowScheduler",
]


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def yaml_safe_load(string):
    return yaml.YAML(typ='safe', pure=True).load(string)


class ScriptEditor(object):
    """Simple editor that simplifies the writing of shell scripts"""
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
    """This object handle the submission of the tasks contained in a |Flow|."""
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
            if do_exit:
                break
            if count > 0:
                time.sleep(sleep_time)

            tasks = self.fetch_tasks_to_run()

            # I don't know why but we receive duplicated tasks.
            if any(task in launched for task in tasks):
                logger.critical("numtasks %d already in launched list:\n%s" % (len(tasks), launched))

            # Preventive test.
            tasks = [t for t in tasks if t not in launched]

            if not tasks:
                continue

            for task in tasks:
                fired = task.start()
                if fired:
                    launched.append(task)
                    num_launched += 1

                if num_launched >= max_nlaunch > 0:
                    logger.info('num_launched >= max_nlaunch, going back to sleep')
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


class PyFlowScheduler(object):
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
           If the mail cannot be sent, the scheduler will shutdown automatically.
           This check prevents the scheduler from being trapped in an infinite loop.
    """
    # Configuration file.
    YAML_FILE = "scheduler.yml"
    USER_CONFIG_DIR = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")

    Error = PyFlowSchedulerError

    @classmethod
    def autodoc(cls):
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
            killjobs_if_errors: "yes" if the scheduler should try to kill all the runnnig jobs
                before exiting due to an error. (DEFAULT: "yes")
        """
        # Options passed to the apscheduler scheduler.
        self.sched_options = AttrDict(
            weeks=kwargs.pop("weeks", 0),
            days=kwargs.pop("days", 0),
            hours=kwargs.pop("hours", 0),
            minutes=kwargs.pop("minutes", 0),
            seconds=kwargs.pop("seconds", 0),
        )
        if all(not v for v in self.sched_options.values()):
            raise self.Error("Wrong set of options passed to the scheduler.")

        self.mailto = kwargs.pop("mailto", None)
        self.verbose = int(kwargs.pop("verbose", 0))
        self.use_dynamic_manager = as_bool(kwargs.pop("use_dynamic_manager", False))
        self.max_njobs_inqueue = kwargs.pop("max_njobs_inqueue", 200)
        self.max_ncores_used = kwargs.pop("max_ncores_used", None)
        self.contact_resource_manager = as_bool(kwargs.pop("contact_resource_manager", False))

        self.remindme_s = float(kwargs.pop("remindme_s", 1 * 24 * 3600))
        self.max_num_pyexcs = int(kwargs.pop("max_num_pyexcs", 0))
        self.max_num_abierrs = int(kwargs.pop("max_num_abierrs", 0))
        self.safety_ratio = int(kwargs.pop("safety_ratio", 5))
        #self.max_etime_s = kwargs.pop("max_etime_s", )
        self.max_nlaunches = kwargs.pop("max_nlaunches", -1)
        self.debug = kwargs.pop("debug", 0)
        self.fix_qcritical = as_bool(kwargs.pop("fix_qcritical", False))
        self.rmflow = as_bool(kwargs.pop("rmflow", False))
        self.killjobs_if_errors = as_bool(kwargs.pop("killjobs_if_errors", True))

        if kwargs:
            raise self.Error("Unknown arguments `%s`" % str(kwargs))

        if has_sched_v3:
            logger.warning("Using scheduler v >= 3.0.0")
            from apscheduler.schedulers.blocking import BlockingScheduler
            self.sched = BlockingScheduler()
        else:
            from apscheduler.scheduler import Scheduler
            self.sched = Scheduler(standalone=True)

        self.nlaunch = 0
        self.num_reminders = 1

        # Used to keep track of the exceptions raised while the scheduler is running
        self.exceptions = deque(maxlen=self.max_num_pyexcs + 10)

        # Used to push additional info during the execution.
        self.history = deque(maxlen=100)

    @classmethod
    def from_file(cls, filepath):
        """Read the configuration parameters from a Yaml file."""
        with open(filepath, "rt") as fh:
            return cls(**yaml_safe_load(fh))

    @classmethod
    def from_string(cls, s):
        """Create an istance from string s containing a YAML dictionary."""
        stream = StringIO(s)
        stream.seek(0)
        return cls(**yaml_safe_load(stream))

    @classmethod
    def from_user_config(cls):
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

    @property
    def pid(self):
        """The pid of the process associated to the scheduler."""
        try:
            return self._pid
        except AttributeError:
            self._pid = os.getpid()
            return self._pid

    @property
    def pid_file(self):
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

    @property
    def num_excs(self):
        """Number of exceptions raised so far."""
        return len(self.exceptions)

    def get_delta_etime(self):
        """Returns a `timedelta` object representing with the elapsed time."""
        return timedelta(seconds=(time.time() - self.start_time))

    def add_flow(self, flow):
        """
        Add a flow to the scheduler.
        """
        if hasattr(self, "_flow"):
            raise self.Error("Only one flow can be added to the scheduler.")

        # Check if we are already using a scheduler to run this flow
        flow.check_pid_file()
        flow.set_spectator_mode(False)

        # Build dirs and files (if not yet done)
        flow.build()

        with open(flow.pid_file, "wt") as fh:
            fh.write(str(self.pid))

        self._pid_file = flow.pid_file
        self._flow = flow

    def start(self):
        """
        Starts the scheduler in a new thread. Returns 0 if success.
        In standalone mode, this method will block until there are no more scheduled jobs.
        """
        self.history.append("Started on %s" % time.asctime())
        self.start_time = time.time()

        if has_sched_v3:
            self.sched.add_job(self.callback, "interval", **self.sched_options)
        else:
            self.sched.add_interval_job(self.callback, **self.sched_options)

        flow = self.flow

        errors = flow.look_before_you_leap()
        if errors:
            self.exceptions.append(errors)
            return 1

        # Try to run the job immediately.
        # If something goes wrong return without initializing the scheduler.
        self._runem_all()

        if self.exceptions:
            self.cleanup()
            self.send_email(msg="Error while trying to run the flow for the first time!\n %s" % self.exceptions)
            return 1

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

        # Allow to change the manager at run-time
        if self.use_dynamic_manager:
            from .tasks import TaskManager
            new_manager = TaskManager.from_user_config()
            for work in flow:
                work.set_manager(new_manager)

        nqjobs = 0
        if self.contact_resource_manager: # and flow.TaskManager.qadapter.QTYPE == "shell":
            # This call is expensive and therefore it's optional (must be activate in manager.yml)
            nqjobs = flow.get_njobs_in_queue()
            if nqjobs is None:
                nqjobs = 0
                if flow.manager.has_queue:
                    logger.warning('Cannot get njobs_inqueue')
        else:
            # Here we just count the number of tasks in the flow who are running.
            # This logic breaks down if there are multiple schedulers running
            # but it's easy to implement without having to contact the resource manager.
            nqjobs = (len(list(flow.iflat_tasks(status=flow.S_RUN))) +
                      len(list(flow.iflat_tasks(status=flow.S_SUB))))

        if nqjobs >= self.max_njobs_inqueue:
            print("Too many jobs in the queue: %s. No job will be submitted." % nqjobs)
            flow.check_status(show=False)
            return

        if self.max_nlaunches == -1:
            max_nlaunch = self.max_njobs_inqueue - nqjobs
        else:
            max_nlaunch = min(self.max_njobs_inqueue - nqjobs, self.max_nlaunches)

        # check status.
        flow.check_status(show=False)

        # This check is not perfect, we should make a list of tasks to sumbit
        # and select only the subset so that we don't exceeed mac_ncores_used
        # Many sections of this code should be rewritten.
        #if self.max_ncores_used is not None and flow.ncores_used > self.max_ncores_used:
        if self.max_ncores_used is not None and flow.ncores_allocated > self.max_ncores_used:
            print("Cannot exceed max_ncores_used %s" % self.max_ncores_used, ", ncores_allocated:", flow.ncores_allocated)
            return

        # Try to restart the unconverged tasks
        for task in flow.unconverged_tasks:
            try:
                logger.info("Trying to restart task: `%s`" % repr(task))
                fired = task.restart()
                if fired:
                    self.nlaunch += 1
                    max_nlaunch -= 1
                    if max_nlaunch == 0:
                        logger.info("Restart: too many jobs in the queue, returning")
                        flow.pickle_dump()
                        return

            except task.RestartError:
                excs.append(straceback())

        # Temporarily disabled by MG because I don't know if fix_critical works after the
        # introduction of the new qadapters
        # reenabled by MsS disable things that do not work at low level
        # fix only prepares for restarting, and sets to ready
        if self.fix_qcritical:
            nfixed = flow.fix_queue_critical()
            if nfixed: print("Fixed %d QCritical error(s)" % nfixed)

        nfixed = flow.fix_abicritical()
        if nfixed: print("Fixed %d AbiCritical error(s)" % nfixed)

        # update database
        flow.pickle_dump()

        # Submit the tasks that are ready.
        try:
            nlaunch = PyLauncher(flow).rapidfire(max_nlaunch=max_nlaunch, sleep_time=10)
            self.nlaunch += nlaunch
            if nlaunch:
                cprint("[%s] Number of launches: %d" % (time.asctime(), nlaunch), "yellow")

        except Exception:
            excs.append(straceback())

        # check status.
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
            #if self.killjobs_if_errors:
            #    print("Exception in callback, will cancel all tasks")
            #    try:
            #        for task in self.flow.iflat_tasks():
            #            task.cancel()
            #    except Exception:
            #        pass

            self.shutdown(msg="Exception raised in callback!\n" + s)

    def _callback(self):
        """The actual callback."""
        if self.debug:
            # Show the number of open file descriptors
            print(">>>>> _callback: Number of open file descriptors: %s" % get_open_fds())

        self._runem_all()

        flow = self.flow
        all_ok = flow.all_ok

        # Mission accomplished. Shutdown the scheduler.
        if all_ok:
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
                # Cannot send mail, shutdown now!
                msg += ("\nThe scheduler tried to send an e-mail to remind the user\n" +
                        " but send_email returned %d. Error is not critical though!" % retcode)
                print(msg)
                #err_lines.append(msg)

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

        # Something wrong. Quit
        if err_lines:
            # Cancel all jobs.
            if self.killjobs_if_errors:
                cprint("killjobs_if_errors set to 'yes'. Killing jobs before exiting.", "yellow")
                try:
                    num_cancelled = 0
                    for task in flow.iflat_tasks():
                        num_cancelled += task.cancel()
                    cprint("Killed %d tasks" % num_cancelled, "yellow")
                except Exception as exc:
                    cprint("Exception while trying to kill jobs:\n%s" % str(exc), "red")

            self.shutdown("\n".join(err_lines))

        return len(self.exceptions)

    def cleanup(self):
        """Cleanup routine: remove the pid file and save the pickle database"""
        try:
            os.remove(self.pid_file)
        except OSError as exc:
            logger.critical("Could not remove pid_file: %s", exc)

        # Save the final status of the flow.
        self.flow.pickle_dump()

    def shutdown(self, msg):
        """Shutdown the scheduler."""
        flow = self.flow

        try:
            self.cleanup()
            self.history.append("Completed on: %s" % time.asctime())
            self.history.append("Elapsed time: %s" % self.get_delta_etime())

            if self.debug:
                print(">>>>> shutdown: Number of open file descriptors: %s" % get_open_fds())

            retcode = self.send_email(msg)
            if self.debug:
                print("send_mail retcode", retcode)

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

            if flow.all_ok:
                app("Flow completed successfully")
            else:
                app("Flow %s didn't complete successfully" % repr(flow.workdir))
                app("use `abirun.py FLOWDIR debug` to analyze the problem.")
                app("Shutdown message:\n%s" % msg)

            print("")
            print("\n".join(lines))
            print("")

            if flow.all_ok:
                print("Calling flow.finalize() ...")
                flow.finalize()
                #print("finalized:", self.flow.finalized)
                if self.rmflow:
                    app("Flow directory will be removed...")
                    try:
                        flow.rmtree()
                    except Exception:
                        logger.warning("Ignoring exception while trying to remove flow dir.")

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

    def send_email(self, msg, tag=None):
        """
        Send an e-mail before completing the shutdown.
        Returns 0 if success.
        """
        try:
            return self._send_email(msg, tag)
        except Exception:
            self.exceptions.append(straceback())
            return -2

    def _send_email(self, msg, tag):
        if self.mailto is None:
            return -1

        header = msg.splitlines()
        app = header.append
        app("Submitted on: %s" % time.ctime(self.start_time))
        app("Completed on: %s" % time.asctime())
        app("Elapsed time: %s" % str(self.get_delta_etime()))

        flow = self.flow
        app("Number of errored tasks: %d" % flow.num_errored_tasks)
        app("Number of unconverged tasks: %d" % flow.num_unconverged_tasks)

        # Add the status of the flow.
        strio = StringIO()
        strio.writelines("\n".join(header) + 4 * "\n")
        flow.show_status(stream=strio)

        if self.exceptions:
            # Report the list of exceptions.
            strio.writelines(self.exceptions)

        if tag is None:
            tag = " [ALL OK]" if flow.all_ok else " [WARNING]"

        return sendmail(subject=flow.name + tag, text=strio.getvalue(), mailto=self.mailto)


def sendmail(subject, text, mailto, sender=None):
    """
    Sends an e-mail with unix sendmail.

    Args:
        subject: String with the subject of the mail.
        text: String with the body of the mail.
        mailto: String or list of string with the recipients.
        sender: string with the sender address.
            If sender is None, username@hostname is used.

    Returns:
        Exit status
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

    sendmail = which("sendmail")
    if sendmail is None: return -1
    # msg is string not bytes so must use universal_newlines
    p = Popen([sendmail, "-t"], stdin=PIPE, stderr=PIPE, universal_newlines=True)

    outdata, errdata = p.communicate(msg)
    return len(errdata)


def __test_sendmail():
    retcode = sendmail("sendmail_test", text="hello\nworld", mailto="nobody@nowhere.com")
    print("Retcode", retcode)
    assert retcode == 0
