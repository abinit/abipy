"""Configuration file for pytest."""

from __future__ import annotations

import copy
import io
import os
from types import SimpleNamespace

import pytest
from monty.collections import AttrDict
from monty.string import marquee
from ruamel.yaml import YAML

from abipy import flowtk
from abipy.flowtk.qadapters import ShellAdapter, kill_process_group
from abipy.integration_tests.helpers import collect_task_errors
from abipy.tools.iotools import yaml_safe_load_path

# Are we running on travis?

# Create the list of Manager configurations used for the integration tests.
# Note that the items in _manager_confs must be hashable hence we cannot use dictionaries directly
# To bypass this problem we operate on dictionaries to generate the different configuration
# and then we convert the dictionary to string with yaml.dump. This string will be passed
# to Manager.from_string in fwp. base_conf looks like:

USER_CONFIG_DIR = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")
# USER_CONFIG_DIR = os.path.dirname(__file__)

# Read the base configuration from file
# with open(os.path.join(USER_CONFIG_DIR, "manager.yml")) as fh:
#    base_conf = yaml.safe_load(fh)

base_conf = yaml_safe_load_path(os.path.join(USER_CONFIG_DIR, "manager.yml"))


def _to_yaml(d: dict) -> str:
    """
    Convert dictionary to YAML string in block style.
    Aliases are disabled so that the anchors of manager.yml (e.g. &hardware) are expanded.
    """
    yaml = YAML(typ="safe", pure=True)
    yaml.default_flow_style = False
    yaml.width = 120
    yaml.representer.ignore_aliases = lambda *args: True
    stream = io.StringIO()
    yaml.dump(d, stream)
    return stream.getvalue()


# Build list of configurations.
# By default, autoparal is disabled and the tasks are executed with a fixed number of MPI processes
# so that the parallel distribution does not depend on the machine (autoparal may select values of
# npfft that are not compatible with the FFT mesh or lead to deadlocks).
# Use ABIPY_ITEST_MPI_PROCS to change the number of MPI processes and
# ABIPY_ITEST_AUTOPARAL=1 to run the tests with autoparal as well.
ITEST_MPI_PROCS = int(os.environ.get("ABIPY_ITEST_MPI_PROCS", 2))
_autoparal_list = [0, 1] if os.environ.get("ABIPY_ITEST_AUTOPARAL", "0") == "1" else [0]

_manager_confs, _manager_ids = [], []

for autoparal in _autoparal_list:
    newd = copy.deepcopy(base_conf)
    if "policy" not in newd:
        newd["policy"] = {}
    newd["policy"]["autoparal"] = autoparal
    # A task is set to S_ERROR if the output file has not been modified for frozen_timeout.
    # The default (1 hour) is too large for the small systems used in the integration tests.
    # Keep it below the pytest timeout so that the scheduler detects the problem first.
    newd["policy"].setdefault("frozen_timeout", "0-0:2:0")
    # mpi_procs is used only if autoparal is disabled.
    mpi_procs = None if autoparal else ITEST_MPI_PROCS
    _manager_confs.append((_to_yaml(newd), mpi_procs))
    # Short id used in the test names e.g. itest_foo[mpi2-paral_kgb1].
    _manager_ids.append("autoparal" if autoparal else f"mpi{mpi_procs}")


# Options of the scheduler that override the values given in the user configuration file.
# Tolerate a few python exceptions (e.g. netcdf file read while Abinit is still writing it).
# The exceptions are still reported as warnings by assert_flow_ok.
SCHEDULER_OVERRIDES = dict(max_num_pyexcs=2)


def _apply_scheduler_overrides(sched):
    for k, v in SCHEDULER_OVERRIDES.items():
        setattr(sched, k, v)
    return sched


@pytest.fixture(params=_manager_confs, ids=_manager_ids)
def fwp(tmpdir, request):
    """
    Parameters used to initialize Flows.

    This fixture allows us to change the |TaskManager| so that we can easily test different configurations.
    """
    conf, mpi_procs = request.param
    manager = flowtk.TaskManager.from_string(conf)
    if mpi_procs is not None:
        manager = manager.new_with_fixed_mpi_omp(mpi_procs, 1)

    return SimpleNamespace(
        # Temporary working directory
        workdir=str(tmpdir),
        manager=manager,
        scheduler=_apply_scheduler_overrides(
            flowtk.PyFlowScheduler.from_file(os.path.join(USER_CONFIG_DIR, "scheduler.yml"))
        ),
        abinit_build=flowtk.AbinitBuild(),
    )


# Variables defining the MPI distribution when paral_kgb == 1.
_NP_VARS = ("np_spkpt", "npkpt", "npband", "npfft", "npspinor")


def _add_abinit_autoparal(abinit_input, input_string: str) -> str:
    """
    When the abipy autoparal is disabled, paral_kgb 1 inputs do not specify the MPI distribution
    and Abinit aborts if nprocs > 1. In this case we add the Abinit variable autoparal 1
    so that Abinit selects np_spkpt, npband, npfft at runtime.
    The distribution depends only on the number of MPI processes hence the tests are reproducible.
    """
    if abinit_input.get("paral_kgb", 0) != 1 or "autoparal" in abinit_input:
        return input_string
    # Abinit accepts autoparal only for GS and DFPT (e.g. GW has its own MPI distribution).
    if abinit_input.get("optdriver", 0) not in (0, 1):
        return input_string
    if any(v in abinit_input for v in _NP_VARS):
        return input_string
    return input_string + "\n# Added by the integration tests (conftest.py)\nautoparal 1\n"


@pytest.fixture(autouse=True)
def _itest_guard(monkeypatch):
    """
    - Apply SCHEDULER_OVERRIDES to the schedulers created with flow.make_scheduler().
    - Keep track of the processes launched by the ShellAdapter and kill the ones
      that are still running at the end of the test (e.g. after a timeout or a failure).
    - Add autoparal 1 to the Abinit input files with paral_kgb 1 and without np* variables,
      see _add_abinit_autoparal.
    """
    orig_make_input = flowtk.AbinitTask.make_input

    def make_input(self, *args, **kwargs):
        s = orig_make_input(self, *args, **kwargs)
        return _add_abinit_autoparal(self.input, s)

    monkeypatch.setattr(flowtk.AbinitTask, "make_input", make_input)

    orig_from_user_config = flowtk.PyFlowScheduler.from_user_config.__func__

    def from_user_config(cls):
        return _apply_scheduler_overrides(orig_from_user_config(cls))

    monkeypatch.setattr(flowtk.PyFlowScheduler, "from_user_config", classmethod(from_user_config))

    processes = []
    orig_submit = ShellAdapter._submit_to_queue

    def _submit_to_queue(self, script_file):
        results = orig_submit(self, script_file)
        processes.append(results.process)
        return results

    monkeypatch.setattr(ShellAdapter, "_submit_to_queue", _submit_to_queue)

    yield

    for process in processes:
        if process.poll() is None:
            print(f"Killing leftover process {process.pid} ({process.args})")
            kill_process_group(process.pid)
            process.wait()


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Add the flow workdir and the tail of the files of the failed Abinit tasks to the report."""
    outcome = yield
    report = outcome.get_result()
    if report.when != "call" or not report.failed:
        return

    workdir = getattr(item.funcargs.get("fwp"), "workdir", None)
    if workdir is None:
        return

    report.sections.append(("Flow workdir", workdir))
    errors = collect_task_errors(workdir)
    if errors:
        report.sections.append(("Abinit task errors", errors))


# Use tuples instead of dict because pytest require objects to be hashable.
_tvars_list = [
    # (("paral_kgb", 0),),
    (("paral_kgb", 1),),
]


@pytest.fixture(params=_tvars_list, ids=["-".join(f"{k}{v}" for k, v in t) for t in _tvars_list])
def tvars(request):
    """
    Abinit variables passed to the test functions.

    This fixture allows us change the variables in the input files
    so that we can easily test different scenarios e.g. runs with or without paral_kgb == 1
    """
    return AttrDict({k: v for k, v in request.param})


def pytest_addoption(parser):
    """Add extra command line options."""
    parser.addoption(
        "--loglevel",
        default="ERROR",
        type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG",
    )
    # parser.addoption('--manager', default=None, help="TaskManager file (defaults to the manager.yml found in cwd"


def pytest_report_header(config):
    """Write the initial header."""
    lines = []
    app = lines.append
    app("\n" + marquee("Begin integration tests for AbiPy + abinit", mark="="))
    app("\tAssuming the environment is properly configured:")
    app("\tIn particular, the abinit executable must be in $PATH.")
    app("\tChange manager.yml according to your platform.")
    app("\tNumber of TaskManager configurations used: %d" % len(_manager_confs))
    if not config.pluginmanager.hasplugin("timeout"):
        app("\tWARNING: pytest-timeout is not installed. Tests may hang forever if Abinit hangs.")

    if config.option.verbose > 0:
        for mid, (s, mpi_procs) in zip(_manager_ids, _manager_confs, strict=True):
            app(80 * "=")
            app(f"TaskManager [{mid}]")
            if mpi_procs is not None:
                app(f"autoparal disabled, tasks executed with mpi_procs: {mpi_procs}")
            app(80 * "=")
            app(s)
    app("")

    # Print info on Abinit build
    abinit_build = flowtk.AbinitBuild()
    print()
    print(abinit_build)
    print()
    if not config.option.verbose:
        print("Use --verbose for additional info")
    else:
        print(abinit_build.info)

    # Initialize logging
    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging

    numeric_level = getattr(logging, config.option.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % config.option.loglevel)
    logging.basicConfig(level=numeric_level)

    return lines
