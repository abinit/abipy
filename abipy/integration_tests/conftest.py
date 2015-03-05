"""Configuration file for pytest."""
from __future__ import print_function, division, unicode_literals

import os
import pytest
import yaml
import copy
import abipy.abilab as abilab

from monty.collections import AttrDict

# Create the list of TaskManager configurations used for the integration tests.
# Note that the items in _manager_confs must be hashable hence we cannot use dictionaries directly
# To bypass this problem we operate on dictionaries to generate the different configuration
# and then we convert the dictionary to string with yaml.dump. This string will be passed
# to TasksManager.from_string in fwp. base_conf looks like:


# Read the base configuration from file
with open(os.path.join(os.path.dirname(__file__), "manager.yml")) as fh:
    base_conf = yaml.load(fh)

# Build list of configurations.
_manager_confs = []

for autoparal in [1]: #, 1]:
    newd = copy.deepcopy(base_conf)
    newd["policy"]["autoparal"] = autoparal
    _manager_confs.append(newd)

_manager_confs = [yaml.dump(d) for d in _manager_confs]
#_manager_confs = [base_conf]


@pytest.fixture(params=_manager_confs)
def fwp(tmpdir, request):
    """
    Parameters used to initialize Flows.

    This fixture allows us to change the :class:`TaskManager`
    so that we can easily test different configurations.
    """
    # Temporary working directory
    fwp.workdir = str(tmpdir)

    # Create the TaskManager.
    fwp.manager = abilab.TaskManager.from_string(request.param)

    fwp.scheduler = abilab.PyFlowScheduler.from_file(os.path.join(os.path.dirname(__file__), "scheduler.yml"))

    return fwp


# Use tuples instead of dict because pytest require objects to be hashable.
_tvars_list = [
    (("paral_kgb", 0),),
    (("paral_kgb", 1),),
]

@pytest.fixture(params=_tvars_list)
def tvars(request):
    """
    Abinit variables passed to the test functions.

    This fixture allows us change the variables in the input files
    so that we can easily test different scenarios e.g. runs with or without
    paral_kgb==1
    """
    return AttrDict({k: v for k, v in request.param})


def pytest_addoption(parser):
    """Add extra command line options."""
    parser.addoption('--loglevel', default="ERROR", type=str,
                     help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")


def pytest_report_header(config):
    """Write the initial header."""
    lines = ["\n*** Integration tests for abipy + abinit + pymatgen ***\n"]
    app = lines.append

    app("Assuming the enviroment is properly configured:")
    app("In particular, we assume that the abinit executable is in $PATH and can be executed.")
    app("Change taskmanager.yml according to your platform.")
    app("Number of taskmanager configurations: %d" % len(_manager_confs))

    if config.option.verbose > 0:
        for i, s in enumerate(_manager_confs):
            app(80 * "=")
            app("TaskManager #%d" % i)
            app(s)
            app(80 * "=")

    app("")

    # Initialize logging
    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, config.option.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % config.option.loglevel)
    logging.basicConfig(level=numeric_level)

    return lines


#def pytest_runtest_logreport(report):
#    """Reporting hook"""
#    if report.outcome == "failed":
#        print("noedid", report.nodeid, "failed")

