"""Configuration file for pytest."""
import os
import pytest
import yaml
import copy
import abipy.flowtk as flowtk

from monty.collections import AttrDict
from monty.string import marquee

# Are we running on travis?
on_travis = os.environ.get("TRAVIS", "False") == "true"

# Create the list of Manager configurations used for the integration tests.
# Note that the items in _manager_confs must be hashable hence we cannot use dictionaries directly
# To bypass this problem we operate on dictionaries to generate the different configuration
# and then we convert the dictionary to string with yaml.dump. This string will be passed
# to Manager.from_string in fwp. base_conf looks like:

USER_CONFIG_DIR = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")
#USER_CONFIG_DIR = os.path.dirname(__file__)

# Read the base configuration from file
with open(os.path.join(USER_CONFIG_DIR, "manager.yml")) as fh:
    base_conf = yaml.safe_load(fh)

# Build list of configurations.
_manager_confs = []

for autoparal in [1]: #, 1]:
    newd = copy.deepcopy(base_conf)
    if "policy" not in newd: newd["policy"] = {}
    newd["policy"]["autoparal"] = autoparal
    _manager_confs.append(newd)

_manager_confs = [yaml.dump(d) for d in _manager_confs]
#_manager_confs = [base_conf]


@pytest.fixture(params=_manager_confs)
def fwp(tmpdir, request):
    """
    Parameters used to initialize Flows.

    This fixture allows us to change the |TaskManager| so that we can easily test different configurations.
    """
    # Temporary working directory
    fwp.workdir = str(tmpdir)

    # Create the TaskManager.
    fwp.manager = flowtk.TaskManager.from_string(request.param)
    fwp.scheduler = flowtk.PyFlowScheduler.from_file(os.path.join(USER_CONFIG_DIR, "scheduler.yml"))
    fwp.on_travis = on_travis
    fwp.abinit_build = flowtk.AbinitBuild()

    return fwp


# Use tuples instead of dict because pytest require objects to be hashable.
if on_travis:
    _tvars_list = [
        (("paral_kgb", 0),),
        #(("paral_kgb", 1),),
    ]

else:
    _tvars_list = [
        #(("paral_kgb", 0),),
        (("paral_kgb", 1),),
    ]


@pytest.fixture(params=_tvars_list)
def tvars(request):
    """
    Abinit variables passed to the test functions.

    This fixture allows us change the variables in the input files
    so that we can easily test different scenarios e.g. runs with or without paral_kgb == 1
    """
    return AttrDict({k: v for k, v in request.param})


def pytest_addoption(parser):
    """Add extra command line options."""
    parser.addoption('--loglevel', default="ERROR", type=str,
                     help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")
    #parser.addoption('--manager', default=None, help="TaskManager file (defaults to the manager.yml found in cwd"


def pytest_report_header(config):
    """Write the initial header."""
    lines = []
    app = lines.append
    app("\n" + marquee("Begin integration tests for AbiPy + abinit", mark="="))
    app("\tAssuming the environment is properly configured:")
    app("\tIn particular, the abinit executable must be in $PATH.")
    app("\tChange manager.yml according to your platform.")
    app("\tNumber of TaskManager configurations used: %d" % len(_manager_confs))

    if config.option.verbose > 0:
        for i, s in enumerate(_manager_confs):
            app(80 * "=")
            app("TaskManager #%d" % i)
            app(s)
            app(80 * "=")
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
        raise ValueError('Invalid log level: %s' % config.option.loglevel)
    logging.basicConfig(level=numeric_level)

    return lines
