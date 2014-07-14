"""Configuration file for pytest."""
from __future__ import print_function, division

import pytest
import yaml
import copy
import abipy.abilab as abilab

# Create the list of TaskManager configurations used for the integration tests.
# Note that the items in manager_confs must be hashable hence we cannot use dictionaries directly
# To bypass this problem we operate on dictionaries to generate the different configuration
# and then we convert the dictionaries to string with yaml.dump. This string will be passed
# to TasksManager.from_string in fwp. base_conf looks like:

#manager_confs = [
#    """
#    qtype: shell
#    mpi_runner: mpirun
#    pre_run:
#        - "source ~/Coding/Abinit/bzr_archives/env.sh"
#    policy:
#        autoparal: 0
#        max_ncpus: 1
#    """,
#]

# Read the base configuration from file
with open("taskmanager.yaml") as fh:
    base_conf = yaml.load(fh)

# Build list of configurations.
manager_confs = []

for autoparal in [0]: #, 1]:
    max_ncpus = 1 if autoparal == 0 else 2

    newd = copy.deepcopy(base_conf)
    newd["policy"]["autoparal"] = autoparal
    newd["policy"]["max_ncpus"] = max_ncpus
    manager_confs.append(newd)

manager_confs = [yaml.dump(d) for d in manager_confs]
#manager_confs = [base_conf]


@pytest.fixture(params=manager_confs)
def fwp(tmpdir, request):
    """
    Parameters used to initialize Flows.
    This fixture allows one to change the TaskManager
    and therefore we can easily test different configurations.
    """
    # Temporary working directory
    fwp.workdir = str(tmpdir)

    # Create the TaskManager.
    fwp.manager = abilab.TaskManager.from_string(request.param)

    return fwp


def pytest_report_header(config):
    """Write the initial header."""
    lines = ["\n*** Integration tests for abipy+abinit+pymatgen ***\n"]
    app = lines.append

    app("Assuming the enviroment is properly configured:")
    app("In particular, we assume that abinit is in $PATH and can be executed.")
    app("Change taskmanager.yaml according to your platform.")
    app("Number of tasksmanager configurations used: %d" % len(manager_confs))

    if config.option.verbose > 0:
        for i, s in enumerate(manager_confs):
            app(80 * "=")
            app("TaskManager #%d" % i)
            app(s)
            app(80 * "=")

    app("")
    return lines


#def pytest_runtest_logreport(report):
#    """Reporting hook"""
#    if report.outcome == "failed":
#        print("noedid", report.nodeid, "failed")

