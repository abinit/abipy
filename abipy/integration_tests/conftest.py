"""Configuration file for pytest."""
from __future__ import print_function, division

import pytest
import abipy.abilab as abilab


# TaskManager configurations
manager_confs = [
    """
    qtype: shell
    mpi_runner: mpirun
    pre_run:
        - "source ~/Coding/Abinit/bzr_archives/env.sh"
    policy:
        autoparal: 0
        max_ncpus: 1
    """,
]


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
