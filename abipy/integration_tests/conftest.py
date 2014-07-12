"""Configuration file for pytest."""
from __future__ import print_function, division

import pytest
import abipy.abilab as abilab


@pytest.fixture
def fwp(tmpdir):
    """
    Parameters used to initialize Flows.
    This fixture allows us to change the TaskManager to test
    different configurations.
    """
    # Temporary working directory
    fwp.workdir = str(tmpdir)

    # Create the TaskManager.
    fwp.manager = abilab.TaskManager.from_user_config()

    return fwp