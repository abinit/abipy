# coding: utf-8
"""Test manager files."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest


class ManagerTest(AbipyTest):

    def test_managers(self):
        """Trying to read all managers files in abipy/data/managers."""
        root = os.path.join(abidata.dirpath, "managers")
        yaml_paths = [os.path.join(root, f) for f in os.listdir(root) if f.endswith(".yml") and "_manager" in f]
        assert yaml_paths
        for p in yaml_paths:
            manager = flowtk.TaskManager.from_file(p)
            print(manager)
            shell = manager.to_shell_manager(mpi_procs=2)
        #assert 0

    def test_schedulers(self):
        """Trying to read all scheduler files in abipy/data/managers."""
        root = os.path.join(abidata.dirpath, "managers")
        yaml_paths = [os.path.join(root, f) for f in os.listdir(root) if f.endswith(".yml") and "_scheduler" in f]
        assert yaml_paths
        for p in yaml_paths:
            sched = flowtk.PyFlowScheduler.from_file(p)
            print(sched)

        #assert 0