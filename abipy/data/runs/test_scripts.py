#!/usr/bin/env python
"""
This script runs all the python scripts located in this directory
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import tempfile

from subprocess import call
from abipy.core.testing import AbipyTest
import abipy.flowtk as flowtk


root = os.path.dirname(__file__)


class TestScripts(AbipyTest):
    #def test_all_scripts(self):
    #    """Testing all scripts in abipy/data/runs"""
    #    root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    #    retcode = call(os.path.join(root, "_runemall.py"))
    #    if retcode != 0: print(retcode, "scripts in ~abipy/data/runs exited with non-zero return code")
    #    assert retcode == 0

    def test_build_flow_method_in_scripts(self):
        """Testing build_flow method in all scripts in abipy/data/runs"""
        parser = flowtk.build_flow_main_parser()

        import importlib
        count, errors = 0, []
        for fname in os.listdir(root):
            if not (fname.endswith(".py") and fname.startswith("run_")): continue
            count += 1
            s = "abipy.data.runs." + fname.replace(".py", "")
            module = importlib.import_module(s)
            # flow will be produced in a temporary workdir.
            workdir = tempfile.mkdtemp(prefix='flow_' + os.path.basename(fname))
            options = parser.parse_args(["--workdir", workdir])
            # Instantiate the manager.
            options.manager = flowtk.TaskManager.as_manager(options.manager)
            try:
                flow = module.build_flow(options)
                assert flow is not None
                flow.build_and_pickle_dump()
            except Exception as exc:
                errors.append(str(exc))

        print("Tested ", count, "scripts")
        assert not errors
        assert count > 0
