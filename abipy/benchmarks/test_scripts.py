#!/usr/bin/env python
"""
This script runs all the python scripts located in this directory
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest
from abipy.benchmarks import build_bench_main_parser, bench_monkey_patch_options


root = os.path.abspath(os.path.join(os.path.dirname(__file__)))


class TestScripts(AbipyTest):

    #def test_all_scripts(self):
    #    """Testing all scripts in abipy/benckmarks"""
    #    from subprocess import call
    #    retcode = call(os.path.join(root, "_runemall.py"))
    #    if retcode != 0: print(retcode, "scripts in ~abipy/data/runs exited with non-zero return code")
    #    assert retcode == 0

    def test_build_flow_method_in_scripts(self):
        """Testing build_flow method in all scripts in abipy/benchmarks"""
        parser = build_bench_main_parser()

        import importlib
        import tempfile
        count, errors = 0, []
        for fname in os.listdir(root):
            if not fname.endswith(".py") or fname.startswith("_") or fname.startswith("test_"):
                continue
            count += 1
            s = "abipy.benchmarks." + fname.replace(".py", "")
            module = importlib.import_module(s)
            # flow will be produced in a temporary workdir.
            workdir = tempfile.mkdtemp(prefix='bench_' + os.path.basename(fname))
            options = parser.parse_args(["--workdir", workdir])
            bench_monkey_patch_options(options)
            # Instantiate the manager.
            options.manager = flowtk.TaskManager.as_manager(options.manager)
            try:
                flow = module.build_flow(options)
                assert flow is not None
                flow.build_and_pickle_dump()
            except Exception:
                errors.append("file %s\n %s" % (s, self.straceback()))

        print("Tested ", count, "scripts")
        assert count > 0
        if errors:
            for i, e in enumerate(errors):
                print(80 * "*")
                print(i, e)
                print(80 * "*")
        assert not errors
