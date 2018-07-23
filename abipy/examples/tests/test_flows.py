"""
This script runs all the flows in the `flow` directory.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import tempfile
import warnings

from subprocess import call
from abipy.core.testing import AbipyTest
import abipy.flowtk as flowtk

root = os.path.join(os.path.dirname(__file__), "..", "flows")


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
            print("Testing:", fname)
            count += 1
            s = "abipy.examples.flows." + fname.replace(".py", "")
            module = importlib.import_module(s)

            if hasattr(module, "exclude_py_versions"):
                if sys.version[0:3] in module.exclude_py_versions:
                    warnings.warn("%s excludes python versions %s" % (s, str(module.exclude_py_versions)))
                    continue

            # flow will be produced in a temporary workdir.
            workdir = tempfile.mkdtemp(prefix='flow_' + os.path.basename(fname))
            options = parser.parse_args(["--workdir", workdir])
            # Instantiate the manager.
            options.manager = flowtk.TaskManager.as_manager(options.manager)

            # Check if flow has requirements on the Abinit version.
            if hasattr(module, "minimum_abinit_version"):
                if not options.manager.abinit_build.version_ge(module.minimum_abinit_version):
                    warnings.warn("%s requires %s but Abinit version: %s" %
                          (s, module.minimum_abinit_version, options.manager.abinit_build.version))
                    continue

            try:
                flow = module.build_flow(options)
                assert flow is not None
                flow.build_and_pickle_dump()
                flow.show_status()
                #flow.make_scheduler().start()
            except Exception:
                errors.append("file %s\n %s" % (s, self.straceback()))

        print("Tested ", count, "scripts")
        assert count > 0
        if errors:
            for i, e in enumerate(errors):
                print(80 * "*", file=sys.stderr)
                print(i, e, file=sys.stderr)
                print(80 * "*", file=sys.stderr)

        assert not errors