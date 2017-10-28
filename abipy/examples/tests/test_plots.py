"""
This script runs all the python scripts located in this directory
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import unittest

from abipy.core.testing import AbipyTest

root = os.path.dirname(__file__)

class TestPlots(AbipyTest):

    #def test_plots_with_runtests(self):
    #    """Running _runplots script."""
    #    if not self.has_matplotlib():
    #        raise unittest.SkipTest("matplotlib is not installed")
    #    script = os.path.abspath(os.path.join(root, "..", "_runplots.py"))
    #    assert os.path.exists(script)
    #    assert os.system(script + " --backend Agg")

    def test_plots_with_exec(self):
        """Running plot script with exec."""
        if os.environ.get("TRAVIS"):
            raise unittest.SkipTest("Skipping plot examples on TRAVIS")
        # Travis issue
        """
ERROR: Failure: SyntaxError (unqualified exec is not allowed in function 'test_plots_with_exec' because it contains a nested function with free variables (test_plots.py, line 51))
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/travis/miniconda/envs/test-environment/lib/python2.7/site-packages/nose/loader.py", line 418, in loadTestsFromName
    addr.filename, addr.module)
  File "/home/travis/miniconda/envs/test-environment/lib/python2.7/site-packages/nose/importer.py", line 47, in importFromPath
    return self.importFromDir(dir_path, fqname)
  File "/home/travis/miniconda/envs/test-environment/lib/python2.7/site-packages/nose/importer.py", line 94, in importFromDir
    mod = load_module(part_fqname, fh, filename, desc)
SyntaxError: unqualified exec is not allowed in function 'test_plots_with_exec' because it contains a nested function with free variables (test_plots.py, line 51)
"""

        if not self.has_matplotlib():
            raise unittest.SkipTest("matplotlib is not installed")

        #0 file /home/travis/build/abinit/abipy/abipy/examples/tests/../plot/plot_qpconverge.py
        # Traceback (most recent call last):
        #  File "/home/travis/build/abinit/abipy/abipy/examples/tests/test_plots.py", line 29, in test_plots
        #    exec(fh.read())
        #  File "<string>", line 16, in <module>
        #  File "<string>", line 16, in <listcomp>
        #NameError: name 'abidata' is not defined
        if sys.version[0:3] >= '3.4':
            raise unittest.SkipTest("Weird behaviour of exec when Python version >= 3.4")

        import matplotlib.pyplot as plt

        plot_dir = os.path.join(root, "..", "plot")
        count, errors = 0, []
        for fname in os.listdir(plot_dir):
            if not (fname.endswith(".py") and fname.startswith("plot_")): continue
            count += 1
            path = os.path.join(plot_dir, fname)
            try:
                with open(path, "rt") as fh:
                    exec(fh.read())
            except Exception:
                errors.append("file %s\n %s" % (path, self.straceback()))
            plt.close("all")

        print("Tested ", count, "plot scripts")
        assert count > 0
        if errors:
            #for i, e in enumerate(errors):
            #    print(80 * "*")
            #    print(i, e)
            #    print(80 * "*")
            raise RuntimeError("\n\n".join(str(e) for e in errors))
