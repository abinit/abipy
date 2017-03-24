"""
This script runs all the python scripts located in this directory
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import unittest

from abipy.core.testing import AbipyTest

root = os.path.dirname(__file__)

class TestPlots(AbipyTest):

    def test_plots(self):
        if not self.has_matplotlib():
            raise unittest.SkipTest("sympy is not installed")

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
            for i, e in enumerate(errors):
                print(80 * "*")
                print(i, e)
                print(80 * "*")
        assert not errors
