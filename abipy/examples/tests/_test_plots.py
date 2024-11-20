"""
This script runs all the python scripts located in this directory
"""
import sys
import os
import unittest

from abipy.core.testing import AbipyTest

root = os.path.dirname(__file__)

class TestPlots(AbipyTest):

    def test_plots_with_exec(self):
        """
        Running plot scripts in example/plots directory with exec.
        """
        # Travis issue
        #if os.environ.get("TRAVIS"):
        #    raise unittest.SkipTest("Skipping plot examples on TRAVIS")

        if not self.has_matplotlib():
            raise unittest.SkipTest("matplotlib is not installed")

        import matplotlib.pyplot as plt
        from abipy.tools.plotting import set_plotly_default_show
        ply_show = False
        ply_show = True
        print("Setting plotly_default_show to: ", ply_show)
        set_plotly_default_show(ply_show)

        plot_dir = os.path.join(root, "..", "plot")
        count, errors = 0, []
        for fname in os.listdir(plot_dir):
            if not (fname.endswith(".py") and fname.startswith("plot_")): continue
            count += 1
            path = os.path.join(plot_dir, fname)
            print("About to execute:", path)
            try:
                with open(path, "rt") as fh:
                    exec(fh.read(), {}, {})
            except Exception:
                errors.append("file %s\n %s" % (path, self.straceback()))
            plt.close("all")

        print("Tested ", count, "plot scripts")
        assert count > 0
        if errors:
            raise RuntimeError("\n\n".join(str(e) for e in errors))
