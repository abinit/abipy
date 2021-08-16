# coding: utf-8

import unittest

from abipy.core.testing import AbipyTest
from abipy.flowtk.launcher import ScriptEditor, PyFlowScheduler


class ScriptEditorTest(AbipyTest):

    def test_api(self):
        """base tests for ScriptEditor"""
        se = ScriptEditor()
        se.shebang()
        se.declare_var("FOO", "BAR")
        se.add_emptyline()
        se.add_comment("This is a comment")
        se.declare_vars({"FOO1": "BAR1"})
        se.load_modules(["module1", "module2"])
        s = se.get_script_str()
        assert "export FOO=BAR" in s
        assert "export FOO1=BAR1" in s
        assert "module load module1" in s


class PyFlowSchedulerTest(AbipyTest):

    def test_pyflowscheduler_api(self):

        assert "weeks:" in PyFlowScheduler.autodoc()

        with self.assertRaises(PyFlowScheduler.Error):
            PyFlowScheduler()
        with self.assertRaises(PyFlowScheduler.Error):
            PyFlowScheduler(foo="bar")

        sched = PyFlowScheduler(seconds=2)
        assert str(sched)
        assert sched.sched_options.seconds == 2
        assert sched.flow is None
        assert sched.num_excs == 0
        assert not sched.rmflow
        assert sched.get_delta_etime()


