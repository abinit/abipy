# coding: utf-8

from abipy.core.testing import AbipyTest
from abipy.flowtk.launcher import ScriptEditor, PyFlowScheduler, MultiFlowScheduler


def test_script_editor():
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
        """Testing PyFlowScheduler API."""

        assert "weeks:" in PyFlowScheduler.autodoc()

        with self.assertRaises(PyFlowScheduler.Error):
            PyFlowScheduler()
        with self.assertRaises(PyFlowScheduler.Error):
            PyFlowScheduler(foo="bar")

        sched = PyFlowScheduler(seconds=2)
        assert str(sched)
        assert sched.sched_options.seconds == 2
        assert sched.flow is None
        assert int(sched.pid) > 0
        assert sched.num_excs == 0
        assert not sched.rmflow

        #sched.start()
        #assert sched.get_delta_etime()

    def test_multiflowscheduler_api(self):

        sched = MultiFlowScheduler(seconds=2, sqldb_path="foobar.db")
        assert str(sched)
        assert sched.sched_options.seconds == 2
        assert not sched.get_incoming_flows()
        assert not sched.get_incoming_flows()
        #assert sched.get_flow_status_by_id(1) == (None, None)
        assert not sched.groupby_status()

        #import threading
        #thread threading.Thread(target=sched.start, demon=False)
        #thread.start()
        #assert sched.get_delta_etime()
        #assert thread.is_alive()

