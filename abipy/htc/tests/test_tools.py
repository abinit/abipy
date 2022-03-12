"""
Tests for abipy.htc.tools module
"""

from abipy.core.testing import AbipyTest
from abipy.htc.tools import find_free_port, port_is_open, pid_exists


class TestTools(AbipyTest):

    def test_port_tools(self):
        """Testing port tools."""
        free_port = find_free_port()
        assert free_port > 0
        assert not port_is_open(free_port)

    def test_pid_tools(self):
        """Testing pid tools."""
        import subprocess
        process = subprocess.Popen(["sleep", "120"])
        assert pid_exists(process.pid)
        process.terminate()
        #process.kill()
        process.wait()
        assert not pid_exists(process.pid)