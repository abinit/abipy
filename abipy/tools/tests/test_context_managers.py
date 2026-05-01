"""Tests for context_managers module."""

import time

from abipy.core.testing import AbipyTest
from abipy.tools.context_managers import Timeout, Timer, temporary_change_attributes


class TestContextManagers(AbipyTest):
    """Test context_managers."""

    def test_api(self):
        """Testing tonumber"""

        class Something:
            def __init__(self, x, y):
                self.x = x
                self.y = y

        with Timer("Creating Something"):
            s = Something(1, 2)

        assert s.x == 1 and s.y == 2

        with temporary_change_attributes(s, x=4, y=5):
            assert s.x == 4 and s.y == 5

        assert s.x == 1 and s.y == 2

        with self.assertRaises(ValueError), Timeout(seconds=0, message="Timeout"):
            pass

        with self.assertRaises(TimeoutError), Timeout(seconds=2, message="Timeout"):
            time.sleep(6)
