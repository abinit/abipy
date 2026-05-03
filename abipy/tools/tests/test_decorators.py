"""Tests for duck module."""

import abipy.tools.decorators as decs
from abipy.core.testing import AbipyTest
from abipy.tools import duck


class DecoratorsTest(AbipyTest):
    def test_return_straceback_ifexc(self):
        """Testing return_straceback_ifexc."""

        def f(a, b):
            return a + b

        with self.assertRaises(TypeError):
            f("hello", 1)

        newf = decs.return_straceback_ifexc(f)
        assert duck.is_string(newf("hello", 1))
