# coding: utf-8
"""Tests for duck module."""
import numpy as np

from abipy.core.testing import AbipyTest
import abipy.tools.decorators as decs
import abipy.tools.duck as duck


class DecoratorsTest(AbipyTest):

    def test_return_straceback_ifexc(self):
        """Testing return_straceback_ifexc."""
        def f(a, b):
            return  a + b

        with self.assertRaises(TypeError):
            f("hello", 1)

        newf = decs.return_straceback_ifexc(f)
        assert duck.is_string(newf("hello", 1))
