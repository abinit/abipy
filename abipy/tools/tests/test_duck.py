# coding: utf-8
"""Tests for duck module."""
from __future__ import division, print_function, absolute_import, unicode_literals

from abipy.core.testing import AbipyTest

import numpy as np

from abipy.tools import duck


class DuckTest(AbipyTest):

    def test_is_string(self):
        """Testing is_string."""
        assert duck.is_string("hello")
        assert not duck.is_string(1)
        assert not duck.is_string([1, 2, 3])

    def test_is_intlike(self):
        """Testing is_intlike."""
        assert duck.is_intlike(1)
        assert duck.is_intlike(1.0)
        assert not duck.is_intlike(1.3)
        assert not duck.is_intlike("1")
        assert not duck.is_intlike([1, 2, 3])

        for t in (np.int32, np.int64, np.float64):
            assert duck.is_intlike(t(123))

        assert not duck.is_intlike(np.complex(123))

    def test_is_listlike(self):
        """Testing is_listlike."""
        #if isinstance(branch, (list, tuple, np.ndarray)):
        assert duck.is_listlike([1, 2, 3])
        assert duck.is_listlike((1,))
        assert not duck.is_listlike({1, 2, 3})
        assert duck.is_listlike(np.zeros(3))
        d = {"hello": 1, "world": 2}
        #assert not duck.is_listlike(d)
        #assert duck.is_listlike(range(1, 3))
