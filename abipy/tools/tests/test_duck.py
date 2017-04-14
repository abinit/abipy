# coding: utf-8
"""Tests for duck module."""
from __future__ import division, print_function, absolute_import, unicode_literals

from abipy.core.testing import AbipyTest

import numpy as np

from abipy.tools import duck


class DuckTest(AbipyTest):

    def test_intlike(self):
        """Testing is_intlike."""
        assert duck.is_intlike(1)
        assert duck.is_intlike(1.0)
        assert not duck.is_intlike(1.3)
        assert not duck.is_intlike("1")
        assert not duck.is_intlike([1, 2, 3])

        for t in (np.int32, np.int64, np.float64):
            assert duck.is_intlike(t(123))

        assert not duck.is_intlike(np.complex(123))

