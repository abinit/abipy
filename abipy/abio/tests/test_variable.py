"""Tests for variable module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.abio.variable import InputVariable


class TestInputVariable(AbipyTest):

    def test_inputvariable(self):
        """Testing InputVariable."""
        v = InputVariable(name="ecut", value=5)
        assert v.name == "ecut"
        assert not v.units
        assert str(v) == " ecut 5"
