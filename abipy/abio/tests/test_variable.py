"""Tests for variable module."""

from abipy.abio.variable import InputVariable
from abipy.core.testing import AbipyTest


class TestInputVariable(AbipyTest):
    def test_inputvariable(self):
        """Testing InputVariable."""
        v = InputVariable(name="ecut", value=5)
        assert v.name == "ecut"
        assert not v.units
        assert str(v) == " ecut 5"
