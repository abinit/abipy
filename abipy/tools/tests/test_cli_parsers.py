"""Tests for cli_parsers module."""

import pytest

from abipy.core.testing import AbipyTest
from abipy.tools import cli_parsers as cli


class TestCliParsers(AbipyTest):
    def test_range_from_str(self):
        assert cli.range_from_str(None) == None
        assert cli.range_from_str("3") == range(3)
        assert cli.range_from_str("2:3") == range(2, 3)
        assert cli.range_from_str("2:3:2") == range(2, 3, 2)
        with pytest.raises(ValueError):
            cli.range_from_str("2:3:2:3")
