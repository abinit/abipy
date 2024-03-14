"""Tests for enums module"""

import pytest

from abipy.core.testing import AbipyTest
from abipy.abio.enums import GWR_TASK, WFK_TASK, RUNL


class TestEnums(AbipyTest):

    def test_api(self):
        assert RUNL.GSTATE == 0
        assert str(RUNL.GSTATE) == "0"
        assert RUNL.validate(0) is None
        with pytest.raises(ValueError):
            RUNL.validate(-666)
        assert RUNL.GWR == 6

        assert GWR_TASK.HDIAGO == "HDIAGO"
        assert GWR_TASK.HDIAGO != "HDIAG"
        assert GWR_TASK.validate("HDIAGO_FULL") is None
        with pytest.raises(ValueError):
            GWR_TASK.validate("foobar")
        with pytest.raises(ValueError):
            GWR_TASK.validate(1)
