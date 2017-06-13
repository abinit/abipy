"""Tests for core.restapi module"""
from __future__ import print_function, division, unicode_literals, absolute_import

from abipy.core.testing import AbipyTest
from abipy.core import restapi


class TestMpRestApi(AbipyTest):
    """Test interfaces with the Materials Project REST API."""
    def test_mprester(self):
        with restapi.get_mprester() as rest:
            print(rest.Error)
