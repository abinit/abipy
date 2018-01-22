# coding: utf-8
"""Tests for tools __init__.py module."""
from __future__ import division, print_function, absolute_import, unicode_literals

import abipy.tools as at
from abipy.core.testing import AbipyTest


class TestTools(AbipyTest):

    def test_getattrd(self):
        """Testing getattrd."""
        with self.assertRaises(AttributeError): at.getattrd(int, 'a')
        assert at.getattrd(int, 'a', default=None) is None
        assert at.getattrd(int, '__class__.__name__') == "type"
        assert at.hasattrd(int, '__class__.__name__')
        assert not at.hasattrd(int, 'foobar.__name__')
