# coding: utf-8
"""Tests for textools module."""
from abipy.tools.text import *
from abipy.core.testing import AbipyTest


class TestTools(AbipyTest):
    """Test texttools."""

    def test_tonumber(self):
        """Testing tonumber"""
        assert tonumber('123') == 123
        assert tonumber('1.23') == 1.23
        assert tonumber('12D3') == 12000
        assert tonumber('12E-3') == 0.012

    def test_nums_and_text(self):
        """Testing nums_and_test"""
        line = "   intxc =         0  ionmov =         0    iscf =         7 xclevel =         2"
        numbers = [0.0, 0.0, 7.0, 2.0]
        text = ' intxc = ionmov = iscf = xclevel ='
        self.assertEqual(nums_and_text(line), (numbers, text.replace(", u'", ", '")))
        line = "WF disk file :     45.918 Mbytes ; DEN or POT disk file :      0.490 Mbytes."
        numbers = [45.918, 0.49]
        text = ' WF disk file : Mbytes ; DEN or POT disk file : Mbytes.'
        self.assertEqual(nums_and_text(line), (numbers, text.replace(", u'", ", '")))
        line = "            acell      1.0000000000d+00  1.0000000000D+00  1.0000000000E+00 Bohr"
        numbers = [1.0, 1.0, 1.0]
        text = ' acell Bohr'
        self.assertEqual(nums_and_text(line), (numbers, text.replace(", u'", ", '")))
