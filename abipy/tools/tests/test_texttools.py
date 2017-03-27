from __future__ import print_function, division

from abipy.tools.text import *
from abipy.core.testing import *


class TestTools(AbipyTest):
    """Test texttools."""
    def test_tonumber(self):
        stnum = '123'
        self.assertEqual(tonumber(stnum), 123)
        stnum = '1.23'
        self.assertEqual(tonumber(stnum), 1.23)
        stnum = '12D3'
        self.assertEqual(tonumber(stnum), 12000)
        stnum = '12E-3'
        self.assertEqual(tonumber(stnum), 0.012)

    def test_nums_and_text(self):
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
