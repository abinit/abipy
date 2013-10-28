"""Tests for phonons"""
from __future__ import print_function, division

import abipy.data as data

from abipy.phonons import PhononBands
from abipy.core.testing import *


class PhononBandsTest(AbipyTest):

    def test_base(self):
        """Base tests for PhononBands"""
        filename = data.ref_file("trf2_5.out_PHBST.nc")

        phbands = PhononBands.from_file(filename)
        print(phbands)

        self.serialize_with_pickle(phbands, protocols=[-1], test_eq=False)

        self.assertEqual(phbands.minfreq, 0.0)
        #self.assertEqual(phbands.maxfreq, 30)

        #dos = phbands.get_dos()

if __name__ == "__main__":
    import unittest
    unittest.main()
