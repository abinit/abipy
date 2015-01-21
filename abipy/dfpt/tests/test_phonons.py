"""Tests for phonons"""
from __future__ import print_function, division

import tempfile
import abipy.data as abidata

from abipy.dfpt.phonons import PhononBands
from abipy.core.testing import *


class PhononBandsTest(AbipyTest):

    def test_base(self):
        """Base tests for PhononBands"""
        filename = abidata.ref_file("trf2_5.out_PHBST.nc")

        phbands = PhononBands.from_file(filename)
        print(phbands)

        self.serialize_with_pickle(phbands, protocols=[-1], test_eq=False)

        self.assertEqual(phbands.minfreq, 0.0)
        #self.assertEqual(phbands.maxfreq, 30)

        #dos = phbands.get_dos()

        # Test XYZ vib
        _, filename = tempfile.mkstemp(text=True)
        phbands.create_xyz_vib(iqpt=0, filename=filename, max_supercell=[4,4,4])


if __name__ == "__main__":
    import unittest
    unittest.main()
