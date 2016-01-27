"""Tests for phonons"""
from __future__ import print_function, division

import tempfile
import abipy.data as abidata

from abipy.dfpt.phonons import PhononBands, PhononDos
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

        # Test XYZ vib
        _, filename = tempfile.mkstemp(text=True)
        phbands.create_xyz_vib(iqpt=0, filename=filename, max_supercell=[4,4,4])

        #dos = phbands.get_phdos()
        #print(dos)


class PhononDosTest(AbipyTest):

    def test_api(self):
        """Testing PhononDos API with fake data."""
        dos = PhononDos(mesh=[1,2,3], values=[4,5,6])
        assert dos.mesh.tolist() == [1,2,3] and dos.h == 1 and dos.values.tolist() == [4,5,6]
        print(dos)
        dos.idos
        h = dos.get_harmonic_thermo(1, 10)
        assert h is not None

        if self.has_matplotlib():
            dos.plot(show=False)


if __name__ == "__main__":
    import unittest
    unittest.main()
