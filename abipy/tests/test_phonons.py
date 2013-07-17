"""Tests for phonons"""
from __future__ import print_function, division

from abipy.phonons import PhononBands
from abipy.tests import *


class PhononBandsTest(AbipyTest):

    def test_base(self):
        """Base tests for PhononBands"""
        filename = get_reference_file("trf2_5.out_PHBST.nc")

        phbands = PhononBands.from_file(filename)
        print(phbands)

        phbands.show_qpoints()

        self.assertEqual(phbands.minfreq, 0.0)
        #self.assertEqual(phbands.maxfreq, 30)

        #dos = phbands.get_dos()
        #phbands.plot()
        #phbands.plot_with_dos(dos)
        #mu = dos.find_mu(8, atol=1.e-4)
        #self.assert_almost_equal(mu, 6.1469984736625003)
