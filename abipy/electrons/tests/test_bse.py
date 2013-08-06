"""Tests for electrons.bse module"""
import unittest

from abipy.electrons.bse import *
from abipy.tests import *


class TestMDF_Reader(AbipyTest):

    def test_MDF_reading(self):
        """Test MDF_Reader."""
        mdf_file = get_reference_file("tbs_4o_DS2_MDF.nc")

        with MDF_Reader(mdf_file) as r:
            self.assertTrue(len(r.wmesh) == r.read_dimvalue("number_of_frequencies"))
            self.assertTrue(len(r.qpoints) == r.read_dimvalue("number_of_qpoints"))

            exc_mdf = r.read_exc_mdf()
            exc_mdf.raw_print()

            rpanlf_mdf = r.read_rpanlf_mdf()
            gwnlf_mdf = r.read_gwnlf_mdf()

            #exc_mdf.plot()
            plotter = MDF_Plotter()

            plotter.add_mdf("EXC", exc_mdf)
            plotter.add_mdf("KS-RPA", rpanlf_mdf)
            plotter.add_mdf("GW-RPA", gwnlf_mdf)
            #plotter.plot()

if __name__ == "__main__":
    unittest.main()
