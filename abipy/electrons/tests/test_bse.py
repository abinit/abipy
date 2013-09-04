"""Tests for electrons.bse module"""
from __future__ import print_function, division

import abipy.data as data

from abipy.electrons.bse import *
from abipy.core.testing import *

class TestMDF_Reader(AbipyTest):

    def test_MDF_reading(self):
        """Test MDF_Reader."""
        mdf_file = data.ref_file("tbs_4o_DS2_MDF.nc")

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

    def test_TSR(self):
        """Test the computation of Tensor"""
        mdf_file = MDF_File(data.ref_file("tbs_4o_DS2_MDF.nc"))
 
        exc_tsr = mdf_file.get_tensor("exc")
        rpa_tsr = mdf_file.get_tensor("rpa")
        gw_tsr = mdf_file.get_tensor("gwrpa")

        #exc_tsr.plot()

if __name__ == "__main__":
    import unittest
    unittest.main()
