"""Tests for electrons.gw module"""
from __future__ import print_function, division

import collections
import numpy as np
import abipy.data as data

from abipy.abilab import abiopen
from abipy.electrons.gw import *
from abipy.electrons.gw import SIGRES_Reader
from abipy.core.testing import *

class TestQPList(AbipyTest):

    def setUp(self):
        self.sigres = sigres= abiopen(data.ref_file("tgw1_9o_DS4_SIGRES.nc"))
        self.qplist = sigres.get_qplist(spin=0, kpoint=sigres.gwkpoints[0])

    def test_qplist(self):
        """Test QPList object."""
        qplist = self.qplist
        self.assertTrue(isinstance(qplist, collections.Iterable))
        
        print(qplist)
        qplist_copy = qplist.copy()
        self.assertTrue(qplist_copy == qplist)

        qpl_e0sort = qplist.sort_by_e0()
        qpl_e0sort.get_e0mesh()

        with self.assertRaises(ValueError):
            qplist.get_e0mesh()

        with self.assertRaises(ValueError):
            qplist.merge(qpl_e0sort)

        other_qplist = self.sigres.get_qplist(spin=0, kpoint=self.sigres.gwkpoints[1])

        qpl_merge = qplist.merge(other_qplist)

        for qp in qplist:
            self.assertTrue(qp in qpl_merge)

        for qp in other_qplist:
            self.assertTrue(qp in qpl_merge)

        # Test QPState object.
        qp = qplist[0]
        print(qp)
        print(qp.tips)

        self.assertAlmostEqual(qp.e0, -5.04619941555265)
        self.assertAlmostEqual(qp.qpe.real, -4.76022137474714)
        self.assertAlmostEqual(qp.qpe.imag, -0.011501666037697)
        self.assertAlmostEqual(qp.sigxme, -16.549383605401)


#class TestSigresReader(AbipyTest):
#    def test_base(self):
#        """Test SIGRES Reader."""
#        with SIGRES_Reader(data.ref_file("tgw1_9o_DS4_SIGRES.nc")) as r:
#            #params = r.read_params()


class TestSigresFile(AbipyTest):

    def test_readall(self):
        for path in data.SIGRES_NCFILES:
            sigres = abiopen(path)

    def test_base(self):
        """Test SIGRES File."""
        sigres = abiopen(data.ref_file("tgw1_9o_DS4_SIGRES.nc"))

        self.assertTrue(sigres.nsppol == 1)

        # Markers are initialied in __init__
        self.assertTrue(sigres.ebands.markers)

        # In this run IBZ = kptgw
        self.assertTrue(len(sigres.ibz) == 6)
        self.assertTrue(sigres.gwkpoints == sigres.ibz)

        kptgw_coords = np.reshape([
            -0.25, -0.25, 0,
            -0.25, 0.25, 0,
            0.5, 0.5, 0,
            -0.25, 0.5, 0.25,
            0.5, 0, 0,
            0, 0, 0 
        ], (-1,3))

        self.assert_almost_equal(sigres.ibz.frac_coords, kptgw_coords)

        qpgaps = [3.53719151871085, 4.35685250045637, 4.11717896881632, 
                  8.71122659251508, 3.29693118466282, 3.125545059031]

        self.assert_almost_equal(sigres.qpgaps, np.reshape(qpgaps, (1,6)))


if __name__ == "__main__":
    import unittest
    unittest.main()
