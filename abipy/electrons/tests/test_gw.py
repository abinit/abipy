"""Tests for electrons.gw module"""
from __future__ import print_function, division

import collections

from abipy import abiopen
from abipy.electrons.gw import *
from abipy.tests import *

class TestQPList(AbipyTest):

    def setUp(self):
        self.sigres = sigres= abiopen(get_reference_file("tgw1_9o_DS4_SIGRES.nc"))
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

        # Test QP object.
        qp = qplist[0]
        print(qp)
        print(qp.tips)

        self.assertAlmostEqual(qp.e0, -5.04619941555265)
        self.assertAlmostEqual(qp.qpe.real, -4.76022137474714)
        self.assertAlmostEqual(qp.qpe.imag, -0.011501666037697)
        self.assertAlmostEqual(qp.sigxme, -16.549383605401)


#class TestSigresFile(AbipyTest):
#    def SetUp(self):
#        self.sigres = abiopen(get_reference_file("tgw1_9o_DS4_SIGRES.nc"))
#
#    def test_readwrite_GWFile(self):
#        """Test read & write GWFile"""
#        pass


if __name__ == "__main__":
    import unittest
    unittest.main()
