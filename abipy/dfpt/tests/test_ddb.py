"""Tests for phonons"""
from __future__ import print_function, division

import os
import numpy as np

from abipy.core.testing import *
from abipy.dfpt.ddb import DdbFile


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')

class DdbTest(AbipyTest):

    def test_ddb_methods(self):
        ddb_fname = os.path.join(test_dir, "AlAs_1qpt_DDB")
        with DdbFile(ddb_fname) as ddb:
            # Test qpoints.
            self.assert_equal(ddb.qpoints, np.reshape([0.25, 0, 0], (-1,3)))

            # Test header
            h = ddb.header
            #print(h)
            assert h.version == 100401
            self.assert_equal(h.ecut, 3)
            self.assert_equal(h.occ, 4*[2])
            assert h.xred.shape == (h.natom, 3) and h.kpt.shape == (h.nkpt, 3)
            print(h.znucl)

            self.assert_equal(h.symrel[1].T.ravel(), [0, -1, 1, 0, -1, 0, 1, -1, 0])
            self.assert_equal(h.symrel[2].T.ravel(), [-1, 0, 0, -1, 0, 1, -1, 1, 0])

            # Test structure
            struct = ddb.structure
            print(struct)
            assert struct.formula == "Al1 As1"

            #ncfile = ddb.get_phmodes_at_qpoint()
            #print(ncfile)
            #ncfile = ddb.get_phbands_and_dos(ngqpt=(4,4,4))
            #print(ncfile)


if __name__ == "__main__": 
    import unittest
    unittest.main()
