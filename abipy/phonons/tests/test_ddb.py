"""Tests for phonons"""
from __future__ import print_function, division

import os
import numpy as np

from abipy.phonons.ddb import DdbFile
from abipy.core.testing import *


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

            # Test structure
            struct = ddb.structure
            print(struct)
            assert struct.formula == "Al1 As1"


if __name__ == "__main__": 
    import unittest
    unittest.main()
