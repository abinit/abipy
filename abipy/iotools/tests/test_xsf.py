#!/usr/bin/env python
"""Tests for xsf module"""
from __future__ import print_function, division

import tempfile
import numpy as np
import abipy.data as data 

from abipy.core.testing import *
from abipy.iotools.xsf import *

class TestXsfUtils(AbipyTest):
    """Unit tests for Density."""

    def setUp(self):
        self.mgb2 = data.structure_from_ucell("MgB2")

    def test_xsf_write_structure(self):
        """Testing crystalline structures in the XSF format."""
        tmp_file = tempfile.TemporaryFile(mode="w+")

        xsf_write_structure(tmp_file, self.mgb2)

        xsf_string = \
"""CRYSTAL
# Primitive lattice vectors in Angstrom
PRIMVEC 1
 2.67255439607878 1.54300000000000 0.00000000000000
 -2.67255439607878 1.54300000000000 0.00000000000000
 0.00000000000000 0.00000000000000 3.52300000000000
# Cartesian coordinates in Angstrom.
PRIMCOORD 1
 3 1
 12     0.00000000000000     0.00000000000000     0.00000000000000
  5    -0.89085146535959     1.54300000000000     1.76150000000000
  5     0.89085146535959     1.54300000000000     1.76150000000000
"""
        tmp_file.seek(0)
        self.assertMultiLineEqual(tmp_file.read(), xsf_string)

    def test_xsf_write_data(self):
        """Testing XSF file with datasets."""
        # 2 x 3 x 2 grid without pbc stored in fortran mode.
        data = np.reshape(np.arange(12), (2,3,2)).T
        tmp_file = tempfile.TemporaryFile(mode="w+")
        xsf_write_data(tmp_file, self.mgb2, data, add_replicas=True)

        xsf_string = \
"""BEGIN_BLOCK_DATAGRID_3D
 data
 BEGIN_DATAGRID_3Dgrid#1
3 4 3
0.000000 0.000000 0.000000
2.672554 1.543000 0.000000
-2.672554 1.543000 0.000000
0.000000 0.000000 3.523000
0.000000 1.000000 0.000000
2.000000 3.000000 2.000000
4.000000 5.000000 4.000000
0.000000 1.000000 0.000000

6.000000 7.000000 6.000000
8.000000 9.000000 8.000000
10.000000 11.000000 10.000000
6.000000 7.000000 6.000000

0.000000 1.000000 0.000000
2.000000 3.000000 2.000000
4.000000 5.000000 4.000000
0.000000 1.000000 0.000000

 END_DATAGRID_3D
END_BLOCK_DATAGRID_3D
"""
        tmp_file.seek(0)
        self.assertMultiLineEqual(tmp_file.read(), xsf_string)

        # Complex array will raise TypeError since we should specify the type. 
        cplx_data = np.array(data, dtype=np.complex)

        with self.assertRaises(TypeError):
            xsf_write_data(tmp_file, self.mgb2, cplx_data)

        tmp_file.seek(0)
        xsf_write_data(tmp_file, self.mgb2, cplx_data, cplx_mode="re")
        tmp_file.seek(0)
        self.assertMultiLineEqual(tmp_file.read(), xsf_string)

    #def test_bxsf_write(self):
    #    bxsf_write(file, bands3d, structure, fermie):


if __name__ == "__main__":
    import unittest
    unittest.main()
