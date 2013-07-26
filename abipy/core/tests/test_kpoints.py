#!/usr/bin/env python
"""Tests for kpoints.kpoints module."""

import numpy as np

from pymatgen.core.lattice import Lattice
from abipy.core.kpoints import *
from abipy.tests import *

class TestWrapWS(AbipyTest):

    def test_wrap_to_ws(self):
        """Testing wrap_to_ws"""
        self.assertAlmostEqual(wrap_to_ws( 0.5), 0.5)
        self.assertAlmostEqual(wrap_to_ws(-0.5), 0.5)
        self.assertAlmostEqual(wrap_to_ws( 0.2), 0.2)
        self.assertAlmostEqual(wrap_to_ws(-0.3),-0.3)
        self.assertAlmostEqual(wrap_to_ws( 0.7),-0.3)
        self.assertAlmostEqual(wrap_to_ws( 2.3), 0.3)
        self.assertAlmostEqual(wrap_to_ws(-1.2),-0.2)
        #self.assertAlmostEqual( wrap_to_ws([0.5,2.3,-1.2]), np.array([0.5,0.3,-0.2]) )


class TestWrapBZ(AbipyTest):

    def test_wrap_to_bz(self):
        """Testing wrap_to_bz"""
        self.assertAlmostEqual(wrap_to_bz( 0.0), 0.0)
        self.assertAlmostEqual(wrap_to_bz( 1.0), 0.0)
        self.assertAlmostEqual(wrap_to_bz( 0.2), 0.2)
        self.assertAlmostEqual(wrap_to_bz(-0.2), 0.8)
        self.assertAlmostEqual(wrap_to_bz( 3.2), 0.2)
        self.assertAlmostEqual(wrap_to_bz(-3.2), 0.8)

class TestKpoint(AbipyTest):
    """Unit tests for Kpoints."""

    def setUp(self):
        self.lattice = Lattice([0.5,0.5,0,0,0.5,0,0,0,0.4])

    def test_askpoints(self):
        """Test askpoints."""
        lattice = self.lattice
        kpts = askpoints([1,2,3], lattice)

        newkpts = askpoints(kpts, lattice)
        self.assertTrue(kpts is newkpts)

        kpts = askpoints([1,2,3,4,5,6], lattice)
        self.assertTrue(len(kpts) == 2)
        self.assertTrue(kpts[0] == Kpoint([1,2,3], lattice))
        self.assertTrue(kpts[1] == Kpoint([4,5,6], lattice))

    def test_kpoint_algebra(self):
        """Test k-point algebra."""
        lattice = self.lattice
        gamma = Kpoint([0,0,0], lattice)
        pgamma = Kpoint([1,0,1], lattice)
        X = Kpoint([0.5,0,0], lattice)
        print(X)

        self.assert_almost_equal(X.versor().norm, 1.0)

        self.assertTrue(X[0] == 0.5)
        self.assertListEqual(pgamma[:2].tolist(), [1,0])

        self.assertEqual(gamma, pgamma)
        self.assertEqual(gamma + pgamma, gamma)
        self.assertEqual(pgamma + X, X)
        self.assertNotEqual(gamma, X)

        self.assertEqual(X.norm, (gamma + X).norm)

        self.assertEqual(X.norm, (gamma + X).norm)
        self.assertEqual(X.norm, np.sqrt(np.sum(X.cart_coords**2)))


if __name__ == "__main__":
    import unittest
    unittest.main()
