#!/usr/bin/env python
"""Tests for kpoints.kpoints module."""
from __future__ import print_function, division

import itertools
import unittest
import numpy as np
import abipy.data as data

from pymatgen.core.lattice import Lattice
from abipy.core.kpoints import (wrap_to_ws, wrap_to_bz, Kpoint, KpointList, KpointsReader,
                                as_kpoints, rc_list, kmesh_from_mpdivs,)
from abipy.core.testing import *


class TestWrapWS(AbipyTest):

    def test_wrap_to_ws(self):
        """Testing wrap_to_ws"""
        self.assert_almost_equal(wrap_to_ws( 0.5), 0.5)
        self.assert_almost_equal(wrap_to_ws(-0.5), 0.5)
        self.assert_almost_equal(wrap_to_ws( 0.2), 0.2)
        self.assert_almost_equal(wrap_to_ws(-0.3),-0.3)
        self.assert_almost_equal(wrap_to_ws( 0.7),-0.3)
        self.assert_almost_equal(wrap_to_ws( 2.3), 0.3)
        self.assert_almost_equal(wrap_to_ws(-1.2),-0.2)
        self.assert_almost_equal(wrap_to_ws(np.array([0.5,2.3,-1.2])), np.array([0.5,0.3,-0.2]))


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
    """Unit tests for Kpoint object."""

    def setUp(self):
        self.lattice = Lattice([0.5,0.5,0,0,0.5,0,0,0,0.4])

    def test_kpoint_algebra(self):
        """Test k-point algebra."""
        lattice = self.lattice
        gamma = Kpoint([0, 0, 0], lattice)
        pgamma = Kpoint([1, 0, 1], lattice)
        X = Kpoint([0.5, 0, 0], lattice)
        K = Kpoint([1/3, 1/3, 1/3], lattice)
        print(X)

        # TODO
        #assert np.all(np.array(X) == X.frac_coords)

        self.serialize_with_pickle(X, protocols=[-1])
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

        self.assertTrue(hash(gamma) == hash(pgamma))

        if hash(K) != hash(X):
            self.assertTrue(K != X)

        # test on_border
        self.assertFalse(gamma.on_border)
        self.assertTrue(X.on_border)
        self.assertFalse(K.on_border)


class TestKpointList(AbipyTest):
    """Unit tests for KpointList."""

    def setUp(self):
        self.lattice = Lattice([0.5,0.5,0,0,0.5,0,0,0,0.4])

    def test_askpoints(self):
        """Test askpoints."""
        lattice = self.lattice
        kpts = as_kpoints([1, 2, 3], lattice)

        self.serialize_with_pickle(kpts, protocols=[-1])

        newkpts = as_kpoints(kpts, lattice)
        self.assertTrue(kpts is newkpts)

        kpts = as_kpoints([1, 2, 3, 4, 5, 6], lattice)
        self.assertTrue(len(kpts) == 2)
        self.assertTrue(kpts[0] == Kpoint([1, 2, 3], lattice))
        self.assertTrue(kpts[1] == Kpoint([4, 5, 6], lattice))

    def test_kpointlist(self):
        """Test KpointList."""
        lattice = self.lattice

        frac_coords = [0, 0, 0, 1/2, 1/2, 1/2, 1/3, 1/3, 1/3]
        weights = [0.1, 0.2, 0.7]

        klist = KpointList(lattice, frac_coords, weights=weights)

        self.serialize_with_pickle(klist, protocols=[-1])
        self.assertMSONable(klist, test_if_subclass=False)

        self.assertTrue(klist.sum_weights() == 1)
        self.assertTrue(len(klist) == 3)

        for i, kpoint in enumerate(klist):
            self.assertTrue(kpoint in klist)
            self.assertTrue(klist.count(kpoint) == 1)
            self.assertTrue(klist.find(kpoint) == i)

        # Changing the weight of the Kpoint object should change the weights of klist.
        for kpoint in klist: kpoint.set_weight(1.0)
        self.assertTrue(np.all(klist.weights == 1.0))

        frac_coords = [0, 0, 0, 1/2, 1/3, 1/3]

        other_klist = KpointList(lattice, frac_coords)

        # Test __add__
        add_klist = klist + other_klist

        for k in itertools.chain(klist, other_klist):
            self.assertTrue(k in add_klist)

        self.assertTrue(add_klist.count([0,0,0]) == 2)

        # Remove duplicated k-points.
        add_klist = add_klist.remove_duplicated()
        self.assertTrue(add_klist.count([0,0,0]) == 1)
        self.assertTrue(len(add_klist) == 4)
        self.assertTrue(add_klist == add_klist.remove_duplicated())


class TestKpointsReader(AbipyTest):

    @unittest.expectedFailure
    def test_reading(self):
        """Test the reading of Kpoints from netcdf files."""

        filenames = [
            "si_scf_GSR.nc",
            "si_nscf_GSR.nc",
            "si_scf_WFK-etsf.nc",
        ]

        for fname in filenames:
            filepath = data.ref_file(fname)
            print("About to read file: %s" % filepath)

            with KpointsReader(filepath) as r:
                kpoints = r.read_kpoints()

                if "_scf" in fname:
                    # expecting a homogeneous sampling.
                    self.assertTrue(kpoints.is_homogeneous)
                    self.assertTrue(kpoints.sum_weights() == 1.0)

                elif "_nscf" in fname:
                    # expecting a path in k-space.
                    self.assertTrue(kpoints.is_path)


class KmeshTest(AbipyTest):
    def test_rc_list(self):
        """Testing rc_list."""
        # Special case mp=1
        rc = rc_list(mp=1, sh=0.0, pbc=False, order="unit_cell")
        self.assert_equal(rc, [0.0])

        rc = rc_list(mp=1, sh=0.0, pbc=True, order="unit_cell")
        self.assert_equal(rc, [0.0, 1.0])

        rc = rc_list(mp=1, sh=0.0, pbc=False, order="bz")
        self.assert_equal(rc, [0.0])

        rc = rc_list(mp=1, sh=0.0, pbc=True, order="bz")
        self.assert_equal(rc, [0.0, 1.0])

        # Even mp
        rc = rc_list(mp=2, sh=0, pbc=False, order="unit_cell")
        self.assert_equal(rc, [0., 0.5])

        rc = rc_list(mp=2, sh=0, pbc=True, order="unit_cell")
        self.assert_equal(rc, [0., 0.5,  1.])

        rc = rc_list(mp=2, sh=0, pbc=False, order="bz")
        self.assert_equal(rc, [-0.5, 0.0])

        rc = rc_list(mp=2, sh=0, pbc=True, order="bz")
        self.assert_equal(rc, [-0.5,  0.,  0.5])

        rc = rc_list(mp=2, sh=0.5, pbc=False, order="unit_cell")
        self.assert_equal(rc, [0.25, 0.75])

        rc = rc_list(mp=2, sh=0.5, pbc=True, order="unit_cell")
        self.assert_equal(rc, [0.25,  0.75, 1.25])

        rc = rc_list(mp=2, sh=0.5, pbc=False, order="bz")
        self.assert_equal(rc, [-0.25,  0.25])

        rc = rc_list(mp=2, sh=0.5, pbc=True, order="bz")
        self.assert_equal(rc, [-0.25,  0.25,  0.75])

        # Odd mp
        rc = rc_list(mp=3, sh=0, pbc=False, order="unit_cell")
        self.assert_almost_equal(rc, [0.,  0.33333333,  0.66666667])

        rc = rc_list(mp=3, sh=0, pbc=True, order="unit_cell")
        self.assert_almost_equal(rc, [ 0.,  0.33333333,  0.66666667,  1.])

        rc = rc_list(mp=3, sh=0, pbc=False, order="bz")
        self.assert_almost_equal(rc, [-0.33333333,  0.,  0.33333333])

        rc = rc_list(mp=3, sh=0, pbc=True, order="bz")
        self.assert_almost_equal(rc, [-0.33333333,  0.,  0.33333333,  0.66666667])

        rc = rc_list(mp=3, sh=0.5, pbc=False, order="unit_cell")
        self.assert_almost_equal(rc, [ 0.16666667, 0.5, 0.83333333])

        rc = rc_list(mp=3, sh=0.5, pbc=True, order="unit_cell")
        self.assert_almost_equal(rc, [ 0.16666667, 0.5,  0.83333333, 1.16666667])

        rc = rc_list(mp=3, sh=0.5, pbc=False, order="bz")
        self.assert_almost_equal(rc, [-0.5, -0.16666667,  0.16666667])

        rc = rc_list(mp=3, sh=0.5, pbc=True, order="bz")
        self.assert_almost_equal(rc, [-0.5, -0.16666667,  0.16666667,  0.5])

    def test_unshifted_kmesh(self):
        """Testing the generation of unshifted kmeshes."""
        mpdivs, shifts = [1,2,3], [0,0,0]

        # No shift, no pbc.
        kmesh = kmesh_from_mpdivs(mpdivs, shifts, order="unit_cell")

        ref_string = \
"""[[ 0.          0.          0.        ]
 [ 0.          0.          0.33333333]
 [ 0.          0.          0.66666667]
 [ 0.          0.5         0.        ]
 [ 0.          0.5         0.33333333]
 [ 0.          0.5         0.66666667]]"""
        self.assertMultiLineEqual(str(kmesh), ref_string)

        # No shift, with pbc.
        pbc_kmesh = kmesh_from_mpdivs(mpdivs, shifts, pbc=True, order="unit_cell")

        ref_string = \
"""[[ 0.          0.          0.        ]
 [ 0.          0.          0.33333333]
 [ 0.          0.          0.66666667]
 [ 0.          0.          1.        ]
 [ 0.          0.5         0.        ]
 [ 0.          0.5         0.33333333]
 [ 0.          0.5         0.66666667]
 [ 0.          0.5         1.        ]
 [ 0.          1.          0.        ]
 [ 0.          1.          0.33333333]
 [ 0.          1.          0.66666667]
 [ 0.          1.          1.        ]
 [ 1.          0.          0.        ]
 [ 1.          0.          0.33333333]
 [ 1.          0.          0.66666667]
 [ 1.          0.          1.        ]
 [ 1.          0.5         0.        ]
 [ 1.          0.5         0.33333333]
 [ 1.          0.5         0.66666667]
 [ 1.          0.5         1.        ]
 [ 1.          1.          0.        ]
 [ 1.          1.          0.33333333]
 [ 1.          1.          0.66666667]
 [ 1.          1.          1.        ]]"""
        self.assertMultiLineEqual(str(pbc_kmesh), ref_string)

        # No shift, no pbc, bz order
        bz_kmesh = kmesh_from_mpdivs(mpdivs, shifts, pbc=False, order="bz")

        ref_string = \
"""[[ 0.         -0.5        -0.33333333]
 [ 0.         -0.5         0.        ]
 [ 0.         -0.5         0.33333333]
 [ 0.          0.         -0.33333333]
 [ 0.          0.          0.        ]
 [ 0.          0.          0.33333333]]"""
        self.assertMultiLineEqual(str(bz_kmesh), ref_string)

        # No shift, pbc, bz order
        bz_kmesh = kmesh_from_mpdivs(mpdivs, shifts, pbc=True, order="bz")

        ref_string = \
"""[[ 0.         -0.5        -0.33333333]
 [ 0.         -0.5         0.        ]
 [ 0.         -0.5         0.33333333]
 [ 0.         -0.5         0.66666667]
 [ 0.          0.         -0.33333333]
 [ 0.          0.          0.        ]
 [ 0.          0.          0.33333333]
 [ 0.          0.          0.66666667]
 [ 0.          0.5        -0.33333333]
 [ 0.          0.5         0.        ]
 [ 0.          0.5         0.33333333]
 [ 0.          0.5         0.66666667]
 [ 1.         -0.5        -0.33333333]
 [ 1.         -0.5         0.        ]
 [ 1.         -0.5         0.33333333]
 [ 1.         -0.5         0.66666667]
 [ 1.          0.         -0.33333333]
 [ 1.          0.          0.        ]
 [ 1.          0.          0.33333333]
 [ 1.          0.          0.66666667]
 [ 1.          0.5        -0.33333333]
 [ 1.          0.5         0.        ]
 [ 1.          0.5         0.33333333]
 [ 1.          0.5         0.66666667]]"""
        self.assertMultiLineEqual(str(bz_kmesh), ref_string)
