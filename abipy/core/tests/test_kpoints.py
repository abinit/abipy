#!/usr/bin/env python
"""Tests for kpoints.kpoints module."""
from __future__ import print_function, division

import itertools
import unittest
import numpy as np
import abipy.data as abidata

from pymatgen.core.lattice import Lattice
from abipy import abilab
from abipy.core.kpoints import (wrap_to_ws, wrap_to_bz, issamek, Kpoint, KpointList, IrredZone, Kpath, KpointsReader,
    has_timrev_from_kptopt, KSamplingInfo, as_kpoints, rc_list, kmesh_from_mpdivs, map_grid2ibz,
    set_atol_kdiff, set_spglib_tols, kpath_from_bounds_and_ndivsm)  #Ktables,
from abipy.core.testing import AbipyTest


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


class TestHelperFunctions(AbipyTest):

    def test_wrap_to_bz(self):
        """Testing wrap_to_bz"""
        self.assertAlmostEqual(wrap_to_bz( 0.0), 0.0)
        self.assertAlmostEqual(wrap_to_bz( 1.0), 0.0)
        self.assertAlmostEqual(wrap_to_bz( 0.2), 0.2)
        self.assertAlmostEqual(wrap_to_bz(-0.2), 0.8)
        self.assertAlmostEqual(wrap_to_bz( 3.2), 0.2)
        self.assertAlmostEqual(wrap_to_bz(-3.2), 0.8)

    def test_is_diagonal(self):
        """Testing is_diagonal"""
        from abipy.core.kpoints import is_diagonal
        assert is_diagonal(np.eye(3, dtype=np.int))
        assert is_diagonal(np.eye(3, dtype=np.float))
        a = np.eye(3, dtype=np.float)
        atol = 1e-12
        a[1, 2] = atol
        assert is_diagonal(a, atol=atol)
        assert not is_diagonal(a, atol=atol / 2)

    def test_has_timrev_from_kptopt(self):
        """Testing has_timrev_from_kptopt."""
        assert has_timrev_from_kptopt(1)
        assert not has_timrev_from_kptopt(4)
        assert has_timrev_from_kptopt(-7)

    def test_kptopt2str(self):
        """Testing kptopt2str."""
        from abipy.core.kpoints import kptopt2str
        for kptopt in [-5, 0, 1, 2, 3, 4]:
            assert kptopt2str(kptopt, verbose=1 if kptopt != 1 else 0)

    def test_kpath_from_bounds_and_ndivsm(self):
        """Testing kpath_from_bounds_and_ndivsm."""
        structure = abilab.Structure.as_structure(abidata.cif_file("si.cif"))
        with self.assertRaises(ValueError):
            kpath_from_bounds_and_ndivsm([(0, 0, 0)], 5, structure)
        with self.assertRaises(ValueError):
            kpath_from_bounds_and_ndivsm([(0, 0, 0), (0, 0, 0)], 5, structure)

        path = kpath_from_bounds_and_ndivsm([(0, 0, 0), (0.5, 0, 0)], 5, structure)
        self.assert_equal(path, [[0.0, 0.0, 0.0 ],
                                 [0.1, 0.0, 0.0 ],
                                 [0.2, 0.0, 0.0 ],
                                 [0.3, 0.0, 0.0 ],
                                 [0.4, 0.0, 0.0 ],
                                 [0.5, 0.0, 0.0 ]])

class TestKpoint(AbipyTest):
    """Unit tests for Kpoint object."""

    def setUp(self):
        self.lattice = Lattice([0.5,0.5,0,0,0.5,0,0,0,0.4])

        # Test API to set tolerances.

        # Change _ATOL_KDIFF
        atol_default = 1e-8
        assert set_atol_kdiff(1e-6) == atol_default
        set_atol_kdiff(atol_default)

        # Change spglib tolerances.
        symprec_default, angle_tolerance_default = 1e-5, -1.0
        s, a = set_spglib_tols(1e-4, -2.0)
        assert s == symprec_default
        assert a == angle_tolerance_default
        set_spglib_tols(symprec_default, angle_tolerance_default)

    def test_kpoint_algebra(self):
        """Test k-point algebra."""
        lattice = self.lattice
        gamma = Kpoint([0, 0, 0], lattice)
        pgamma = Kpoint([1, 0, 1], lattice)
        X = Kpoint([0.5, 0, 0], lattice)
        K = Kpoint([1/3, 1/3, 1/3], lattice)
        repr(X); str(X)
        assert X.to_string(verbose=2)

        assert gamma.is_gamma()
        assert not pgamma.is_gamma()
        assert pgamma.is_gamma(allow_umklapp=True)
        assert not X.is_gamma()

        # TODO
        #assert np.all(np.array(X) == X.frac_coords)

        self.serialize_with_pickle(X, protocols=[-1])
        self.assert_almost_equal(X.versor().norm, 1.0)

        other_gamma = Kpoint.gamma(lattice, weight=1)
        assert other_gamma == [0, 0, 0]
        assert other_gamma == gamma
        assert gamma.versor() == gamma

        X_outside = Kpoint([1.5, 0, 0], lattice)
        assert X_outside.wrap_to_ws() == X
        assert X_outside.wrap_to_ws() == [0.5, 0, 0]

        X_outside = Kpoint([0.7, 0, 0], lattice)
        assert X_outside.wrap_to_bz() == [-0.3, 0, 0]

        assert X[0] == 0.5
        self.assert_equal(pgamma[:2].tolist(), [1,0])

        assert gamma == pgamma
        assert gamma + pgamma == gamma
        assert pgamma + X == X
        assert gamma != X
        # TODO
        #assert gamma != 0

        assert X.norm == (gamma + X).norm
        assert X.norm ==  (gamma + X).norm
        assert X.norm == np.sqrt(np.sum(X.cart_coords**2))
        # TODO
        #assert X != 0.5

        assert hash(gamma) == hash(pgamma)
        if hash(K) != hash(X):
            assert K != X

        # test on_border
        assert not gamma.on_border
        assert X.on_border
        assert not K.on_border


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
        assert kpts is newkpts

        kpts = as_kpoints([1, 2, 3, 4, 5, 6], lattice)
        assert len(kpts) == 2
        assert kpts[0] == Kpoint([1, 2, 3], lattice)
        assert kpts[1] == Kpoint([4, 5, 6], lattice)

    def test_kpointlist(self):
        """Test KpointList."""
        lattice = self.lattice

        frac_coords = [0, 0, 0, 1/2, 1/2, 1/2, 1/3, 1/3, 1/3]
        weights = [0.1, 0.2, 0.7]

        klist = KpointList(lattice, frac_coords, weights=weights)
        repr(klist); str(klist)

        self.serialize_with_pickle(klist, protocols=[-1])
        self.assertMSONable(klist, test_if_subclass=False)

        self.assert_equal(klist.frac_coords.flatten(), frac_coords)
        self.assert_equal(klist.get_cart_coords(), np.reshape([k.cart_coords for k in klist], (-1, 3)))
        assert klist.sum_weights() == 1
        assert len(klist) == 3

        for i, kpoint in enumerate(klist):
            assert kpoint in klist
            assert klist.count(kpoint) == 1
            assert klist.find(kpoint) == i

        # Changing the weight of the Kpoint object should change the weights of klist.
        for kpoint in klist: kpoint.set_weight(1.0)
        assert np.all(klist.weights == 1.0)

        # Test find_closest
        iclose, kclose, dist = klist.find_closest([0, 0, 0])
        assert iclose == 0 and dist == 0.

        iclose, kclose, dist = klist.find_closest(Kpoint([0.001, 0.002, 0.003], klist.reciprocal_lattice))
        assert iclose == 0
        self.assert_almost_equal(dist, 0.001984943324127921)

        # Compute mapping k_index --> (k + q)_index, g0
        k2kqg = klist.get_k2kqg_map((0, 0, 0))
        assert all(ikq == ik for ik, (ikq, g0) in k2kqg.items())
        k2kqg = klist.get_k2kqg_map((1/2, 1/2, 1/2))
        assert len(k2kqg) == 2
        assert k2kqg[0][0] == 1 and np.all(k2kqg[0][1] == 0)
        assert k2kqg[1][0] == 0 and np.all(k2kqg[1][1] == 1)

        frac_coords = [0, 0, 0, 1/2, 1/3, 1/3]
        other_klist = KpointList(lattice, frac_coords)

        # Test __add__
        add_klist = klist + other_klist

        for k in itertools.chain(klist, other_klist):
            assert k in add_klist

        assert add_klist.count([0,0,0]) == 2

        # Remove duplicated k-points.
        add_klist = add_klist.remove_duplicated()
        assert add_klist.count([0,0,0]) == 1
        assert len(add_klist) == 4
        assert add_klist == add_klist.remove_duplicated()


class TestIrredZone(AbipyTest):

    def test_irredzone_api(self):
        """Testing IrredZone API."""
        structure = abilab.Structure.as_structure(abidata.cif_file("si.cif"))

        ibz = IrredZone.from_ngkpt(structure, ngkpt=[4, 4, 4], shiftk=[0.0, 0.0, 0.0], verbose=2)
        repr(ibz); str(ibz)
        assert ibz.to_string(verbose=2)
        assert ibz.is_ibz
        assert len(ibz) == 8
        assert ibz.ksampling.kptopt == 1
        self.assert_equal(ibz.ksampling.mpdivs, [4, 4, 4])

        ibz = IrredZone.from_kppa(structure, kppa=1000, shiftk=[0.5, 0.5, 0.5], kptopt=1, verbose=1)
        assert ibz.is_ibz
        assert len(ibz) == 60
        self.assert_equal(ibz.ksampling.mpdivs, [8, 8, 8])


class TestKpath(AbipyTest):

    def test_kpath_api(self):
        """Testing Kpath API."""
        structure = abilab.Structure.as_structure(abidata.cif_file("si.cif"))

        knames = ["G", "X", "L", "G"]
        kpath = Kpath.from_names(structure, knames, line_density=5)
        repr(kpath); str(kpath)
        assert kpath.to_string(verbose=2, title="Kpath")
        assert not kpath.is_ibz and kpath.is_path
        assert kpath[0].is_gamma and kpath[-1].is_gamma
        #assert len(kpath.ds) == len(self) - 1
        #assert kpath.ksampling.kptopt == 1
        #self.assert_equal(kpath.ksampling.mpdivs, [4, 4, 4])

        assert len(kpath.ds) == len(kpath) - 1
        assert len(kpath.versors) == len(kpath) - 1
        assert len(kpath.lines) == len(knames) - 1
        self.assert_almost_equal(kpath.frac_bounds, structure.get_kcoords_from_names(knames))
        self.assert_almost_equal(kpath.cart_bounds, structure.get_kcoords_from_names(knames, cart_coords=True))

        r = kpath.find_points_along_path(kpath.get_cart_coords())
        assert len(r.ikfound) == len(kpath)
        self.assert_equal(r.ikfound,
            [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  0])

        #kpath = IrredZone.from_kppa(structure, kppa=1000, shiftk=[0.5, 0.5, 0.5], kptopt=1, verbose=1)
        #assert not kpath.is_ibz and kpath.is_path
        #assert len(kpath) == 60
        #self.assert_equal(kpath.ksampling.mpdivs, [8, 8, 8])


class TestKpointsReader(AbipyTest):

    def test_reading(self):
        """Test the reading of Kpoints from netcdf files."""
        filenames = [
            "si_scf_GSR.nc",
            "si_nscf_GSR.nc",
            "si_scf_WFK.nc",
        ]

        for fname in filenames:
            filepath = abidata.ref_file(fname)
            print("About to read file: %s" % filepath)

            with KpointsReader(filepath) as r:
                kpoints = r.read_kpoints()
                repr(kpoints); str(kpoints)

                if "_scf" in fname:
                    # expecting a homogeneous sampling.
                    assert not kpoints.is_path
                    assert kpoints.is_ibz
                    assert kpoints.sum_weights() == 1.0
                    assert kpoints.ksampling.kptopt == 1
                    mpdivs, shifts = kpoints.mpdivs_shifts
                    assert np.all(mpdivs == [8, 8, 8])
                    assert len(shifts) == 1 and np.all(shifts[0] == [0, 0, 0])

                elif "_nscf" in fname:
                    # expecting a path in k-space.
                    assert kpoints.is_path
                    assert not kpoints.is_ibz
                    assert kpoints.ksampling.kptopt == -2
                    mpdivs, shifts = kpoints.mpdivs_shifts
                    assert mpdivs is None
                    assert len(kpoints.lines) == abs(kpoints.ksampling.kptopt)

            # Test pickle and json
            self.serialize_with_pickle(kpoints)
            self.assertMSONable(kpoints, test_if_subclass=False)


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
        def rm_spaces(s):
            return " ".join(s.split()).replace("[ ", "[")

        mpdivs, shifts = [1, 2, 3], [0, 0, 0]

        # No shift, no pbc.
        kmesh = kmesh_from_mpdivs(mpdivs, shifts, order="unit_cell")

        ref_string = \
"""[[ 0.          0.          0.        ]
 [ 0.          0.          0.33333333]
 [ 0.          0.          0.66666667]
 [ 0.          0.5         0.        ]
 [ 0.          0.5         0.33333333]
 [ 0.          0.5         0.66666667]]"""
        self.assertMultiLineEqual(rm_spaces(str(kmesh)), rm_spaces(ref_string))

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
        self.assertMultiLineEqual(rm_spaces(str(pbc_kmesh)), rm_spaces(ref_string))

        # No shift, no pbc, bz order
        bz_kmesh = kmesh_from_mpdivs(mpdivs, shifts, pbc=False, order="bz")

        ref_string = \
"""[[ 0.         -0.5        -0.33333333]
 [ 0.         -0.5         0.        ]
 [ 0.         -0.5         0.33333333]
 [ 0.          0.         -0.33333333]
 [ 0.          0.          0.        ]
 [ 0.          0.          0.33333333]]"""
        self.assertMultiLineEqual(rm_spaces(str(bz_kmesh)), rm_spaces(ref_string))

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
        self.assertMultiLineEqual(rm_spaces(str(bz_kmesh)), rm_spaces(ref_string))


class TestKsamplingInfo(AbipyTest):

    def test_ksampling(self):
        """Test KsamplingInfo API."""
        with self.assertRaises(ValueError):
            KSamplingInfo(foo=1)

        # from_mpdivs constructor
        mpdivs, shifts = [2, 3, 4], [0.5, 0.5, 0.5]
        kptopt = 1
        ksi = KSamplingInfo.from_mpdivs(mpdivs, shifts, kptopt)
        repr(ksi); str(ksi)
        self.assert_equal(ksi.mpdivs, mpdivs)
        self.assert_equal(ksi.kptrlatt, np.diag(mpdivs))
        self.assert_equal(ksi.shifts.flatten(), shifts)
        assert ksi.shifts.shape == (1, 3)
        assert ksi.kptopt == kptopt
        assert ksi.is_mesh
        assert ksi.has_diagonal_kptrlatt
        assert not ksi.is_path

        # from kptrlatt constructor
        kptrlatt = np.diag(mpdivs)
        ksi = KSamplingInfo.from_kptrlatt(kptrlatt, shifts, kptopt)
        repr(ksi); str(ksi)
        assert ksi.kptrlatt.shape == (3, 3)
        self.assert_equal(ksi.kptrlatt, np.diag(mpdivs))
        self.assert_equal(ksi.mpdivs, np.diag(ksi.kptrlatt))
        self.assert_equal(ksi.shifts.flatten(), shifts)
        assert ksi.kptopt == kptopt
        assert ksi.is_mesh
        assert ksi.has_diagonal_kptrlatt
        assert not ksi.is_path

        # kptrlatt with non-zero off-diagonal elements.
        shifts = [0.5, 0.5, 0.5]
        kptrlatt = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        kptopt = 1
        ksi = KSamplingInfo.from_kptrlatt(kptrlatt, shifts, kptopt)
        repr(ksi); str(ksi)
        assert ksi.mpdivs is None
        assert not ksi.has_diagonal_kptrlatt
        assert not ksi.is_path

        # from_kbounds constructor
        kbounds = [0, 0, 0, 1, 1, 1]
        ksi = KSamplingInfo.from_kbounds(kbounds)
        repr(ksi); str(ksi)
        assert (ksi.mpdivs, ksi.kptrlatt, ksi.kptrlatt_orig, ksi.shifts, ksi.shifts_orig) == 5 * (None,)
        assert ksi.kptopt == -1
        assert ksi.kptrlatt is None
        assert not ksi.is_mesh
        assert not ksi.has_diagonal_kptrlatt
        assert ksi.is_path

        assert ksi is KSamplingInfo.as_ksampling(ksi)

        ksi_from_dict = KSamplingInfo.as_ksampling({k: v for k, v in ksi.items()})
        assert ksi_from_dict.kptopt == ksi.kptopt

        ksi_none = KSamplingInfo.as_ksampling(None)
        repr(ksi_none); str(ksi_none)
        assert ksi_none.kptopt == 0
        assert not ksi_none.is_mesh
        assert not ksi_none.is_path


class TestKmappingTools(AbipyTest):

    def setUp(self):
        with abilab.abiopen(abidata.ref_file("mgb2_kmesh181818_FATBANDS.nc")) as ncfile:
            self.mgb2 = ncfile.structure
            assert ncfile.ebands.kpoints.is_ibz
            self.kibz = [k.frac_coords for k in ncfile.ebands.kpoints]
            self.has_timrev = True
            #self.has_timrev = has_timrev_from_kptopt(kptopt)
            self.ngkpt = [18, 18, 18]

    def test_map_grid2ibz(self):
        """Testing map_grid2ibz."""
        bz2ibz = map_grid2ibz(self.mgb2, self.kibz, self.ngkpt, self.has_timrev, pbc=False)

        bz = []
        nx, ny, nz = self.ngkpt
        for ix, iy, iz in itertools.product(range(nx), range(ny), range(nz)):
            bz.append([ix/nz, iy/ny, iz/nz])
        bz = np.reshape(bz, (-1, 3))

        abispg = self.mgb2.abi_spacegroup

        nmax = 54
        errors = []
        for ik_bz, kbz in enumerate(bz[:nmax]):
            ik_ibz = bz2ibz[ik_bz]
            ki = self.kibz[ik_ibz]
            for symmop in abispg.fm_symmops:
                krot = symmop.rotate_k(ki)
                if issamek(krot, kbz):
                    break
            else:
                errors.append((ik_bz, kbz))

        assert not errors

    #def test_with_from_structure_with_symrec(self):
    #    """Generate Ktables from a structure with Abinit symmetries."""
    #    self.mgb2 = self.get_abistructure.mgb2("mgb2_kpath_FATBANDS.nc")
    #    assert self.mgb2.abi_spacegroup is not None
    #    mesh = [4, 4, 4]
    #    k = Ktables(self.mgb2, mesh, is_shift=None, has_timrev=True)
    #    repr(k); str(k)
    #    k.print_bz2ibz()

    #def test_with_structure_without_symrec(self):
    #    """Generate Ktables from a structure without Abinit symmetries."""
    #    assert self.mgb2.abi_spacegroup is None
    #    k = Ktables(self.mgb2, mesh, is_shift, has_timrev)
    #    repr(k); str(k)
    #    k.print_bz2ibz()
